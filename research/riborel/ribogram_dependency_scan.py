#!/usr/bin/env python3
"""RiboGRAM substitution-based nucleotide-dependency maps.

This implements the saturation-mutagenesis score introduced by Tomaz da
Silva et al. (Nature Genetics 2025), adapted to RiboGRAM's A/C/G/U masked
language-model head.  For every source position, all three alternative bases
are substituted without masking.  For every target position and target base,
the absolute change in log2 odds relative to the reference sequence is
calculated.  The directed source-to-target dependency is the maximum over the
three substitutions and four target bases.

The program deliberately writes a dense directed matrix for validation and
representative-sequence scans.  A later atlas stage should sparsify maps using
distance- and class-matched null thresholds rather than one global cutoff.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
from contextlib import nullcontext
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import torch


BASES = "ACGU"


@dataclass
class ScanStats:
    sequence_id: str
    length: int
    forward_sequences: int
    mutation_batches: int
    max_batch_size_used: int
    elapsed_seconds: float
    sequences_per_second: float
    nucleotide_predictions_per_second: float
    dependency_min: float
    dependency_median: float
    dependency_p99: float
    dependency_max: float
    mean_source_influence: float
    mean_target_sensitivity: float
    device: str
    gpu_name: str
    torch_version: str
    score_definition: str


def normalize_rna(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper().replace("T", "U")
    if not normalized:
        raise ValueError("sequence is empty")
    invalid = sorted(set(normalized).difference(BASES))
    if invalid:
        raise ValueError(f"sequence contains unsupported symbols: {invalid}")
    return normalized


def read_single_fasta(path: Path) -> Tuple[str, str]:
    identifier: Optional[str] = None
    pieces: List[str] = []
    with path.open("rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    raise ValueError("FASTA must contain exactly one record")
                identifier = line[1:].split()[0]
            else:
                if identifier is None:
                    raise ValueError("FASTA sequence appears before its header")
                pieces.append(line)
    if identifier is None:
        raise ValueError("FASTA contains no record")
    return identifier, normalize_rna("".join(pieces))


def biological_logits(output: Dict[str, torch.Tensor], length: int) -> torch.Tensor:
    """Return A/C/G/U logits at biological positions as float32."""
    # RiboGRAM base vocabulary: <cls>, <pad>, <eos>, <unk>, A, C, G, U, <mask>
    return output["logits"][:, 1 : length + 1, 4:8].float()


def log2_odds_from_logits(logits: torch.Tensor, eps: float) -> torch.Tensor:
    # Renormalizing the four biological symbols is equivalent to applying a
    # softmax directly to their logits.  This excludes special-token mass.
    probabilities = torch.softmax(logits, dim=-1)
    probabilities = probabilities.clamp(min=eps, max=1.0 - eps)
    return torch.log2(probabilities) - torch.log2(1.0 - probabilities)


def choose_variant_batch_size(
    length_with_specials: int,
    max_batch_size: int,
    attention_element_budget: int,
) -> int:
    """Choose a multiple-of-three batch that respects approximate B*L^2 cost."""
    by_attention = max(3, attention_element_budget // (length_with_specials**2))
    chosen = max(3, min(max_batch_size, by_attention))
    return max(3, (chosen // 3) * 3)


def compute_dependency_map(
    model: torch.nn.Module,
    alphabet: object,
    converter: object,
    sequence_id: str,
    sequence: str,
    *,
    device: torch.device,
    max_batch_size: int,
    attention_element_budget: int,
    eps: float,
) -> Tuple[np.ndarray, np.ndarray, ScanStats]:
    sequence = normalize_rna(sequence)
    length = len(sequence)
    _, _, reference_tokens = converter([(sequence_id, sequence)])
    reference_tokens = reference_tokens.to(device)
    length_with_specials = int(reference_tokens.shape[1])
    variant_batch_size = choose_variant_batch_size(
        length_with_specials,
        max_batch_size=max_batch_size,
        attention_element_budget=attention_element_budget,
    )
    sources_per_batch = max(1, variant_batch_size // 3)

    base_token_ids = torch.tensor(
        [alphabet.get_idx(base) for base in BASES],
        dtype=torch.long,
        device=device,
    )
    reference_base_ids = reference_tokens[0, 1 : length + 1]

    dependency = torch.zeros((length, length), dtype=torch.float32, device="cpu")
    variant_influence = torch.zeros((length, 3), dtype=torch.float32, device="cpu")
    mutation_batches = 0
    max_batch_size_used = 0
    started = time.perf_counter()

    def autocast_context():
        return (
            torch.autocast(device_type="cuda", dtype=torch.float16)
            if device.type == "cuda"
            else nullcontext()
        )

    with torch.inference_mode():
        with autocast_context():
            reference_output = model(reference_tokens)
        reference_log_odds = log2_odds_from_logits(
            biological_logits(reference_output, length), eps=eps
        )[0]
        del reference_output

        for source_start in range(0, length, sources_per_batch):
            source_stop = min(length, source_start + sources_per_batch)
            source_positions = torch.arange(
                source_start, source_stop, dtype=torch.long, device=device
            )
            source_reference = reference_base_ids[source_positions]

            rows: List[int] = []
            alternatives: List[int] = []
            for local_index, reference_id in enumerate(source_reference.tolist()):
                for alternative_id in base_token_ids.tolist():
                    if alternative_id != reference_id:
                        rows.append(local_index)
                        alternatives.append(alternative_id)

            local_sources = torch.tensor(rows, dtype=torch.long, device=device)
            alternative_ids = torch.tensor(alternatives, dtype=torch.long, device=device)
            batch_sources = source_positions[local_sources]
            batch_tokens = reference_tokens.repeat(len(batch_sources), 1)
            batch_tokens[
                torch.arange(len(batch_sources), device=device), batch_sources + 1
            ] = alternative_ids

            with autocast_context():
                output = model(batch_tokens)
            mutant_log_odds = log2_odds_from_logits(
                biological_logits(output, length), eps=eps
            )
            # [variant, target, target-base] -> [variant, target]
            variant_effect = (mutant_log_odds - reference_log_odds).abs().amax(dim=-1)
            n_sources = source_stop - source_start
            source_effect = variant_effect.reshape(n_sources, 3, length).amax(dim=1)
            dependency[source_start:source_stop] = source_effect.cpu()

            # Per-substitution influence is the mean effect on all other targets.
            per_variant = variant_effect.reshape(n_sources, 3, length)
            for local_source in range(n_sources):
                absolute_source = source_start + local_source
                if length > 1:
                    summed = per_variant[local_source].sum(dim=-1)
                    summed -= per_variant[local_source, :, absolute_source]
                    variant_influence[absolute_source] = (summed / (length - 1)).cpu()

            mutation_batches += 1
            max_batch_size_used = max(max_batch_size_used, len(batch_sources))
            del output, mutant_log_odds, variant_effect, per_variant, batch_tokens

    dependency.fill_diagonal_(0.0)
    elapsed = time.perf_counter() - started
    dependency_np = dependency.numpy()
    off_diagonal = dependency_np[~np.eye(length, dtype=bool)] if length > 1 else np.zeros(1)
    source_influence = (
        dependency_np.sum(axis=1) / max(length - 1, 1)
    )
    target_sensitivity = (
        dependency_np.sum(axis=0) / max(length - 1, 1)
    )
    gpu_name = (
        torch.cuda.get_device_name(device) if device.type == "cuda" else "CPU"
    )
    forward_sequences = 1 + 3 * length
    stats = ScanStats(
        sequence_id=sequence_id,
        length=length,
        forward_sequences=forward_sequences,
        mutation_batches=mutation_batches,
        max_batch_size_used=max_batch_size_used,
        elapsed_seconds=elapsed,
        sequences_per_second=forward_sequences / max(elapsed, 1e-9),
        nucleotide_predictions_per_second=(forward_sequences * length) / max(elapsed, 1e-9),
        dependency_min=float(np.min(off_diagonal)),
        dependency_median=float(np.median(off_diagonal)),
        dependency_p99=float(np.quantile(off_diagonal, 0.99)),
        dependency_max=float(np.max(off_diagonal)),
        mean_source_influence=float(np.mean(source_influence)),
        mean_target_sensitivity=float(np.mean(target_sensitivity)),
        device=str(device),
        gpu_name=gpu_name,
        torch_version=torch.__version__,
        score_definition=(
            "max over 3 source substitutions and 4 target bases of absolute "
            "change in log2 odds; directed source-to-target; self edges zero"
        ),
    )
    return dependency_np, variant_influence.numpy(), stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--sequence")
    input_group.add_argument("--fasta", type=Path)
    parser.add_argument("--sequence-id", default="sequence")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--max-batch-size", type=int, default=192)
    parser.add_argument("--attention-element-budget", type=int, default=5_000_000)
    parser.add_argument("--eps", type=float, default=1e-6)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.max_batch_size < 3:
        raise ValueError("--max-batch-size must be at least 3")
    if args.attention_element_budget < 3:
        raise ValueError("--attention-element-budget must be positive")
    if not 0.0 < args.eps < 0.1:
        raise ValueError("--eps must be between 0 and 0.1")

    sys.path.insert(0, str(args.source_dir.resolve()))
    from ribogram import Architecture, load_model_and_alphabet

    if args.fasta is not None:
        sequence_id, sequence = read_single_fasta(args.fasta)
    else:
        sequence_id, sequence = args.sequence_id, normalize_rna(args.sequence)

    device = torch.device(args.device)
    if device.type == "cuda":
        torch.cuda.set_device(device)
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA is unavailable")

    alphabet, converter, model = load_model_and_alphabet(
        args.checkpoint,
        architecture=Architecture.RIBOGRAM,
    )
    model = model.to(device).eval()
    dependency, variant_influence, stats = compute_dependency_map(
        model,
        alphabet,
        converter,
        sequence_id,
        sequence,
        device=device,
        max_batch_size=args.max_batch_size,
        attention_element_budget=args.attention_element_budget,
        eps=args.eps,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        sequence_id=np.asarray(sequence_id),
        sequence=np.asarray(sequence),
        dependency=dependency.astype(np.float16),
        variant_influence=variant_influence.astype(np.float16),
    )
    stats_path = args.output.with_suffix(args.output.suffix + ".json")
    stats_path.write_text(
        json.dumps(asdict(stats), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(asdict(stats), indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
