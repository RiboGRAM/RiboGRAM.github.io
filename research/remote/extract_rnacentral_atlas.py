#!/usr/bin/env python3
"""Stream one RNAcentral pilot shard through RiboGRAM.

The extractor writes sequence embeddings, released 13/88-class predictions,
and compact per-nucleotide attention tracks.  It never writes dense LxL maps.
The attention implementation follows the RiboGRAM manuscript: head-averaged
APC received attention for layers 7 and 12, plus symmetrized/APC-adjusted
layer-12 head-20 scores for training-free base-pair inference.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import time
from contextlib import nullcontext
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Sequence, Tuple

import numpy as np
import torch
from numpy.lib.format import open_memmap


TAXA_PATTERN = re.compile(r"\bfound in (\d+) tax(?:on|a)\b", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--head13", type=Path, required=True)
    parser.add_argument("--head88", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--records", type=int, required=True)
    parser.add_argument("--nt", type=int, required=True)
    parser.add_argument("--token-budget", type=int, default=16384)
    parser.add_argument("--max-batch-size", type=int, default=128)
    parser.add_argument("--buffer-records", type=int, default=4096)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--progress-seconds", type=float, default=30.0)
    args = parser.parse_args()
    if args.records < 1 or args.nt < 1:
        parser.error("--records and --nt must be positive")
    if args.token_budget < 1 or args.max_batch_size < 1:
        parser.error("batch limits must be positive")
    return args


def atomic_json(path: Path, payload: Mapping[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def read_fasta(path: Path) -> Iterator[Tuple[str, str, str, int, int]]:
    header = ""
    chunks: List[str] = []

    def build() -> Tuple[str, str, str, int, int]:
        sequence = "".join(chunks).upper().replace("T", "U")
        fields = header[1:].split()
        identifier = fields[0]
        rna_type = fields[1] if len(fields) > 1 else "unknown"
        match = TAXA_PATTERN.search(header)
        taxa_count = int(match.group(1)) if match else 0
        ambiguous = sum(base not in "ACGU" for base in sequence)
        return identifier, sequence, rna_type, taxa_count, ambiguous

    with path.open("rt", encoding="ascii", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield build()
                header = line
                chunks = []
            else:
                chunks.append(line)
        if header:
            yield build()


def buffered_sorted_records(
    path: Path, buffer_records: int
) -> Iterator[Tuple[str, str, str, int, int]]:
    buffer: List[Tuple[str, str, str, int, int]] = []
    for record in read_fasta(path):
        buffer.append(record)
        if len(buffer) >= buffer_records:
            buffer.sort(key=lambda item: (len(item[1]), item[0]))
            yield from buffer
            buffer.clear()
    if buffer:
        buffer.sort(key=lambda item: (len(item[1]), item[0]))
        yield from buffer


def dynamic_batches(
    records: Iterable[Tuple[str, str, str, int, int]],
    token_budget: int,
    max_batch_size: int,
) -> Iterator[List[Tuple[str, str, str, int, int]]]:
    batch: List[Tuple[str, str, str, int, int]] = []
    maximum = 0
    for record in records:
        length_with_specials = len(record[1]) + 2
        candidate_maximum = max(maximum, length_with_specials)
        if batch and (
            candidate_maximum * (len(batch) + 1) > token_budget
            or len(batch) >= max_batch_size
        ):
            yield batch
            batch = []
            maximum = 0
        batch.append(record)
        maximum = max(maximum, length_with_specials)
    if batch:
        yield batch


def biological_mask(tokens: torch.Tensor, alphabet) -> torch.Tensor:
    return (
        tokens.ne(alphabet.padding_idx)
        & tokens.ne(alphabet.cls_idx)
        & tokens.ne(alphabet.eos_idx)
    )


class ClassificationEnsemble:
    def __init__(self, checkpoint: Path, head_class, device: torch.device):
        try:
            bundle = torch.load(str(checkpoint), map_location="cpu", weights_only=False)
        except TypeError:
            bundle = torch.load(str(checkpoint), map_location="cpu")
        self.class_names = list(bundle["class_names"])
        self.members = []
        for artifact in bundle["members"]:
            config = artifact["head"]
            head = head_class(
                input_dim=int(artifact["input_dim"]),
                hidden_dim=int(config["hidden"]),
                depth=int(config["depth"]),
                dropout=float(config["dropout"]),
                num_classes=int(artifact["num_classes"]),
            )
            head.load_state_dict(artifact["model_state_dict"], strict=True)
            self.members.append(
                (
                    head.to(device).eval(),
                    tuple(artifact["feature_names"]),
                    artifact["center"].to(device),
                    artifact["scale"].to(device),
                )
            )

    def predict(self, feature_map: Mapping[str, torch.Tensor]) -> torch.Tensor:
        probabilities = []
        for head, names, center, scale in self.members:
            features = torch.cat([feature_map[name] for name in names], dim=1)
            standardized = (features - center) / scale.clamp_min(1e-5)
            probabilities.append(head(standardized).softmax(dim=-1))
        return torch.stack(probabilities).mean(dim=0)


def apc_received(attention: torch.Tensor, valid: torch.Tensor) -> torch.Tensor:
    """Return manuscript-normalized received attention for [B,H,L,L]."""
    pair_mask = valid[:, :, None] & valid[:, None, :]
    lengths = valid.sum(dim=1).clamp_min(1).float()
    matrix = attention.float().mean(dim=1) * pair_mask
    row_mean = matrix.sum(dim=-1, keepdim=True) / lengths[:, None, None]
    column_mean = matrix.sum(dim=-2, keepdim=True) / lengths[:, None, None]
    global_mean = matrix.sum(dim=(-2, -1), keepdim=True) / lengths.square()[:, None, None]
    corrected = (
        matrix - row_mean * column_mean / global_mean.clamp_min(1e-12)
    ).clamp_min(0.0) * pair_mask
    row_total = corrected.sum(dim=-1, keepdim=True)
    normalized = torch.where(
        row_total > 0,
        corrected / row_total.clamp_min(1e-12),
        torch.zeros_like(corrected),
    )
    received = normalized.sum(dim=-2) * valid
    fraction = received / received.sum(dim=-1, keepdim=True).clamp_min(1e-12)
    return fraction * lengths[:, None]


def pairing_top_partner(
    attention: torch.Tensor, valid: torch.Tensor
) -> Tuple[torch.Tensor, torch.Tensor]:
    """Return top symmetrized/APC L12H20 partner and score per nucleotide."""
    pair_mask = valid[:, :, None] & valid[:, None, :]
    lengths = valid.sum(dim=1).clamp_min(1).float()
    matrix = attention.float()
    matrix = 0.5 * (matrix + matrix.transpose(-1, -2))
    matrix = matrix * pair_mask
    row_mean = matrix.sum(dim=-1, keepdim=True) / lengths[:, None, None]
    column_mean = matrix.sum(dim=-2, keepdim=True) / lengths[:, None, None]
    global_mean = matrix.sum(dim=(-2, -1), keepdim=True) / lengths.square()[:, None, None]
    corrected = (
        matrix - row_mean * column_mean / global_mean.clamp_min(1e-12)
    ).clamp_min(0.0) * pair_mask
    diagonal = torch.eye(
        corrected.shape[-1], dtype=torch.bool, device=corrected.device
    )[None, :, :]
    candidates = corrected.masked_fill(~pair_mask | diagonal, float("-inf"))
    score, partner = candidates.max(dim=-1)
    score = torch.where(valid & torch.isfinite(score), score, torch.zeros_like(score))
    partner = torch.where(valid, partner, torch.full_like(partner, 65535))
    return score, partner


def create_or_open_arrays(
    output_dir: Path, records: int, nt: int, resume: bool
) -> Dict[str, np.memmap]:
    specifications = {
        "embeddings": ((records, 1280), np.float16),
        "ids": ((records,), "S16"),
        "rna_types": ((records,), "S64"),
        "lengths": ((records,), np.uint16),
        "taxa_counts": ((records,), np.uint32),
        "ambiguous_counts": ((records,), np.uint16),
        "ncrna13_id": ((records,), np.uint8),
        "ncrna13_confidence": ((records,), np.float16),
        "ncrna88_id": ((records,), np.uint8),
        "ncrna88_confidence": ((records,), np.float16),
        "offsets": ((records + 1,), np.uint64),
        "attention_tracks": ((nt, 3), np.float16),
        "pair_partner": ((nt,), np.uint16),
    }
    arrays: Dict[str, np.memmap] = {}
    for name, (shape, dtype) in specifications.items():
        path = output_dir / f"{name}.npy"
        mode = "r+" if resume else "w+"
        arrays[name] = open_memmap(path, mode=mode, dtype=dtype, shape=shape)
    return arrays


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    success_path = args.output_dir / "SUCCESS"
    if success_path.exists():
        print(f"Already complete: {success_path}", flush=True)
        return 0

    sys.path.insert(0, str(args.source_dir))
    from infer import NcRNAClassificationHead, pooled_ncrna_features
    from ribogram.pretrained import load_model_and_alphabet

    device = torch.device(args.device)
    if device.type == "cuda":
        torch.cuda.set_device(device)
        torch.backends.cuda.matmul.allow_tf32 = True
    print(f"[{time.strftime('%F %T')}] loading model on {device}", flush=True)
    alphabet, converter, model = load_model_and_alphabet(args.checkpoint)
    model = model.to(device).eval()
    classifier13 = ClassificationEnsemble(
        args.head13, NcRNAClassificationHead, device
    )
    classifier88 = ClassificationEnsemble(
        args.head88, NcRNAClassificationHead, device
    )

    class_metadata = {
        "ncrna13": classifier13.class_names,
        "ncrna88": classifier88.class_names,
        "attention_tracks": [
            "layer7_head_mean_apc_received",
            "layer12_head_mean_apc_received",
            "layer12_head20_symmetrized_apc_top_score",
        ],
        "pair_partner_indexing": "zero_based; 65535 means missing",
        "pairing_threshold_from_manuscript": 0.4,
        "embedding": "RiboGRAM final layer mean over biological nucleotides",
    }
    atomic_json(args.output_dir / "schema.json", class_metadata)

    progress_path = args.output_dir / "progress.json"
    record_cursor = 0
    nt_cursor = 0
    if args.resume and progress_path.exists():
        progress = json.loads(progress_path.read_text(encoding="utf-8"))
        record_cursor = int(progress["records_complete"])
        nt_cursor = int(progress["nt_complete"])
    elif args.resume:
        raise RuntimeError("--resume requested but progress.json is absent")

    arrays = create_or_open_arrays(
        args.output_dir, args.records, args.nt, args.resume
    )
    if not args.resume:
        arrays["offsets"][0] = 0

    ordered_records = dynamic_batches(
        buffered_sorted_records(args.fasta, args.buffer_records),
        args.token_budget,
        args.max_batch_size,
    )
    seen_records = 0
    started = time.time()
    last_progress = started

    inference_context = torch.inference_mode()
    with inference_context:
        for batch in ordered_records:
            batch_size = len(batch)
            if seen_records + batch_size <= record_cursor:
                seen_records += batch_size
                continue
            if seen_records < record_cursor:
                batch = batch[record_cursor - seen_records :]
                batch_size = len(batch)
                seen_records = record_cursor

            labels = [item[0] for item in batch]
            sequences = [item[1] for item in batch]
            records_for_converter = list(zip(labels, sequences))
            _, _, tokens = converter(records_for_converter)
            tokens = tokens.to(device, non_blocking=True)
            mask = biological_mask(tokens, alphabet)
            lengths = mask.sum(dim=1)
            maximum_length = int(lengths.max().item())

            autocast = (
                torch.cuda.amp.autocast(dtype=torch.float16)
                if device.type == "cuda"
                else nullcontext()
            )
            with autocast:
                output = model(
                    tokens,
                    need_head_weights=True,
                    attention_layers=(7, 12),
                )
            if tuple(output["attention_layer_indices"]) != (7, 12):
                raise RuntimeError("Unexpected attention layer order")

            representation = output["representation"]
            feature_map = pooled_ncrna_features(representation, mask)
            embedding = feature_map["mean"]
            probability13 = classifier13.predict(feature_map)
            probability88 = classifier88.predict(feature_map)

            valid = (
                torch.arange(maximum_length, device=device)[None, :]
                < lengths[:, None]
            )
            attention = output["attentions"][:, :, :, 1 : maximum_length + 1, 1 : maximum_length + 1]
            received7 = apc_received(attention[:, 0], valid)
            received12 = apc_received(attention[:, 1], valid)
            pair_score, pair_partner = pairing_top_partner(
                attention[:, 1, 19], valid
            )
            track_tensor = torch.stack(
                (received7, received12, pair_score), dim=-1
            )
            flat_tracks = track_tensor[valid].half().cpu().numpy()
            flat_partners = (
                pair_partner[valid].to(torch.int32).cpu().numpy().astype(np.uint16)
            )
            embedding_np = embedding.half().cpu().numpy()
            lengths_np = lengths.to(torch.int64).cpu().numpy()
            ids13 = probability13.argmax(dim=1)
            conf13 = probability13.gather(1, ids13[:, None]).squeeze(1)
            ids88 = probability88.argmax(dim=1)
            conf88 = probability88.gather(1, ids88[:, None]).squeeze(1)

            next_record = record_cursor + batch_size
            next_nt = nt_cursor + int(lengths_np.sum())
            if next_record > args.records or next_nt > args.nt:
                raise RuntimeError("Observed records/nt exceed manifest allocation")

            record_slice = slice(record_cursor, next_record)
            nt_slice = slice(nt_cursor, next_nt)
            arrays["embeddings"][record_slice] = embedding_np
            arrays["ids"][record_slice] = np.asarray(labels, dtype="S16")
            arrays["rna_types"][record_slice] = np.asarray(
                [item[2] for item in batch], dtype="S64"
            )
            arrays["lengths"][record_slice] = lengths_np.astype(np.uint16)
            arrays["taxa_counts"][record_slice] = np.asarray(
                [item[3] for item in batch], dtype=np.uint32
            )
            arrays["ambiguous_counts"][record_slice] = np.asarray(
                [item[4] for item in batch], dtype=np.uint16
            )
            arrays["ncrna13_id"][record_slice] = ids13.byte().cpu().numpy()
            arrays["ncrna13_confidence"][record_slice] = conf13.half().cpu().numpy()
            arrays["ncrna88_id"][record_slice] = ids88.byte().cpu().numpy()
            arrays["ncrna88_confidence"][record_slice] = conf88.half().cpu().numpy()
            cumulative = np.cumsum(lengths_np, dtype=np.uint64) + np.uint64(nt_cursor)
            arrays["offsets"][record_cursor + 1 : next_record + 1] = cumulative
            arrays["attention_tracks"][nt_slice] = flat_tracks
            arrays["pair_partner"][nt_slice] = flat_partners

            record_cursor = next_record
            nt_cursor = next_nt
            seen_records += batch_size
            now = time.time()
            if now - last_progress >= args.progress_seconds:
                elapsed = max(now - started, 1e-9)
                payload = {
                    "records_complete": record_cursor,
                    "records_expected": args.records,
                    "nt_complete": nt_cursor,
                    "nt_expected": args.nt,
                    "records_per_second": record_cursor / elapsed,
                    "nt_per_second": nt_cursor / elapsed,
                    "updated_unix": time.time(),
                }
                atomic_json(progress_path, payload)
                print(json.dumps(payload, sort_keys=True), flush=True)
                last_progress = now

    if record_cursor != args.records or nt_cursor != args.nt:
        raise RuntimeError(
            f"Manifest mismatch records={record_cursor}/{args.records} "
            f"nt={nt_cursor}/{args.nt}"
        )
    for array in arrays.values():
        array.flush()
    final_progress = {
        "records_complete": record_cursor,
        "records_expected": args.records,
        "nt_complete": nt_cursor,
        "nt_expected": args.nt,
        "elapsed_seconds": time.time() - started,
        "complete": True,
    }
    atomic_json(progress_path, final_progress)
    success_path.write_text(json.dumps(final_progress) + "\n", encoding="utf-8")
    print(json.dumps(final_progress, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
