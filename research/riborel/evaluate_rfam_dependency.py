#!/usr/bin/env python3
"""Evaluate a RiboGRAM dependency map against an Rfam seed structure."""

from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("rt", encoding="utf-8", errors="replace")


def find_stockholm_record(path: Path, accession: str) -> Dict[str, str]:
    family_accession = ""
    family_id = ""
    sequence = None
    consensus = ""
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "//":
                if sequence is not None:
                    return {
                        "family_accession": family_accession,
                        "family_id": family_id,
                        "accession": accession,
                        "aligned_sequence": sequence,
                        "ss_cons": consensus,
                    }
                family_accession = ""
                family_id = ""
                sequence = None
                consensus = ""
                continue
            if line.startswith("#=GF AC"):
                family_accession = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF ID"):
                family_id = line.split(None, 2)[2].strip()
            elif line.startswith("#=GC SS_cons"):
                consensus += line.split(None, 2)[2].strip()
            elif line and not line.startswith("#"):
                fields = line.split()
                if len(fields) >= 2 and fields[0] == accession:
                    sequence = (sequence or "") + fields[1].strip()
    raise KeyError(f"accession not found in Rfam seed: {accession}")


def wuss_pairs(structure: str) -> List[Tuple[int, int]]:
    opening_to_closing = {"(": ")", "[": "]", "{": "}", "<": ">"}
    closing_to_opening = {value: key for key, value in opening_to_closing.items()}
    stacks: Dict[str, List[int]] = {key: [] for key in opening_to_closing}
    letter_stacks: Dict[str, List[int]] = {}
    pairs: List[Tuple[int, int]] = []
    for index, symbol in enumerate(structure):
        if symbol in opening_to_closing:
            stacks[symbol].append(index)
        elif symbol in closing_to_opening:
            opener = closing_to_opening[symbol]
            if stacks[opener]:
                pairs.append((stacks[opener].pop(), index))
        elif "A" <= symbol <= "Z":
            letter_stacks.setdefault(symbol, []).append(index)
        elif "a" <= symbol <= "z":
            opener = symbol.upper()
            if letter_stacks.get(opener):
                pairs.append((letter_stacks[opener].pop(), index))
    return sorted(pairs)


def map_consensus_pairs(aligned_sequence: str, ss_cons: str) -> Tuple[str, List[Tuple[int, int]]]:
    if len(aligned_sequence) != len(ss_cons):
        raise ValueError(
            f"alignment/SS_cons lengths differ: {len(aligned_sequence)} != {len(ss_cons)}"
        )
    column_to_position: Dict[int, int] = {}
    sequence: List[str] = []
    for column, symbol in enumerate(aligned_sequence):
        if symbol not in "-.~":
            column_to_position[column] = len(sequence)
            sequence.append(symbol.upper().replace("T", "U"))
    mapped = []
    for left, right in wuss_pairs(ss_cons):
        if left in column_to_position and right in column_to_position:
            mapped.append((column_to_position[left], column_to_position[right]))
    return "".join(sequence), mapped


def rankdata_average(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=np.float64)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and sorted_values[stop] == sorted_values[start]:
            stop += 1
        average_rank = (start + 1 + stop) / 2.0
        ranks[order[start:stop]] = average_rank
        start = stop
    return ranks


def roc_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    positives = labels.astype(bool)
    n_pos = int(positives.sum())
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = rankdata_average(scores)
    rank_sum = float(ranks[positives].sum())
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    order = np.argsort(-scores, kind="mergesort")
    ranked = labels[order].astype(np.float64)
    n_pos = int(ranked.sum())
    if n_pos == 0:
        return float("nan")
    precision = np.cumsum(ranked) / np.arange(1, len(ranked) + 1)
    return float((precision * ranked).sum() / n_pos)


def precision_at(labels: np.ndarray, scores: np.ndarray, k: int) -> float:
    k = max(1, min(int(k), len(labels)))
    selected = np.argpartition(scores, len(scores) - k)[-k:]
    return float(labels[selected].mean())


def evaluate(
    dependency: np.ndarray,
    sequence: str,
    pairs: Sequence[Tuple[int, int]],
    min_separation: int,
) -> Dict[str, object]:
    length = len(sequence)
    if dependency.shape != (length, length):
        raise ValueError(f"dependency shape {dependency.shape} != {(length, length)}")
    symmetric = np.maximum(dependency, dependency.T)
    rows, columns = np.triu_indices(length, k=max(1, min_separation))
    scores = symmetric[rows, columns].astype(np.float64)
    pair_set = {tuple(sorted(pair)) for pair in pairs if abs(pair[1] - pair[0]) >= min_separation}
    labels = np.asarray(
        [(int(left), int(right)) in pair_set for left, right in zip(rows, columns)],
        dtype=np.int8,
    )
    n_pairs = int(labels.sum())
    base_pair_background = n_pairs / max(len(labels), 1)
    top_pair_count = max(n_pairs, 1)
    top_fraction_k = max(1, int(math.ceil(0.01 * len(labels))))
    canonical = {"AU", "UA", "GC", "CG", "GU", "UG"}
    order = np.argsort(-scores, kind="mergesort")
    top_l = order[: min(length, len(order))]
    canonical_top_l = np.mean(
        [sequence[rows[index]] + sequence[columns[index]] in canonical for index in top_l]
    )
    direction_delta = np.abs(dependency - dependency.T)
    return {
        "length": length,
        "reference_pairs": n_pairs,
        "candidate_pairs": int(len(labels)),
        "pair_background_fraction": float(base_pair_background),
        "roc_auc": roc_auc(labels, scores),
        "average_precision": average_precision(labels, scores),
        "precision_at_reference_pair_count": precision_at(labels, scores, top_pair_count),
        "precision_at_L": precision_at(labels, scores, length),
        "precision_at_top_1_percent": precision_at(labels, scores, top_fraction_k),
        "canonical_fraction_top_L": float(canonical_top_l),
        "median_reference_pair_score": float(np.median(scores[labels.astype(bool)])),
        "median_nonpair_score": float(np.median(scores[~labels.astype(bool)])),
        "pair_score_enrichment": float(
            np.median(scores[labels.astype(bool)]) / max(np.median(scores[~labels.astype(bool)]), 1e-12)
        ),
        "directed_asymmetry_mean": float(
            direction_delta[~np.eye(length, dtype=bool)].mean()
        ),
        "directed_reciprocity_correlation": float(
            np.corrcoef(
                dependency[rows, columns].astype(np.float64),
                dependency[columns, rows].astype(np.float64),
            )[0, 1]
        ),
        "min_sequence_separation": int(min_separation),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rfam-seed", type=Path, required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--dependency", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-separation", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    record = find_stockholm_record(args.rfam_seed, args.accession)
    sequence, pairs = map_consensus_pairs(record["aligned_sequence"], record["ss_cons"])
    archive = np.load(args.dependency)
    scanned_sequence = str(archive["sequence"].item())
    if sequence != scanned_sequence:
        raise ValueError(
            "Rfam ungapped sequence does not match scanned sequence: "
            f"{sequence!r} != {scanned_sequence!r}"
        )
    metrics = evaluate(
        np.asarray(archive["dependency"], dtype=np.float32),
        sequence,
        pairs,
        min_separation=args.min_separation,
    )
    payload = {
        "family_accession": record["family_accession"],
        "family_id": record["family_id"],
        "sequence_accession": args.accession,
        "metrics": metrics,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
