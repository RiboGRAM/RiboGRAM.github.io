#!/usr/bin/env python3
"""Aggregate Rfam structure recovery from substitution dependency maps."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import numpy as np

from evaluate_rfam_dependency import evaluate


def quantiles(values: List[float]) -> Dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    probabilities = [0.0, 0.25, 0.5, 0.75, 1.0]
    result = np.quantile(array, probabilities)
    return {str(p): float(v) for p, v in zip(probabilities, result)}


def summarize(rows: List[Dict[str, object]]) -> Dict[str, object]:
    if not rows:
        return {"records": 0}
    metric_names = [
        "roc_auc",
        "average_precision",
        "precision_at_reference_pair_count",
        "precision_at_L",
        "precision_at_top_1_percent",
        "pair_score_enrichment",
        "directed_asymmetry_mean",
        "directed_reciprocity_correlation",
    ]
    output: Dict[str, object] = {
        "records": len(rows),
        "total_nucleotides": sum(int(row["length"]) for row in rows),
        "total_reference_pairs": sum(int(row["reference_pairs"]) for row in rows),
        "length_quantiles": quantiles([float(row["length"]) for row in rows]),
    }
    for metric in metric_names:
        values = [float(row[metric]) for row in rows if np.isfinite(float(row[metric]))]
        if values:
            output[f"{metric}_quantiles"] = quantiles(values)
            output[f"{metric}_mean"] = float(np.mean(values))
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--min-separation", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = [json.loads(line) for line in args.manifest.read_text(encoding="utf-8").splitlines() if line.strip()]
    rows: List[Dict[str, object]] = []
    missing = []
    failed = []
    for record in manifest:
        family = str(record["family_accession"])
        path = args.results / f"{family}.npz"
        if not path.is_file():
            missing.append(family)
            continue
        try:
            archive = np.load(path)
            sequence = str(archive["sequence"].item())
            metrics = evaluate(
                np.asarray(archive["dependency"], dtype=np.float32),
                sequence,
                [tuple(pair) for pair in record["consensus_pairs"]],
                min_separation=args.min_separation,
            )
            row = {
                "family_accession": family,
                "family_id": record["family_id"],
                "rfam_type": record["rfam_type"],
                "sequence_accession": record["sequence_accession"],
                **metrics,
            }
            rows.append(row)
        except Exception as error:
            failed.append({"family_accession": family, "error": str(error)})
    rows.sort(key=lambda row: str(row["family_accession"]))
    args.rows.parent.mkdir(parents=True, exist_ok=True)
    with args.rows.open("wt", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    by_type: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    by_length: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_type[str(row["rfam_type"])].append(row)
        length = int(row["length"])
        label = "20-80" if length <= 80 else "81-160" if length <= 160 else "161-320" if length <= 320 else "321-512" if length <= 512 else "513-1024"
        by_length[label].append(row)
    payload = {
        "overall": summarize(rows),
        "by_rfam_type": {
            key: summarize(value)
            for key, value in sorted(by_type.items())
            if len(value) >= 2
        },
        "by_length": {key: summarize(value) for key, value in sorted(by_length.items())},
        "missing_count": len(missing),
        "missing_families": missing[:100],
        "failed": failed,
        "min_sequence_separation": args.min_separation,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
