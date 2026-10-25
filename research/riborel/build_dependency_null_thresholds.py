#!/usr/bin/env python3
"""Estimate distance-, length- and GC-matched dependency null thresholds."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np


DISTANCE_EDGES = (2, 4, 8, 16, 32, 64, 128, 256, 513, 1025)
LENGTH_EDGES = (20, 81, 161, 321, 513, 1025)
GC_EDGES = (0.0, 0.35, 0.50, 0.65, 1.0000001)


def interval_label(value: float, edges: Sequence[float], *, decimals: int = 0) -> str:
    for lower, upper in zip(edges[:-1], edges[1:]):
        if lower <= value < upper:
            if decimals:
                return f"{lower:.{decimals}f}-{min(upper, 1.0):.{decimals}f}"
            return f"{int(lower)}-{int(upper - 1)}"
    return "out-of-range"


def length_label(length: int) -> str:
    return interval_label(float(length), LENGTH_EDGES)


def distance_label(distance: int) -> str:
    return interval_label(float(distance), DISTANCE_EDGES)


def gc_label(sequence: str) -> str:
    fraction = sum(base in "GC" for base in sequence) / len(sequence)
    return interval_label(fraction, GC_EDGES, decimals=2)


def key(length: str, gc: str, distance: str) -> str:
    return f"length={length}|gc={gc}|distance={distance}"


def summarize(values: np.ndarray) -> Dict[str, float]:
    probabilities = (0.5, 0.9, 0.95, 0.99, 0.995, 0.999, 1.0)
    quantile = np.quantile(values, probabilities)
    output: Dict[str, float] = {
        "count": int(len(values)),
        "mean": float(np.mean(values)),
    }
    output.update(
        {f"q{probability:g}": float(value) for probability, value in zip(probabilities, quantile)}
    )
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = [
        json.loads(line)
        for line in args.manifest.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    pools: Dict[str, List[np.ndarray]] = defaultdict(list)
    scanned: List[str] = []
    missing: List[str] = []
    total_edges = 0
    for record in manifest:
        family = str(record["family_accession"])
        path = args.results / f"{family}.npz"
        if not path.is_file():
            missing.append(family)
            continue
        archive = np.load(path)
        sequence = str(archive["sequence"].item())
        dependency = np.asarray(archive["dependency"], dtype=np.float32)
        if dependency.shape != (len(sequence), len(sequence)):
            raise ValueError(f"{family}: dependency shape mismatch")
        length_bin = length_label(len(sequence))
        gc_bin = gc_label(sequence)
        source, target = np.nonzero(~np.eye(len(sequence), dtype=bool))
        separation = np.abs(source - target)
        valid = separation >= 2
        source, target, separation = source[valid], target[valid], separation[valid]
        scores = dependency[source, target]
        distance_bins = np.searchsorted(
            np.asarray(DISTANCE_EDGES[1:]), separation, side="right"
        )
        for index in np.unique(distance_bins):
            selected = scores[distance_bins == index]
            example_distance = int(separation[distance_bins == index][0])
            distance_bin = distance_label(example_distance)
            pools[key(length_bin, gc_bin, distance_bin)].append(selected)
            pools[key(length_bin, "*", distance_bin)].append(selected)
            pools[key("*", "*", distance_bin)].append(selected)
            total_edges += int(len(selected))
        scanned.append(family)

    thresholds: Dict[str, Dict[str, float]] = {}
    for cell, arrays in sorted(pools.items()):
        values = np.concatenate(arrays).astype(np.float64, copy=False)
        thresholds[cell] = summarize(values)
    payload = {
        "method": {
            "null_sequences": "exact directed-dinucleotide-preserving Eulerian shuffles",
            "score": "directed max-absolute log2-odds substitution dependency",
            "matching": "length bin, GC bin and absolute source-target separation bin; length-only and global fallbacks included",
            "distance_edges": list(DISTANCE_EDGES),
            "length_edges": list(LENGTH_EDGES),
            "gc_edges": list(GC_EDGES),
        },
        "null_records_in_manifest": len(manifest),
        "null_records_scanned": len(scanned),
        "directed_edges_at_exact_cells": total_edges,
        "missing": missing,
        "thresholds": thresholds,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "null_records_scanned": len(scanned),
                "missing": len(missing),
                "threshold_cells": len(thresholds),
                "directed_edges_at_exact_cells": total_edges,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
