#!/usr/bin/env python3
"""Summarize real and null-calibrated strong RiboGRAM dependency edges."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from build_dependency_null_thresholds import DISTANCE_EDGES, distance_label, gc_label, length_label
from export_dependency_atlas import threshold_for


CANONICAL = {"AU", "UA", "GC", "CG", "GU", "UG"}


def broad_type(value: str) -> str:
    return value.split(";", 1)[0].strip() or "unknown"


def subtype(value: str) -> str:
    parts = [part.strip() for part in value.split(";") if part.strip()]
    return parts[-1] if len(parts) > 1 else (parts[0] if parts else "unknown")


def blank_group() -> Dict[str, object]:
    return {
        "families": 0,
        "candidate_pairs": 0,
        "q99_pairs": 0,
        "q999_pairs": 0,
        "q999_reciprocal_pairs": 0,
        "q999_canonical_pairs": 0,
        "q999_reference_pairs": 0,
        "q999_nonreference_pairs": 0,
        "reference_pairs": 0,
        "reference_pairs_q99": 0,
        "reference_pairs_q999": 0,
        "q999_separation": Counter(),
        "q999_base_combinations": Counter(),
        "q999_direction_lower_to_higher": 0,
    }


def finalize(group: Dict[str, object]) -> Dict[str, object]:
    output = dict(group)
    output["q999_separation"] = dict(sorted(group["q999_separation"].items()))
    output["q999_base_combinations"] = dict(
        group["q999_base_combinations"].most_common()
    )
    candidate = max(int(group["candidate_pairs"]), 1)
    q999 = max(int(group["q999_pairs"]), 1)
    reference = max(int(group["reference_pairs"]), 1)
    output.update(
        q99_pair_fraction=float(group["q99_pairs"] / candidate),
        q999_pair_fraction=float(group["q999_pairs"] / candidate),
        q999_reciprocal_fraction=float(group["q999_reciprocal_pairs"] / q999),
        q999_canonical_fraction=float(group["q999_canonical_pairs"] / q999),
        q999_lower_to_higher_fraction=float(group["q999_direction_lower_to_higher"] / q999),
        reference_recall_q99=float(group["reference_pairs_q99"] / reference),
        reference_recall_q999=float(group["reference_pairs_q999"] / reference),
    )
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--thresholds", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-separation", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = [
        json.loads(line)
        for line in args.manifest.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    thresholds = json.loads(args.thresholds.read_text(encoding="utf-8"))["thresholds"]
    groups: Dict[Tuple[str, str], Dict[str, object]] = defaultdict(blank_group)
    missing: List[str] = []
    for ordinal, record in enumerate(manifest, start=1):
        family = str(record["family_accession"])
        path = args.results / f"{family}.npz"
        if not path.is_file():
            missing.append(family)
            continue
        archive = np.load(path)
        sequence = str(archive["sequence"].item())
        dependency = np.asarray(archive["dependency"], dtype=np.float32)
        length_bin = length_label(len(sequence))
        gc_bin = gc_label(sequence)
        rfam_type = str(record.get("rfam_type", "unknown"))
        group_keys = (
            ("overall", "all"),
            ("broad_type", broad_type(rfam_type)),
            ("subtype", subtype(rfam_type)),
            ("length", length_bin),
        )

        left, right = np.triu_indices(len(sequence), k=max(1, args.min_separation))
        forward = dependency[left, right].astype(np.float64)
        reverse = dependency[right, left].astype(np.float64)
        maximum = np.maximum(forward, reverse)
        separation = right - left
        separation_bins = np.searchsorted(
            np.asarray(DISTANCE_EDGES[1:]), separation, side="right"
        )
        q99 = np.empty(len(left), dtype=np.float64)
        q999 = np.empty(len(left), dtype=np.float64)
        labels = np.empty(len(left), dtype=object)
        for bin_index in np.unique(separation_bins):
            selected = separation_bins == bin_index
            label = distance_label(int(separation[selected][0]))
            null, _ = threshold_for(
                thresholds,
                length_bin=length_bin,
                gc_bin=gc_bin,
                distance_bin=label,
            )
            q99[selected] = float(null["q0.99"])
            q999[selected] = float(null["q0.999"])
            labels[selected] = label
        pass99 = maximum >= q99
        pass999 = maximum >= q999
        both999 = (forward >= q999) & (reverse >= q999)
        lower_dominant = forward >= reverse
        pair_strings = np.asarray(
            [sequence[int(i)] + sequence[int(j)] for i, j in zip(left, right)],
            dtype=object,
        )
        canonical = np.asarray([value in CANONICAL for value in pair_strings], dtype=bool)
        reference_set = {
            tuple(sorted((int(i), int(j))))
            for i, j in record.get("consensus_pairs", [])
            if abs(int(j) - int(i)) >= args.min_separation
        }
        is_reference = np.asarray(
            [(int(i), int(j)) in reference_set for i, j in zip(left, right)], dtype=bool
        )

        for group_key in group_keys:
            group = groups[group_key]
            group["families"] += 1
            group["candidate_pairs"] += int(len(left))
            group["q99_pairs"] += int(pass99.sum())
            group["q999_pairs"] += int(pass999.sum())
            group["q999_reciprocal_pairs"] += int((pass999 & both999).sum())
            group["q999_canonical_pairs"] += int((pass999 & canonical).sum())
            group["q999_reference_pairs"] += int((pass999 & is_reference).sum())
            group["q999_nonreference_pairs"] += int((pass999 & ~is_reference).sum())
            group["reference_pairs"] += int(is_reference.sum())
            group["reference_pairs_q99"] += int((pass99 & is_reference).sum())
            group["reference_pairs_q999"] += int((pass999 & is_reference).sum())
            group["q999_direction_lower_to_higher"] += int((pass999 & lower_dominant).sum())
            for label, count in Counter(labels[pass999].tolist()).items():
                group["q999_separation"][str(label)] += int(count)
            for bases, count in Counter(pair_strings[pass999].tolist()).items():
                group["q999_base_combinations"][str(bases)] += int(count)
        if ordinal % 250 == 0:
            print(json.dumps({"processed": ordinal}), flush=True)

    payload: Dict[str, object] = {
        "method": {
            "strong_edge": "maximum of the two directions exceeds the matched dinucleotide-null threshold",
            "q99": "99th percentile matched by RNA length, GC fraction and sequence separation",
            "q999": "99.9th percentile matched by RNA length, GC fraction and sequence separation",
        },
        "families_in_manifest": len(manifest),
        "missing": missing,
        "overall": finalize(groups[("overall", "all")]),
    }
    for group_type in ("broad_type", "subtype", "length"):
        payload[f"by_{group_type}"] = {
            name: finalize(group)
            for (kind, name), group in sorted(groups.items())
            if kind == group_type
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload["overall"], indent=2, sort_keys=True), flush=True)
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
