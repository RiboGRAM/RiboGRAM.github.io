#!/usr/bin/env python3
"""Export null-calibrated sparse dependency arcs for the web atlas."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from build_dependency_null_thresholds import (
    DISTANCE_EDGES,
    distance_label,
    gc_label,
    key,
    length_label,
)


def threshold_for(
    thresholds: Dict[str, Dict[str, float]],
    *,
    length_bin: str,
    gc_bin: str,
    distance_bin: str,
) -> Tuple[Dict[str, float], str]:
    candidates = (
        key(length_bin, gc_bin, distance_bin),
        key(length_bin, "*", distance_bin),
        key("*", "*", distance_bin),
    )
    for candidate in candidates:
        if candidate in thresholds:
            return thresholds[candidate], candidate
    raise KeyError(
        f"no null threshold for length={length_bin}, gc={gc_bin}, distance={distance_bin}"
    )


def round_score(value: float) -> float:
    return round(float(value), 5)


def write_gzip_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = (json.dumps(payload, separators=(",", ":"), sort_keys=True) + "\n").encode("utf-8")
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as handle:
            handle.write(encoded)
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--thresholds", type=Path, required=True)
    parser.add_argument("--web-dir", type=Path, required=True)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--max-arcs-per-family", type=int, default=250)
    parser.add_argument("--min-separation", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = [
        json.loads(line)
        for line in args.manifest.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    null_payload = json.loads(args.thresholds.read_text(encoding="utf-8"))
    thresholds: Dict[str, Dict[str, float]] = null_payload["thresholds"]
    rows: List[Dict[str, object]] = []
    missing: List[str] = []
    totals = {
        "directed_edges_tested": 0,
        "undirected_pairs_tested": 0,
        "q99_pairs": 0,
        "q999_pairs": 0,
        "reference_pairs": 0,
        "reference_pairs_q99": 0,
        "reference_pairs_q999": 0,
        "nonreference_pairs_q999": 0,
        "exported_arcs": 0,
    }
    for ordinal, record in enumerate(manifest, start=1):
        family = str(record["family_accession"])
        path = args.results / f"{family}.npz"
        if not path.is_file():
            missing.append(family)
            continue
        archive = np.load(path)
        sequence = str(archive["sequence"].item())
        dependency = np.asarray(archive["dependency"], dtype=np.float32)
        length = len(sequence)
        if dependency.shape != (length, length):
            raise ValueError(f"{family}: dependency shape mismatch")

        length_bin = length_label(length)
        gc_bin = gc_label(sequence)
        known_pairs = {
            tuple(sorted((int(left), int(right))))
            for left, right in record["consensus_pairs"]
            if abs(int(right) - int(left)) >= args.min_separation
        }
        left, right = np.triu_indices(length, k=max(1, args.min_separation))
        forward = dependency[left, right].astype(np.float64)
        reverse = dependency[right, left].astype(np.float64)
        maximum = np.maximum(forward, reverse)
        q99 = np.empty(len(left), dtype=np.float64)
        q999 = np.empty(len(left), dtype=np.float64)
        threshold_cells: Dict[str, Dict[str, object]] = {}
        separation = right - left
        separation_bins = np.searchsorted(
            np.asarray(DISTANCE_EDGES[1:]), separation, side="right"
        )
        for bin_index in np.unique(separation_bins):
            selected = separation_bins == bin_index
            example_distance = int(separation[selected][0])
            distance_bin = distance_label(example_distance)
            null, source_cell = threshold_for(
                thresholds,
                length_bin=length_bin,
                gc_bin=gc_bin,
                distance_bin=distance_bin,
            )
            q99[selected] = float(null["q0.99"])
            q999[selected] = float(null["q0.999"])
            threshold_cells[distance_bin] = {
                "q99": round_score(float(null["q0.99"])),
                "q999": round_score(float(null["q0.999"])),
                "null_edges": int(null["count"]),
                "source_cell": source_cell,
            }

        pass99 = maximum >= q99
        pass999 = maximum >= q999
        is_reference = np.asarray(
            [(int(i), int(j)) in known_pairs for i, j in zip(left, right)], dtype=bool
        )
        candidate_indices = np.flatnonzero(pass99)
        # Null-normalized excess gives comparable ordering across distances.
        normalized = maximum[candidate_indices] / np.maximum(q99[candidate_indices], 1e-12)
        order = candidate_indices[np.argsort(-normalized, kind="mergesort")]
        kept = order[: args.max_arcs_per_family]

        arcs = []
        for index in kept:
            tier = 2 if pass999[index] else 1
            arcs.append(
                [
                    int(left[index]),
                    int(right[index]),
                    round_score(forward[index]),
                    round_score(reverse[index]),
                    tier,
                    int(is_reference[index]),
                ]
            )
        reference_total = int(is_reference.sum())
        reference_q99 = int((pass99 & is_reference).sum())
        reference_q999 = int((pass999 & is_reference).sum())
        family_row: Dict[str, object] = {
            "family_accession": family,
            "family_id": record.get("family_id", ""),
            "rfam_type": record.get("rfam_type", "unknown"),
            "sequence_accession": record["sequence_accession"],
            "length": length,
            "reference_pairs": reference_total,
            "reference_pairs_q99": reference_q99,
            "reference_pairs_q999": reference_q999,
            "reference_recall_q99": reference_q99 / max(reference_total, 1),
            "reference_recall_q999": reference_q999 / max(reference_total, 1),
            "q99_pairs": int(pass99.sum()),
            "q999_pairs": int(pass999.sum()),
            "nonreference_pairs_q999": int((pass999 & ~is_reference).sum()),
            "exported_arcs": len(arcs),
            "truncated_export": bool(len(order) > len(kept)),
        }
        rows.append(family_row)

        web_payload = {
            "schema": "ribogram-substitution-dependency-v1",
            "family_accession": family,
            "family_id": record.get("family_id", ""),
            "rfam_type": record.get("rfam_type", "unknown"),
            "sequence_accession": record["sequence_accession"],
            "sequence": sequence,
            "length": length,
            "score_definition": "max over 3 source substitutions and 4 target bases of absolute change in log2 odds",
            "edge_orientation": "[lower_position, higher_position, lower_to_higher, higher_to_lower, null_tier, is_rfam_pair]; positions are zero-based",
            "null_tiers": {"1": ">= matched dinucleotide-null q99", "2": ">= matched dinucleotide-null q99.9"},
            "thresholds_by_separation": threshold_cells,
            "reference_pairs": sorted([list(pair) for pair in known_pairs]),
            "arcs": arcs,
            "summary": family_row,
            "interpretation_boundary": "model substitution dependency; not a demonstrated physical interaction",
        }
        write_gzip_json(args.web_dir / f"{family}.json.gz", web_payload)

        totals["directed_edges_tested"] += int(2 * len(left))
        totals["undirected_pairs_tested"] += int(len(left))
        totals["q99_pairs"] += int(pass99.sum())
        totals["q999_pairs"] += int(pass999.sum())
        totals["reference_pairs"] += reference_total
        totals["reference_pairs_q99"] += reference_q99
        totals["reference_pairs_q999"] += reference_q999
        totals["nonreference_pairs_q999"] += int((pass999 & ~is_reference).sum())
        totals["exported_arcs"] += len(arcs)
        if ordinal % 250 == 0:
            print(json.dumps({"processed": ordinal, "exported": len(rows)}), flush=True)

    args.rows.parent.mkdir(parents=True, exist_ok=True)
    with args.rows.open("wt", encoding="utf-8") as handle:
        for row in sorted(rows, key=lambda item: str(item["family_accession"])):
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    summary = {
        "families_in_manifest": len(manifest),
        "families_exported": len(rows),
        "missing": missing,
        "max_arcs_per_family": args.max_arcs_per_family,
        "min_separation": args.min_separation,
        **totals,
        "reference_recall_q99": totals["reference_pairs_q99"] / max(totals["reference_pairs"], 1),
        "reference_recall_q999": totals["reference_pairs_q999"] / max(totals["reference_pairs"], 1),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
