#!/usr/bin/env python3
"""Summarize RiboGRAM saturation-mutagenesis dependency architecture.

The primary benchmark asks whether a substitution-derived dependency map
recovers Rfam consensus pairs.  In addition to ordinary ROC/AP, this script
computes a distance-matched pair AUC so that recovery cannot be explained
solely by a preference for short-range edges.  It also quantifies direction,
reciprocity, pseudoknot-associated pairs, composition and paired-position
hubness.  Results are descriptive model diagnostics, not evidence of physical
RNA interactions.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from evaluate_rfam_dependency import evaluate


DISTANCE_EDGES = (2, 4, 8, 16, 32, 64, 128, 256, 513, 1025)
LENGTH_EDGES = (20, 81, 161, 321, 513, 1025)
CANONICAL = {"AU", "UA", "GC", "CG", "GU", "UG"}


def finite(value: object) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def quantiles(values: Iterable[float]) -> Dict[str, float]:
    array = np.asarray([float(value) for value in values if finite(value)], dtype=np.float64)
    if not len(array):
        return {}
    probabilities = (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
    result = np.quantile(array, probabilities)
    return {f"{probability:g}": float(value) for probability, value in zip(probabilities, result)}


def bootstrap_mean_ci(values: Iterable[float], *, seed: int, replicates: int = 2000) -> Dict[str, float]:
    array = np.asarray([float(value) for value in values if finite(value)], dtype=np.float64)
    if not len(array):
        return {}
    rng = np.random.default_rng(seed)
    means = np.empty(replicates, dtype=np.float64)
    # Chunking avoids allocating replicates x families all at once.
    for start in range(0, replicates, 128):
        stop = min(replicates, start + 128)
        indices = rng.integers(0, len(array), size=(stop - start, len(array)))
        means[start:stop] = array[indices].mean(axis=1)
    low, high = np.quantile(means, (0.025, 0.975))
    return {
        "mean": float(array.mean()),
        "ci95_low": float(low),
        "ci95_high": float(high),
        "n": int(len(array)),
    }


def rankdata_average(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=np.float64)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and sorted_values[stop] == sorted_values[start]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * (start + stop - 1) + 1.0
        start = stop
    return ranks


def spearman(left: Sequence[float], right: Sequence[float]) -> float:
    x = np.asarray(left, dtype=np.float64)
    y = np.asarray(right, dtype=np.float64)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return float("nan")
    x_rank = rankdata_average(x[mask])
    y_rank = rankdata_average(y[mask])
    if np.std(x_rank) == 0 or np.std(y_rank) == 0:
        return float("nan")
    return float(np.corrcoef(x_rank, y_rank)[0, 1])


def bin_label(value: int, edges: Sequence[int]) -> str:
    for lower, upper in zip(edges[:-1], edges[1:]):
        if lower <= value < upper:
            return f"{lower}-{upper - 1}"
    return f">={edges[-1]}"


def type_parts(rfam_type: str) -> Tuple[str, str, str]:
    parts = [part.strip() for part in rfam_type.split(";") if part.strip()]
    broad = parts[0] if parts else "unknown"
    subtype = parts[-1] if len(parts) > 1 else broad
    hierarchy = " / ".join(parts) if parts else "unknown"
    return broad, subtype, hierarchy


def gc_ruggedness(sequence: str, window: int = 11) -> float:
    gc = np.asarray([base in "GC" for base in sequence], dtype=np.float64)
    if len(gc) <= 1:
        return 0.0
    kernel = np.ones(min(window, len(gc)), dtype=np.float64)
    numerator = np.convolve(gc, kernel, mode="same")
    denominator = np.convolve(np.ones_like(gc), kernel, mode="same")
    smooth = numerator / denominator
    return float(np.mean(np.abs(np.diff(smooth))))


def crossing_pair_mask(pairs: Sequence[Tuple[int, int]]) -> np.ndarray:
    normalized = [tuple(sorted((int(left), int(right)))) for left, right in pairs]
    crossing = np.zeros(len(normalized), dtype=bool)
    for first, (left, right) in enumerate(normalized):
        for second in range(first + 1, len(normalized)):
            other_left, other_right = normalized[second]
            if (left < other_left < right < other_right) or (
                other_left < left < other_right < right
            ):
                crossing[first] = True
                crossing[second] = True
    return crossing


def distance_matched_pair_percentiles(
    symmetric: np.ndarray,
    pairs: Sequence[Tuple[int, int]],
    *,
    min_separation: int,
) -> Tuple[np.ndarray, List[str]]:
    length = symmetric.shape[0]
    normalized_pairs = {
        tuple(sorted((int(left), int(right))))
        for left, right in pairs
        if abs(int(right) - int(left)) >= min_separation
    }
    rows, columns = np.triu_indices(length, k=max(1, min_separation))
    distance = columns - rows
    scores = symmetric[rows, columns].astype(np.float64)
    labels = np.asarray(
        [(int(left), int(right)) in normalized_pairs for left, right in zip(rows, columns)],
        dtype=bool,
    )
    bins = np.searchsorted(np.asarray(DISTANCE_EDGES[1:]), distance, side="right")
    negative_by_bin: Dict[int, np.ndarray] = {}
    for index in np.unique(bins):
        negative_by_bin[int(index)] = np.sort(scores[(bins == index) & ~labels])

    percentiles: List[float] = []
    labels_out: List[str] = []
    for left, right in sorted(normalized_pairs):
        separation = right - left
        index = int(np.searchsorted(np.asarray(DISTANCE_EDGES[1:]), separation, side="right"))
        negatives = negative_by_bin.get(index, np.empty(0, dtype=np.float64))
        if not len(negatives):
            continue
        score = float(symmetric[left, right])
        below = int(np.searchsorted(negatives, score, side="left"))
        at_or_below = int(np.searchsorted(negatives, score, side="right"))
        percentiles.append((below + at_or_below) / (2.0 * len(negatives)))
        labels_out.append(bin_label(separation, DISTANCE_EDGES))
    return np.asarray(percentiles, dtype=np.float64), labels_out


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator / max(denominator, 1e-12))


def json_safe(value: object) -> object:
    """Recursively replace NaN/inf with null so outputs are strict JSON."""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(float(value)) else None
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def analyze_record(record: Dict[str, object], dependency: np.ndarray, sequence: str) -> Dict[str, object]:
    pairs = [tuple(map(int, pair)) for pair in record["consensus_pairs"]]
    length = len(sequence)
    symmetric = np.maximum(dependency, dependency.T)
    standard = evaluate(dependency, sequence, pairs, min_separation=2)

    normalized_pairs = [
        tuple(sorted(pair)) for pair in pairs if abs(pair[1] - pair[0]) >= 2
    ]
    pair_scores = np.asarray([symmetric[left, right] for left, right in normalized_pairs], dtype=np.float64)
    pair_forward = np.asarray([dependency[left, right] for left, right in normalized_pairs], dtype=np.float64)
    pair_reverse = np.asarray([dependency[right, left] for left, right in normalized_pairs], dtype=np.float64)
    direction_signed = (pair_forward - pair_reverse) / (pair_forward + pair_reverse + 1e-12)
    reciprocity = np.minimum(pair_forward, pair_reverse) / (
        np.maximum(pair_forward, pair_reverse) + 1e-12
    )

    matched_percentiles, matched_bins = distance_matched_pair_percentiles(
        symmetric, pairs, min_separation=2
    )
    crossing_mask = crossing_pair_mask(normalized_pairs)
    canonical_mask = np.asarray(
        [sequence[left] + sequence[right] in CANONICAL for left, right in normalized_pairs],
        dtype=bool,
    )

    paired_position = np.zeros(length, dtype=bool)
    for left, right in normalized_pairs:
        paired_position[left] = True
        paired_position[right] = True
    off_diagonal_denominator = max(length - 1, 1)
    outbound = dependency.sum(axis=1, dtype=np.float64) / off_diagonal_denominator
    inbound = dependency.sum(axis=0, dtype=np.float64) / off_diagonal_denominator
    hub = 0.5 * (outbound + inbound)
    unpaired_position = ~paired_position

    broad, subtype, hierarchy = type_parts(str(record.get("rfam_type", "unknown")))
    row: Dict[str, object] = {
        "family_accession": record["family_accession"],
        "family_id": record.get("family_id", ""),
        "description": record.get("description", ""),
        "rfam_type": record.get("rfam_type", "unknown"),
        "broad_type": broad,
        "subtype": subtype,
        "type_hierarchy": hierarchy,
        "sequence_accession": record["sequence_accession"],
        "sequence": sequence,
        "length": length,
        "length_bin": bin_label(length, LENGTH_EDGES),
        "reference_paired_fraction": float(2 * len(normalized_pairs) / length),
        "gc_fraction": float(sum(base in "GC" for base in sequence) / length),
        "gc_ruggedness_11nt": gc_ruggedness(sequence),
        **standard,
        "distance_matched_pair_auc": float(np.mean(matched_percentiles)),
        "distance_matched_pair_percentile_median": float(np.median(matched_percentiles)),
        "distance_matched_pairs_evaluated": int(len(matched_percentiles)),
        "distance_matched_pair_fraction_above_0_99": float(np.mean(matched_percentiles >= 0.99)),
        "pair_direction_5to3_fraction": float(np.mean(pair_forward > pair_reverse)),
        "pair_direction_signed_mean": float(np.mean(direction_signed)),
        "pair_direction_asymmetry_median": float(np.median(np.abs(direction_signed))),
        "pair_reciprocity_median": float(np.median(reciprocity)),
        "canonical_reference_pairs": int(canonical_mask.sum()),
        "canonical_pair_score_median": float(np.median(pair_scores[canonical_mask])) if canonical_mask.any() else float("nan"),
        "noncanonical_pair_score_median": float(np.median(pair_scores[~canonical_mask])) if (~canonical_mask).any() else float("nan"),
        "crossing_reference_pairs": int(crossing_mask.sum()),
        "crossing_pair_score_median": float(np.median(pair_scores[crossing_mask])) if crossing_mask.any() else float("nan"),
        "noncrossing_pair_score_median": float(np.median(pair_scores[~crossing_mask])) if (~crossing_mask).any() else float("nan"),
        "paired_position_hub_median": float(np.median(hub[paired_position])) if paired_position.any() else float("nan"),
        "unpaired_position_hub_median": float(np.median(hub[unpaired_position])) if unpaired_position.any() else float("nan"),
        "paired_position_hub_enrichment": safe_ratio(
            float(np.median(hub[paired_position])) if paired_position.any() else float("nan"),
            float(np.median(hub[unpaired_position])) if unpaired_position.any() else float("nan"),
        ),
        "mean_outbound_dependency": float(np.mean(outbound)),
        "mean_inbound_dependency": float(np.mean(inbound)),
        "max_hub_position_1based": int(np.argmax(hub) + 1),
        "max_hub_score": float(np.max(hub)),
        "distance_matched_bins": {
            label: {
                "pairs": int(sum(item == label for item in matched_bins)),
                "mean_percentile": float(
                    np.mean([
                        value for value, item in zip(matched_percentiles, matched_bins) if item == label
                    ])
                ),
            }
            for label in sorted(set(matched_bins))
        },
    }
    return row


HEADLINE_METRICS = (
    "roc_auc",
    "average_precision",
    "precision_at_reference_pair_count",
    "pair_score_enrichment",
    "distance_matched_pair_auc",
    "distance_matched_pair_fraction_above_0_99",
    "pair_direction_asymmetry_median",
    "pair_reciprocity_median",
    "paired_position_hub_enrichment",
)


def summarize(rows: Sequence[Dict[str, object]], *, seed: int) -> Dict[str, object]:
    if not rows:
        return {"records": 0}
    output: Dict[str, object] = {
        "records": len(rows),
        "total_nucleotides": int(sum(int(row["length"]) for row in rows)),
        "total_reference_pairs": int(sum(int(row["reference_pairs"]) for row in rows)),
        "length_quantiles": quantiles(float(row["length"]) for row in rows),
    }
    for offset, metric in enumerate(HEADLINE_METRICS):
        values = [float(row[metric]) for row in rows if finite(row.get(metric))]
        if values:
            output[f"{metric}_quantiles"] = quantiles(values)
            output[f"{metric}_mean_ci95"] = bootstrap_mean_ci(
                values, seed=seed + 1009 * offset
            )
    output["families_above_random_distance_matched_auc"] = int(
        sum(float(row["distance_matched_pair_auc"]) > 0.5 for row in rows)
    )
    output["families_above_0_75_distance_matched_auc"] = int(
        sum(float(row["distance_matched_pair_auc"]) >= 0.75 for row in rows)
    )
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261026)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = [
        json.loads(line)
        for line in args.manifest.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    rows: List[Dict[str, object]] = []
    missing: List[str] = []
    failed: List[Dict[str, str]] = []
    for ordinal, record in enumerate(manifest, start=1):
        family = str(record["family_accession"])
        path = args.results / f"{family}.npz"
        if not path.is_file():
            missing.append(family)
            continue
        try:
            archive = np.load(path)
            sequence = str(archive["sequence"].item())
            expected = str(record["sequence"])
            if sequence != expected:
                raise ValueError("dependency sequence does not match manifest")
            dependency = np.asarray(archive["dependency"], dtype=np.float32)
            rows.append(analyze_record(record, dependency, sequence))
        except Exception as error:  # preserve a complete audit trail
            failed.append(
                {
                    "family_accession": family,
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            )
        if ordinal % 250 == 0:
            print(json.dumps({"processed": ordinal, "rows": len(rows), "failed": len(failed)}), flush=True)

    rows.sort(key=lambda row: str(row["family_accession"]))
    args.rows.parent.mkdir(parents=True, exist_ok=True)
    with args.rows.open("wt", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(json_safe(row), sort_keys=True, allow_nan=False) + "\n")

    groupings = {
        "by_broad_type": "broad_type",
        "by_subtype": "subtype",
        "by_type_hierarchy": "type_hierarchy",
        "by_length": "length_bin",
    }
    grouped_output: Dict[str, object] = {}
    for output_name, field in groupings.items():
        groups: Dict[str, List[Dict[str, object]]] = defaultdict(list)
        for row in rows:
            groups[str(row[field])].append(row)
        grouped_output[output_name] = {
            key: summarize(value, seed=args.seed + index * 7919)
            for index, (key, value) in enumerate(sorted(groups.items()))
            if len(value) >= 2
        }

    correlation_predictors = (
        "length",
        "gc_fraction",
        "gc_ruggedness_11nt",
        "reference_paired_fraction",
        "pair_reciprocity_median",
    )
    correlation_outcomes = (
        "roc_auc",
        "average_precision",
        "distance_matched_pair_auc",
    )
    correlations: Dict[str, Dict[str, float]] = {}
    for outcome in correlation_outcomes:
        correlations[outcome] = {}
        for predictor in correlation_predictors:
            correlations[outcome][predictor] = spearman(
                [float(row[predictor]) for row in rows],
                [float(row[outcome]) for row in rows],
            )

    eligible = [row for row in rows if int(row["reference_pairs"]) >= 5]
    ranked = sorted(
        eligible,
        key=lambda row: (float(row["distance_matched_pair_auc"]), str(row["family_accession"])),
        reverse=True,
    )
    ranking_fields = (
        "family_accession",
        "family_id",
        "description",
        "rfam_type",
        "sequence_accession",
        "length",
        "reference_pairs",
        "roc_auc",
        "average_precision",
        "distance_matched_pair_auc",
        "pair_score_enrichment",
    )

    def slim(row: Dict[str, object]) -> Dict[str, object]:
        return {field: row[field] for field in ranking_fields}

    payload = {
        "method": {
            "score": "max over 3 source substitutions and 4 target bases of absolute change in log2 odds",
            "orientation": "directed source-to-target; diagonal excluded",
            "reference": "Rfam consensus pairs mapped to one ungapped seed representative per family",
            "distance_matched_metric": "mean within-sequence percentile of every reference-pair score among non-pairs in the same power-of-two separation bin",
            "interpretation_boundary": "model dependency, not a demonstrated physical or molecular interaction",
        },
        "overall": summarize(rows, seed=args.seed),
        **grouped_output,
        "spearman_correlations": correlations,
        "top_30_families_by_distance_matched_auc": [slim(row) for row in ranked[:30]],
        "bottom_30_families_by_distance_matched_auc": [slim(row) for row in reversed(ranked[-30:])],
        "missing_count": len(missing),
        "missing_families": missing[:200],
        "failed": failed,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(
        json.dumps(json_safe(payload), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload["overall"], indent=2, sort_keys=True), flush=True)
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
