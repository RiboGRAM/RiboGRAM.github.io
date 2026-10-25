#!/usr/bin/env python3
"""Evaluate the released RiboGRAM zero-shot pairing atlas against Rfam Seed.

The website stores every Layer-12/Head-20 symmetrized, positive-APC pair with
score >= 0.4 for 92,781 Rfam 15.0 Seed sequences.  This program maps Rfam
SS_cons pairs onto each ungapped Seed sequence, evaluates those sparse calls,
aggregates calls in alignment-column space, and compares exact representative
sequences with the separately generated saturation-mutagenesis atlas.

The attention matrix used here is symmetric.  Its edges are candidate pairing
signals and must not be interpreted as directed perturbational dependencies.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple


Pair = Tuple[int, int]
GAPS = set("-.~")
CANONICAL_WOBBLE = {"AU", "UA", "GC", "CG", "GU", "UG"}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def iter_stockholm(path: Path) -> Iterator[Dict[str, object]]:
    metadata: Dict[str, str] = {}
    sequences: Dict[str, str] = defaultdict(str)
    ss_cons = ""
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "//":
                if metadata.get("AC"):
                    yield {
                        "metadata": dict(metadata),
                        "sequences": dict(sequences),
                        "ss_cons": ss_cons,
                    }
                metadata.clear()
                sequences.clear()
                ss_cons = ""
            elif line.startswith("#=GF "):
                fields = line.split(None, 2)
                if len(fields) == 3 and fields[1] in {"AC", "ID", "DE", "TP", "CL"}:
                    metadata.setdefault(fields[1], fields[2].strip())
            elif line.startswith("#=GC SS_cons"):
                ss_cons += line.split(None, 2)[2].strip()
            elif line and not line.startswith("#"):
                fields = line.split()
                if len(fields) >= 2:
                    sequences[fields[0]] += fields[1]


def wuss_pairs(structure: str) -> List[Pair]:
    opening_to_closing = {"(": ")", "[": "]", "{": "}", "<": ">"}
    closing_to_opening = {value: key for key, value in opening_to_closing.items()}
    stacks: Dict[str, List[int]] = {key: [] for key in opening_to_closing}
    letter_stacks: Dict[str, List[int]] = {}
    pairs: List[Pair] = []
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
    return sorted((min(left, right), max(left, right)) for left, right in pairs)


def map_alignment(aligned: str, consensus_pairs: Sequence[Pair]):
    column_to_position: Dict[int, int] = {}
    position_to_column: List[int] = []
    sequence: List[str] = []
    for column, symbol in enumerate(aligned):
        if symbol not in GAPS:
            column_to_position[column] = len(sequence)
            position_to_column.append(column)
            sequence.append(symbol.upper().replace("T", "U"))
    mapped = {
        (column_to_position[left], column_to_position[right])
        for left, right in consensus_pairs
        if left in column_to_position and right in column_to_position
    }
    return "".join(sequence), mapped, position_to_column


def safe_div(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator else 0.0


def prf(true_positive: int, predicted: int, reference: int) -> Dict[str, float]:
    precision = safe_div(true_positive, predicted)
    recall = safe_div(true_positive, reference)
    return {
        "true_positive": int(true_positive),
        "predicted": int(predicted),
        "reference": int(reference),
        "precision": precision,
        "recall": recall,
        "f1": safe_div(2.0 * precision * recall, precision + recall),
    }


def quantiles(values: Sequence[float]) -> Dict[str, float]:
    ordered = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not ordered:
        return {}
    result = {}
    for probability in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0):
        index = probability * (len(ordered) - 1)
        lower = int(math.floor(index))
        upper = int(math.ceil(index))
        fraction = index - lower
        result[str(probability)] = ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction
    return result


def stem_pairs(pairs: Set[Pair], minimum_run: int) -> Set[Pair]:
    selected: Set[Pair] = set()
    for left, right in pairs:
        if (left - 1, right + 1) in pairs:
            continue
        run: List[Pair] = []
        cursor = (left, right)
        while cursor in pairs:
            run.append(cursor)
            cursor = (cursor[0] + 1, cursor[1] - 1)
        if len(run) >= minimum_run:
            selected.update(run)
    return selected


def family_metrics(rows: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    true_positive = sum(int(row["true_positive"]) for row in rows)
    predicted = sum(int(row["predicted_pairs"]) for row in rows)
    reference = sum(int(row["reference_pairs"]) for row in rows)
    sequence_f1 = [float(row["f1"]) for row in rows]
    return {
        "sequences": len(rows),
        "micro": prf(true_positive, predicted, reference),
        "sequence_f1_quantiles": quantiles(sequence_f1),
        "pair_bearing_fraction": safe_div(sum(int(row["predicted_pairs"]) > 0 for row in rows), len(rows)),
        "predicted_pairs_per_sequence_quantiles": quantiles([int(row["predicted_pairs"]) for row in rows]),
    }


def load_json_gzip(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rfam-seed", type=Path, required=True)
    parser.add_argument("--details-dir", type=Path, required=True)
    parser.add_argument("--dependency-dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--family-rows", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sequence_rows: List[Dict[str, object]] = []
    family_rows: List[Dict[str, object]] = []
    by_broad_type: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    by_subtype: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    by_length: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    global_predicted_stem2: Set[Tuple[str, str, Pair]] = set()
    global_predicted_stem3: Set[Tuple[str, str, Pair]] = set()
    stem_counts = Counter()
    base_counts = Counter()
    missing_detail_families: List[str] = []
    missing_detail_sequences: List[str] = []
    representative_lookup: Dict[Tuple[str, str], Dict[str, object]] = {}
    alignment_consensus_summaries: Dict[str, Dict[str, int]] = {
        str(value): Counter() for value in (0.05, 0.10, 0.25, 0.50)
    }

    families_seen = 0
    for family in iter_stockholm(args.rfam_seed):
        families_seen += 1
        metadata = family["metadata"]
        family_accession = str(metadata["AC"])
        detail_path = args.details_dir / f"{family_accession}.json.gz"
        if not detail_path.is_file():
            missing_detail_families.append(family_accession)
            continue
        detail = load_json_gzip(detail_path)
        point_by_accession = {
            str(point.get("source_accession", "")): point
            for point in detail.get("points", {}).values()
        }
        consensus_pairs = wuss_pairs(str(family["ss_cons"]))
        predicted_column_support: Counter[Pair] = Counter()
        coverage_columns: List[Set[int]] = []
        family_sequence_rows: List[Dict[str, object]] = []
        aligned_by_accession = dict(family["sequences"])

        type_path = [part.strip() for part in str(metadata.get("TP", "unknown")).split(";") if part.strip()]
        broad_type = type_path[0] if type_path else "unknown"
        subtype = type_path[-1] if type_path else "unknown"

        for accession, aligned in aligned_by_accession.items():
            point = point_by_accession.get(accession)
            if point is None:
                missing_detail_sequences.append(f"{family_accession}:{accession}")
                continue
            sequence, reference_pairs, position_to_column = map_alignment(aligned, consensus_pairs)
            site_sequence = str(point.get("sequence", "")).upper().replace("T", "U")
            if sequence != site_sequence:
                raise ValueError(f"sequence mismatch for {family_accession}:{accession}")
            predicted_pairs = {
                (min(int(pair[0]), int(pair[1])) - 1, max(int(pair[0]), int(pair[1])) - 1)
                for pair in point.get("attention", {}).get("pairs", [])
                if int(pair[0]) != int(pair[1])
            }
            if any(left < 0 or right >= len(sequence) for left, right in predicted_pairs):
                raise ValueError(f"pair indexing mismatch for {family_accession}:{accession}")
            true_positive_pairs = predicted_pairs & reference_pairs
            precision = safe_div(len(true_positive_pairs), len(predicted_pairs))
            recall = safe_div(len(true_positive_pairs), len(reference_pairs))
            f1 = safe_div(2.0 * precision * recall, precision + recall)
            length = len(sequence)
            length_bin = (
                "20-80" if length <= 80 else
                "81-160" if length <= 160 else
                "161-320" if length <= 320 else
                "321-512" if length <= 512 else "513-1024"
            )
            predicted_stem2 = stem_pairs(predicted_pairs, 2)
            predicted_stem3 = stem_pairs(predicted_pairs, 3)
            reference_stem2 = stem_pairs(reference_pairs, 2)
            reference_stem3 = stem_pairs(reference_pairs, 3)
            row = {
                "family_accession": family_accession,
                "family_id": metadata.get("ID", ""),
                "sequence_accession": accession,
                "length": length,
                "broad_type": broad_type,
                "subtype": subtype,
                "predicted_pairs": len(predicted_pairs),
                "reference_pairs": len(reference_pairs),
                "true_positive": len(true_positive_pairs),
                "precision": precision,
                "recall": recall,
                "f1": f1,
                "predicted_stem2_pairs": len(predicted_stem2),
                "predicted_stem3_pairs": len(predicted_stem3),
                "reference_stem2_pairs": len(reference_stem2),
                "reference_stem3_pairs": len(reference_stem3),
                "true_positive_predicted_stem2": len(true_positive_pairs & predicted_stem2),
                "true_positive_predicted_stem3": len(true_positive_pairs & predicted_stem3),
                "reference_stem2_recovered": len(reference_stem2 & predicted_pairs),
                "reference_stem3_recovered": len(reference_stem3 & predicted_pairs),
            }
            sequence_rows.append(row)
            family_sequence_rows.append(row)
            by_broad_type[broad_type].append(row)
            by_subtype[subtype].append(row)
            by_length[length_bin].append(row)
            representative_lookup[(family_accession, accession)] = {
                "predicted_pairs": predicted_pairs,
                "reference_pairs": reference_pairs,
            }

            stem_counts["predicted"] += len(predicted_pairs)
            stem_counts["predicted_stem2"] += len(predicted_stem2)
            stem_counts["predicted_stem3"] += len(predicted_stem3)
            stem_counts["tp"] += len(true_positive_pairs)
            stem_counts["tp_stem2"] += len(true_positive_pairs & predicted_stem2)
            stem_counts["tp_stem3"] += len(true_positive_pairs & predicted_stem3)
            stem_counts["reference"] += len(reference_pairs)
            stem_counts["reference_stem2"] += len(reference_stem2)
            stem_counts["reference_stem3"] += len(reference_stem3)
            stem_counts["reference_stem2_recovered"] += len(reference_stem2 & predicted_pairs)
            stem_counts["reference_stem3_recovered"] += len(reference_stem3 & predicted_pairs)

            for left, right in predicted_pairs:
                bases = sequence[left] + sequence[right]
                base_counts["predicted"] += 1
                base_counts["predicted_canonical_wobble"] += bases in CANONICAL_WOBBLE
                if (left, right) in reference_pairs:
                    base_counts["tp"] += 1
                    base_counts["tp_canonical_wobble"] += bases in CANONICAL_WOBBLE
                column_pair = tuple(sorted((position_to_column[left], position_to_column[right])))
                predicted_column_support[column_pair] += 1
            coverage_columns.append(set(position_to_column))

        candidate_columns = set(predicted_column_support) | set(consensus_pairs)
        coverage: Counter[Pair] = Counter()
        for column_pair in candidate_columns:
            left, right = column_pair
            coverage[column_pair] = sum(left in columns and right in columns for columns in coverage_columns)

        family_consensus_results = {}
        for threshold in (0.05, 0.10, 0.25, 0.50):
            predicted_consensus = {
                pair for pair, count in predicted_column_support.items()
                if safe_div(count, coverage[pair]) >= threshold
            }
            reference_consensus = set(consensus_pairs)
            tp = len(predicted_consensus & reference_consensus)
            result = prf(tp, len(predicted_consensus), len(reference_consensus))
            family_consensus_results[str(threshold)] = result
            accumulator = alignment_consensus_summaries[str(threshold)]
            accumulator["true_positive"] += tp
            accumulator["predicted"] += len(predicted_consensus)
            accumulator["reference"] += len(reference_consensus)

        aggregate = family_metrics(family_sequence_rows)
        family_row = {
            "family_accession": family_accession,
            "family_id": metadata.get("ID", ""),
            "broad_type": broad_type,
            "subtype": subtype,
            **aggregate,
            "alignment_consensus": family_consensus_results,
        }
        family_rows.append(family_row)

    dependency_overlap = None
    if args.dependency_dir and args.dependency_dir.is_dir():
        counts = Counter()
        families = 0
        for dependency_path in sorted(args.dependency_dir.glob("*.json.gz")):
            dependency = load_json_gzip(dependency_path)
            key = (str(dependency["family_accession"]), str(dependency["sequence_accession"]))
            zero = representative_lookup.get(key)
            if zero is None:
                continue
            families += 1
            predicted_pairs = set(zero["predicted_pairs"])
            reference_pairs = set(zero["reference_pairs"])
            q999_displayed = {
                (int(arc[0]), int(arc[1]))
                for arc in dependency.get("arcs", []) if int(arc[4]) == 2
            }
            zero_tp = predicted_pairs & reference_pairs
            dep_tp_displayed = q999_displayed & reference_pairs
            counts["reference"] += len(reference_pairs)
            counts["zero_predicted"] += len(predicted_pairs)
            counts["zero_tp"] += len(zero_tp)
            counts["dep_q999_displayed"] += len(q999_displayed)
            counts["dep_tp_displayed"] += len(dep_tp_displayed)
            counts["both_reference"] += len(reference_pairs & predicted_pairs & q999_displayed)
            counts["either_reference_displayed"] += len(reference_pairs & (predicted_pairs | q999_displayed))
            counts["zero_dep_overlap_displayed"] += len(predicted_pairs & q999_displayed)
            summary = dependency.get("summary", {})
            counts["dep_q999_full"] += int(summary.get("q999_pairs", 0))
            counts["dep_reference_q999_full"] += int(summary.get("reference_pairs_q999", 0))

        dependency_overlap = {
            "families_with_exact_representative": families,
            "zero_shot": prf(counts["zero_tp"], counts["zero_predicted"], counts["reference"]),
            "dependency_q999_full": prf(
                counts["dep_reference_q999_full"], counts["dep_q999_full"], counts["reference"]
            ),
            "displayed_arc_overlap_lower_bound": {
                "note": "Atlas dependency arcs are capped at 250 per family; positional overlap values are lower bounds.",
                "zero_shot_and_dependency_pairs": counts["zero_dep_overlap_displayed"],
                "zero_shot_pairs_supported_by_dependency_fraction": safe_div(
                    counts["zero_dep_overlap_displayed"], counts["zero_predicted"]
                ),
                "reference_pairs_recovered_by_both": counts["both_reference"],
                "reference_pairs_recovered_by_either": counts["either_reference_displayed"],
                "reference_recall_both": safe_div(counts["both_reference"], counts["reference"]),
                "reference_recall_either": safe_div(counts["either_reference_displayed"], counts["reference"]),
            },
        }

    consensus_summary = {}
    for threshold, values in alignment_consensus_summaries.items():
        consensus_summary[threshold] = prf(
            values["true_positive"], values["predicted"], values["reference"]
        )

    global_metrics = family_metrics(sequence_rows)
    family_f1 = [float(row["micro"]["f1"]) for row in family_rows]
    summary = {
        "method": {
            "zero_shot_signal": "RiboGRAM Layer-12 Head-20 attention, symmetrized, positive APC, score >= 0.4",
            "reference": "Rfam 15.0 SS_cons pairs mapped to each ungapped Seed sequence",
            "rfam_seed_url": "https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/Rfam.seed.gz",
            "rfam_seed_sha256": sha256(args.rfam_seed),
            "directionality": "none; the stored matrix is explicitly symmetrized",
            "interpretation_boundary": "candidate model pairing signal, not an experimental contact or perturbational dependency",
        },
        "coverage": {
            "rfam_families_seen": families_seen,
            "families_evaluated": len(family_rows),
            "sequences_evaluated": len(sequence_rows),
            "missing_detail_families": missing_detail_families,
            "missing_detail_sequences_count": len(missing_detail_sequences),
            "missing_detail_sequences_sample": missing_detail_sequences[:20],
        },
        "overall": global_metrics,
        "family_f1_quantiles": quantiles(family_f1),
        "by_broad_type": {key: family_metrics(value) for key, value in sorted(by_broad_type.items())},
        "by_subtype": {key: family_metrics(value) for key, value in sorted(by_subtype.items())},
        "by_length": {key: family_metrics(value) for key, value in sorted(by_length.items())},
        "geometry": {
            "predicted_pairs_in_antidiagonal_run_ge_2_fraction": safe_div(stem_counts["predicted_stem2"], stem_counts["predicted"]),
            "predicted_pairs_in_antidiagonal_run_ge_3_fraction": safe_div(stem_counts["predicted_stem3"], stem_counts["predicted"]),
            "precision_all_predicted_pairs": safe_div(stem_counts["tp"], stem_counts["predicted"]),
            "precision_pairs_in_run_ge_2": safe_div(stem_counts["tp_stem2"], stem_counts["predicted_stem2"]),
            "precision_pairs_in_run_ge_3": safe_div(stem_counts["tp_stem3"], stem_counts["predicted_stem3"]),
            "reference_run_ge_2_recall": safe_div(stem_counts["reference_stem2_recovered"], stem_counts["reference_stem2"]),
            "reference_run_ge_3_recall": safe_div(stem_counts["reference_stem3_recovered"], stem_counts["reference_stem3"]),
        },
        "base_pair_composition": {
            "predicted_canonical_or_wobble_fraction": safe_div(base_counts["predicted_canonical_wobble"], base_counts["predicted"]),
            "true_positive_canonical_or_wobble_fraction": safe_div(base_counts["tp_canonical_wobble"], base_counts["tp"]),
        },
        "family_alignment_consensus": consensus_summary,
        "dependency_overlap": dependency_overlap,
        "limitations": [
            "Rfam Seed sequences are not an independent held-out set and may overlap the pretraining corpus.",
            "Only attention cells at or above 0.4 are stored; score-ranking AUC cannot be reconstructed.",
            "The zero-shot matrix is symmetrized, so it cannot support source-to-target direction claims.",
            "Displayed dependency arcs are capped at 250 per family, so positional overlap with dependency q99.9 edges is a lower bound.",
        ],
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.family_rows.parent.mkdir(parents=True, exist_ok=True)
    with args.family_rows.open("wt", encoding="utf-8") as handle:
        for row in family_rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
