#!/usr/bin/env python3
"""Select one structure-covered Rfam seed representative per family.

The output is JSON Lines and includes the ungapped representative sequence,
its mapped consensus base pairs, Rfam metadata, and an approximate saturation
scan cost proportional to L^3.  Records are assigned to cost-balanced shards
with longest-processing-time-first scheduling.
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


BASES = set("ACGU")


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("rt", encoding="utf-8", errors="replace")


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


def iter_stockholm_families(path: Path) -> Iterator[Dict[str, object]]:
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
                continue
            if line.startswith("#=GF "):
                fields = line.split(None, 2)
                if len(fields) == 3 and fields[1] in {"AC", "ID", "DE", "TP", "CL"}:
                    metadata.setdefault(fields[1], fields[2].strip())
            elif line.startswith("#=GC SS_cons"):
                ss_cons += line.split(None, 2)[2].strip()
            elif line and not line.startswith("#"):
                fields = line.split()
                if len(fields) >= 2:
                    sequences[fields[0]] += fields[1]


def ungap_and_map(
    aligned: str,
    consensus_pairs: Sequence[Tuple[int, int]],
) -> Tuple[str, List[Tuple[int, int]], float]:
    column_to_position: Dict[int, int] = {}
    sequence: List[str] = []
    for column, symbol in enumerate(aligned):
        if symbol not in "-.~":
            column_to_position[column] = len(sequence)
            sequence.append(symbol.upper().replace("T", "U"))
    mapped = [
        (column_to_position[left], column_to_position[right])
        for left, right in consensus_pairs
        if left in column_to_position and right in column_to_position
    ]
    coverage = len(column_to_position) / max(len(aligned), 1)
    return "".join(sequence), mapped, coverage


def choose_representative(
    family: Dict[str, object],
    min_length: int,
    max_length: int,
) -> Optional[Dict[str, object]]:
    metadata = family["metadata"]
    ss_cons = str(family["ss_cons"])
    consensus_pairs = wuss_pairs(ss_cons)
    candidates = []
    for accession, aligned in family["sequences"].items():
        if len(aligned) != len(ss_cons):
            continue
        sequence, pairs, alignment_coverage = ungap_and_map(aligned, consensus_pairs)
        if not min_length <= len(sequence) <= max_length:
            continue
        if set(sequence).difference(BASES):
            continue
        structural_coverage = len(pairs) / max(len(consensus_pairs), 1)
        paired_fraction = 2 * len(pairs) / max(len(sequence), 1)
        # Prioritize structure coverage, then mapped pair count, then lower gap
        # content. Accession is the deterministic final tie breaker.
        score = (
            structural_coverage,
            len(pairs),
            alignment_coverage,
            -len(sequence),
            accession,
        )
        candidates.append((score, accession, sequence, pairs, alignment_coverage))
    if not candidates:
        return None
    _, accession, sequence, pairs, alignment_coverage = max(candidates)
    return {
        "family_accession": metadata.get("AC", ""),
        "family_id": metadata.get("ID", ""),
        "description": metadata.get("DE", ""),
        "rfam_type": metadata.get("TP", "unknown"),
        "clan": metadata.get("CL", ""),
        "sequence_accession": accession,
        "sequence": sequence,
        "length": len(sequence),
        "consensus_pairs": pairs,
        "reference_pair_count": len(pairs),
        "reference_paired_fraction": 2 * len(pairs) / len(sequence),
        "alignment_coverage": alignment_coverage,
        "estimated_cost": int(len(sequence) ** 3),
    }


def balanced_shards(records: List[Dict[str, object]], count: int) -> List[int]:
    loads = [0] * count
    assignments = [0] * len(records)
    order = sorted(
        range(len(records)),
        key=lambda index: (int(records[index]["estimated_cost"]), records[index]["family_accession"]),
        reverse=True,
    )
    for index in order:
        shard = min(range(count), key=lambda value: (loads[value], value))
        assignments[index] = shard
        loads[shard] += int(records[index]["estimated_cost"])
    return assignments


def stratified_subset(records: List[Dict[str, object]], limit: int) -> List[Dict[str, object]]:
    """Round-robin across Rfam types and length bins for a diverse pilot."""
    if limit <= 0 or limit >= len(records):
        return records
    bins: Dict[Tuple[str, int], List[Dict[str, object]]] = defaultdict(list)
    for record in records:
        length = int(record["length"])
        length_bin = 0 if length <= 80 else 1 if length <= 160 else 2 if length <= 320 else 3 if length <= 512 else 4
        type_path = str(record["rfam_type"]).strip() or "unknown"
        bins[(type_path, length_bin)].append(record)
    for values in bins.values():
        values.sort(key=lambda row: (-int(row["reference_pair_count"]), str(row["family_accession"])))
    keys = sorted(bins)
    selected: List[Dict[str, object]] = []
    while len(selected) < limit and keys:
        next_keys = []
        for key in keys:
            if bins[key] and len(selected) < limit:
                selected.append(bins[key].pop(0))
            if bins[key]:
                next_keys.append(key)
        keys = next_keys
    return selected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rfam-seed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--min-length", type=int, default=20)
    parser.add_argument("--max-length", type=int, default=1024)
    parser.add_argument("--min-reference-pairs", type=int, default=1)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--shards", type=int, default=4)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    records: List[Dict[str, object]] = []
    family_count = 0
    for family in iter_stockholm_families(args.rfam_seed):
        family_count += 1
        representative = choose_representative(
            family,
            min_length=args.min_length,
            max_length=args.max_length,
        )
        if representative is not None and int(representative["reference_pair_count"]) >= args.min_reference_pairs:
            records.append(representative)
    records.sort(key=lambda row: str(row["family_accession"]))
    records = stratified_subset(records, args.limit)
    assignments = balanced_shards(records, args.shards)
    for index, (record, shard) in enumerate(zip(records, assignments)):
        record["record_index"] = index
        record["shard"] = shard

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wt", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")

    type_counts = Counter(str(row["rfam_type"]).split(";")[0].strip() for row in records)
    shard_costs = Counter()
    shard_records = Counter()
    for row in records:
        shard = int(row["shard"])
        shard_costs[shard] += int(row["estimated_cost"])
        shard_records[shard] += 1
    lengths = sorted(int(row["length"]) for row in records)
    summary = {
        "rfam_families_seen": family_count,
        "selected_representatives": len(records),
        "length_min": min(lengths) if lengths else None,
        "length_median": lengths[len(lengths) // 2] if lengths else None,
        "length_max": max(lengths) if lengths else None,
        "total_nucleotides": sum(lengths),
        "total_forward_sequences": sum(1 + 3 * length for length in lengths),
        "broad_type_counts": dict(type_counts.most_common()),
        "shard_records": {str(key): shard_records[key] for key in range(args.shards)},
        "shard_costs": {str(key): shard_costs[key] for key in range(args.shards)},
        "filters": {
            "min_length": args.min_length,
            "max_length": args.max_length,
            "min_reference_pairs": args.min_reference_pairs,
            "limit": args.limit,
        },
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
