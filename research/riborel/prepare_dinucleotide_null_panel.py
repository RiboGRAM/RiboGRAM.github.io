#!/usr/bin/env python3
"""Create a composition-matched null panel for dependency-map calibration.

Each selected Rfam representative is replaced by an exact dinucleotide-count
preserving Eulerian shuffle.  The default allocation intentionally limits the
number of long RNAs because saturation-mutagenesis cost grows approximately as
L^3 while still yielding millions of long-distance null edges.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


LENGTH_TARGETS = (
    (20, 80, 40),
    (81, 160, 40),
    (161, 320, 32),
    (321, 512, 16),
    (513, 1024, 8),
)
GC_EDGES = (0.0, 0.35, 0.50, 0.65, 1.0000001)


def stable_seed(text: str, seed: int) -> int:
    digest = hashlib.sha256(f"{seed}\0{text}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "little")


def dinucleotides(sequence: str) -> Counter:
    return Counter(sequence[index : index + 2] for index in range(len(sequence) - 1))


def dinucleotide_shuffle(sequence: str, *, seed: int) -> str:
    """Return an exact directed-dinucleotide-preserving Eulerian shuffle."""
    if len(sequence) < 3:
        return sequence
    expected = dinucleotides(sequence)
    best = sequence
    for attempt in range(32):
        rng = np.random.default_rng(seed + attempt * 104729)
        adjacency: Dict[str, List[str]] = defaultdict(list)
        for left, right in zip(sequence[:-1], sequence[1:]):
            adjacency[left].append(right)
        for edges in adjacency.values():
            rng.shuffle(edges)

        stack = [sequence[0]]
        reverse_path: List[str] = []
        while stack:
            current = stack[-1]
            if adjacency[current]:
                stack.append(adjacency[current].pop())
            else:
                reverse_path.append(stack.pop())
        shuffled = "".join(reversed(reverse_path))
        if len(shuffled) != len(sequence) or dinucleotides(shuffled) != expected:
            raise RuntimeError("Eulerian shuffle failed its dinucleotide invariant")
        best = shuffled
        if shuffled != sequence:
            break
    return best


def length_bin(length: int) -> str:
    for lower, upper, _ in LENGTH_TARGETS:
        if lower <= length <= upper:
            return f"{lower}-{upper}"
    return "out-of-range"


def gc_bin(sequence: str) -> str:
    fraction = sum(base in "GC" for base in sequence) / len(sequence)
    for lower, upper in zip(GC_EDGES[:-1], GC_EDGES[1:]):
        if lower <= fraction < upper:
            return f"{lower:.2f}-{min(upper, 1.0):.2f}"
    return "unknown"


def broad_type(rfam_type: str) -> str:
    return rfam_type.split(";", 1)[0].strip() or "unknown"


def choose_stratified(records: Sequence[Dict[str, object]], *, seed: int) -> List[Dict[str, object]]:
    selected: List[Dict[str, object]] = []
    for lower, upper, target in LENGTH_TARGETS:
        candidates = [
            record for record in records if lower <= int(record["length"]) <= upper
        ]
        strata: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
        for record in candidates:
            sequence = str(record["sequence"])
            key = (broad_type(str(record.get("rfam_type", "unknown"))), gc_bin(sequence))
            strata[key].append(record)
        for values in strata.values():
            values.sort(
                key=lambda record: stable_seed(str(record["family_accession"]), seed)
            )
        keys = sorted(strata)
        taken = 0
        while keys and taken < target:
            next_keys: List[Tuple[str, str]] = []
            for key in keys:
                if strata[key] and taken < target:
                    selected.append(strata[key].pop(0))
                    taken += 1
                if strata[key]:
                    next_keys.append(key)
            keys = next_keys
    return selected


def balanced_shards(records: Sequence[Dict[str, object]], count: int) -> List[int]:
    loads = [0] * count
    assignments = [0] * len(records)
    order = sorted(
        range(len(records)),
        key=lambda index: (
            int(records[index]["estimated_cost"]),
            str(records[index]["family_accession"]),
        ),
        reverse=True,
    )
    for index in order:
        shard = min(range(count), key=lambda candidate: (loads[candidate], candidate))
        assignments[index] = shard
        loads[shard] += int(records[index]["estimated_cost"])
    return assignments


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261026)
    parser.add_argument("--shards", type=int, default=4)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source = [
        json.loads(line)
        for line in args.manifest.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    chosen = choose_stratified(source, seed=args.seed)
    records: List[Dict[str, object]] = []
    unchanged = 0
    for record in chosen:
        original_family = str(record["family_accession"])
        original_sequence = str(record["sequence"])
        shuffled = dinucleotide_shuffle(
            original_sequence,
            seed=stable_seed(original_family, args.seed),
        )
        unchanged += int(shuffled == original_sequence)
        null_record = dict(record)
        null_record.update(
            family_accession=f"NULL_{original_family}",
            family_id=f"dinucleotide-null:{record.get('family_id', original_family)}",
            description=f"Dinucleotide-preserving null for {original_family}",
            sequence_accession=f"NULL:{record['sequence_accession']}",
            sequence=shuffled,
            consensus_pairs=[],
            reference_pair_count=0,
            reference_paired_fraction=0.0,
            estimated_cost=int(len(shuffled) ** 3),
            null_source_family=original_family,
            null_source_accession=record["sequence_accession"],
            null_method="exact directed-dinucleotide-preserving Eulerian shuffle",
            null_seed=args.seed,
            gc_bin=gc_bin(shuffled),
            length_bin=length_bin(len(shuffled)),
        )
        records.append(null_record)

    assignments = balanced_shards(records, args.shards)
    for index, (record, shard) in enumerate(zip(records, assignments)):
        record["record_index"] = index
        record["shard"] = shard

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wt", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")

    shard_cost = Counter()
    shard_records = Counter()
    for record in records:
        shard_cost[int(record["shard"])] += int(record["estimated_cost"])
        shard_records[int(record["shard"])] += 1
    payload = {
        "source_records": len(source),
        "selected_nulls": len(records),
        "unchanged_forced_sequences": unchanged,
        "selection_targets": [
            {"length_min": lower, "length_max": upper, "target": target}
            for lower, upper, target in LENGTH_TARGETS
        ],
        "selected_by_length": dict(Counter(str(record["length_bin"]) for record in records)),
        "selected_by_broad_type": dict(
            Counter(broad_type(str(record.get("rfam_type", "unknown"))) for record in records)
        ),
        "selected_by_gc": dict(Counter(str(record["gc_bin"]) for record in records)),
        "total_nucleotides": sum(int(record["length"]) for record in records),
        "total_forward_sequences": sum(1 + 3 * int(record["length"]) for record in records),
        "shard_records": {str(index): shard_records[index] for index in range(args.shards)},
        "shard_costs": {str(index): shard_cost[index] for index in range(args.shards)},
        "seed": args.seed,
        "null_method": "exact directed-dinucleotide-preserving Eulerian shuffle",
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
