#!/usr/bin/env python3
"""Run one cost-balanced shard of Rfam representative dependency scans."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import traceback
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import torch

from ribogram_dependency_scan import compute_dependency_map, normalize_rna


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_") or "record"


def atomic_save_npz(path: Path, **arrays) -> None:
    temporary = path.with_name(path.name + ".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    os.replace(temporary, path)


def atomic_write_json(path: Path, payload: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def load_manifest(path: Path, shard: int) -> List[Dict[str, object]]:
    records = []
    with path.open("rt", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            record = json.loads(raw)
            if int(record["shard"]) == shard:
                records.append(record)
    records.sort(key=lambda row: (-int(row["estimated_cost"]), str(row["family_accession"])))
    return records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--max-batch-size", type=int, default=192)
    parser.add_argument("--attention-element-budget", type=int, default=5_000_000)
    parser.add_argument("--eps", type=float, default=1e-6)
    parser.add_argument("--max-records", type=int, default=0)
    parser.add_argument("--fail-fast", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    records = load_manifest(args.manifest, args.shard)
    if args.max_records > 0:
        records = records[: args.max_records]
    if not records:
        raise ValueError(f"manifest contains no records for shard {args.shard}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.output_dir / f"worker_{args.shard:02d}.jsonl"
    progress_path = args.output_dir / f"progress_{args.shard:02d}.json"
    success_path = args.output_dir / f"SHARD_{args.shard:02d}_SUCCESS"

    sys.path.insert(0, str(args.source_dir.resolve()))
    from ribogram import Architecture, load_model_and_alphabet

    device = torch.device(args.device)
    if device.type == "cuda":
        torch.cuda.set_device(device)
    alphabet, converter, model = load_model_and_alphabet(
        args.checkpoint,
        architecture=Architecture.RIBOGRAM,
    )
    model = model.to(device).eval()

    completed = 0
    skipped = 0
    failed = 0
    started = time.time()
    for ordinal, record in enumerate(records, start=1):
        family = safe_name(str(record["family_accession"]))
        output_path = args.output_dir / f"{family}.npz"
        metadata_path = args.output_dir / f"{family}.json"
        if output_path.is_file() and metadata_path.is_file():
            skipped += 1
            continue
        event = {
            "event": "scan",
            "family_accession": record["family_accession"],
            "sequence_accession": record["sequence_accession"],
            "length": record["length"],
            "ordinal": ordinal,
            "records_in_shard": len(records),
            "started_at": time.strftime("%FT%T%z"),
        }
        try:
            sequence = normalize_rna(str(record["sequence"]))
            last_error = None
            result = None
            # Retry OOMs with progressively smaller B*L^2 budgets.
            for divisor in (1, 2, 4, 8):
                budget = max(3, args.attention_element_budget // divisor)
                try:
                    result = compute_dependency_map(
                        model,
                        alphabet,
                        converter,
                        str(record["sequence_accession"]),
                        sequence,
                        device=device,
                        max_batch_size=args.max_batch_size,
                        attention_element_budget=budget,
                        eps=args.eps,
                    )
                    event["attention_element_budget_used"] = budget
                    break
                except torch.cuda.OutOfMemoryError as error:
                    last_error = error
                    torch.cuda.empty_cache()
            if result is None:
                raise last_error or RuntimeError("dependency computation failed")
            dependency, variant_influence, stats = result
            atomic_save_npz(
                output_path,
                sequence_id=np.asarray(str(record["sequence_accession"])),
                sequence=np.asarray(sequence),
                dependency=dependency.astype(np.float16),
                variant_influence=variant_influence.astype(np.float16),
            )
            payload = {
                "manifest": record,
                "scan": asdict(stats),
                "worker": {
                    "shard": args.shard,
                    "device": str(device),
                    "attention_element_budget_used": event["attention_element_budget_used"],
                },
            }
            atomic_write_json(metadata_path, payload)
            event.update(
                status="complete",
                elapsed_seconds=stats.elapsed_seconds,
                dependency_max=stats.dependency_max,
            )
            completed += 1
        except Exception as error:
            failed += 1
            event.update(
                status="failed",
                error_type=type(error).__name__,
                error=str(error),
                traceback=traceback.format_exc(limit=8),
            )
            if args.fail_fast:
                with log_path.open("at", encoding="utf-8") as handle:
                    handle.write(json.dumps(event, sort_keys=True) + "\n")
                raise
        with log_path.open("at", encoding="utf-8") as handle:
            handle.write(json.dumps(event, sort_keys=True) + "\n")
        atomic_write_json(
            progress_path,
            {
                "shard": args.shard,
                "records": len(records),
                "completed_this_run": completed,
                "skipped_existing": skipped,
                "failed": failed,
                "processed": ordinal,
                "elapsed_seconds": time.time() - started,
                "updated_at": time.strftime("%FT%T%z"),
            },
        )
        print(json.dumps(event, sort_keys=True), flush=True)

    final = {
        "shard": args.shard,
        "records": len(records),
        "completed_this_run": completed,
        "skipped_existing": skipped,
        "failed": failed,
        "elapsed_seconds": time.time() - started,
        "finished_at": time.strftime("%FT%T%z"),
    }
    atomic_write_json(success_path, final)
    print(json.dumps(final, indent=2, sort_keys=True), flush=True)
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
