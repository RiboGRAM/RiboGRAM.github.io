#!/usr/bin/env python3
"""Fail-fast integrity checks for the RNAcentral atlas artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    root = args.results
    merged = json.loads((root / "merged_manifest.json").read_text(encoding="utf-8"))
    n, nt = int(merged["records"]), int(merged["nucleotides"])
    shapes = {
        "embeddings": (n, 1280), "ids": (n,), "lengths": (n,),
        "offsets": (n + 1,), "attention_tracks": (nt, 3), "pair_partner": (nt,),
        "ncrna13_id": (n,), "ncrna88_id": (n,),
    }
    for name, expected in shapes.items():
        observed = np.load(root / f"{name}.npy", mmap_mode="r").shape
        if observed != expected:
            raise RuntimeError(f"{name}: {observed} != {expected}")
    offsets = np.load(root / "offsets.npy", mmap_mode="r")
    if int(offsets[0]) != 0 or int(offsets[-1]) != nt:
        raise RuntimeError("offset endpoints are invalid")
    if np.any(np.diff(offsets) <= 0):
        raise RuntimeError("offsets are not strictly increasing")
    lengths = np.load(root / "lengths.npy", mmap_mode="r")
    if not np.array_equal(np.diff(offsets), lengths.astype(np.uint64)):
        raise RuntimeError("lengths disagree with nucleotide offsets")
    ids = np.load(root / "ids.npy", mmap_mode="r")
    if len(np.unique(ids)) != n:
        raise RuntimeError("RNAcentral IDs are not unique")

    rng = np.random.default_rng(20261025)
    sample = rng.choice(n, size=min(n, 10000), replace=False)
    embeddings = np.load(root / "embeddings.npy", mmap_mode="r")
    if not np.isfinite(np.asarray(embeddings[sample], dtype=np.float32)).all():
        raise RuntimeError("non-finite embedding values")
    tracks = np.load(root / "attention_tracks.npy", mmap_mode="r")
    nt_sample = rng.choice(nt, size=min(nt, 100000), replace=False)
    if not np.isfinite(np.asarray(tracks[nt_sample], dtype=np.float32)).all():
        raise RuntimeError("non-finite attention values")

    report = {"records": n, "nucleotides": nt}
    if (root / "KNN_SUCCESS").exists():
        knn = np.load(root / "knn_indices.npy", mmap_mode="r")
        similarity = np.load(root / "knn_similarity.npy", mmap_mode="r")
        if knn.shape[0] != n or similarity.shape != knn.shape:
            raise RuntimeError("kNN arrays have invalid shapes")
        report["knn_invalid_fraction"] = float(np.mean(knn[sample] < 0))
        report["knn_finite_fraction"] = float(np.mean(np.isfinite(similarity[sample])))
    if (root / "CLUSTER_SUCCESS").exists():
        labels = np.load(root / "cluster_id.npy", mmap_mode="r")
        sizes = np.load(root / "cluster_size.npy")
        if len(labels) != n or int(sizes.sum()) != n:
            raise RuntimeError("cluster membership totals are invalid")
        status = np.load(root / "annotation_status.npy", mmap_mode="r")
        if len(status) != n or int(status.max()) > 6:
            raise RuntimeError("annotation status values are invalid")
        report["clusters"] = int(len(sizes))
        report["status_counts"] = np.bincount(status, minlength=7).tolist()
    if (root / "DISPLAY_SUCCESS").exists():
        display = np.load(root / "display_indices.npy")
        xy = np.load(root / "display_xy.npy")
        if xy.shape != (len(display), 2) or not np.isfinite(xy).all():
            raise RuntimeError("display projection is invalid")
        if len(np.unique(display)) != len(display):
            raise RuntimeError("display indices are not unique")
        report["display_points"] = int(len(display))
    (root / "validation_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
