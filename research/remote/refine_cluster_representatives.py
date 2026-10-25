#!/usr/bin/env python3
"""Replace arbitrary Leiden representatives with PCA-centroid medoids."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy import sparse


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--batch", type=int, default=200000)
    args = parser.parse_args()
    root = args.results
    labels = np.load(root / "cluster_id.npy", mmap_mode="r")
    sizes = np.load(root / "cluster_size.npy")
    manifest = json.loads((root / "knn_manifest.json").read_text(encoding="utf-8"))
    values = np.load(root / f"pca{manifest['pca_dimension']}.npy", mmap_mode="r")
    n, dimension = values.shape
    clusters = len(sizes)
    started = time.time()

    rows = np.asarray(labels, dtype=np.int32)
    columns = np.arange(n, dtype=np.int32)
    membership = sparse.csr_matrix(
        (np.ones(n, dtype=np.float32), (rows, columns)), shape=(clusters, n)
    )
    centroids = np.asarray(membership @ values, dtype=np.float32)
    centroids /= np.maximum(sizes[:, None], 1)
    centroids /= np.maximum(np.linalg.norm(centroids, axis=1, keepdims=True), 1e-12)
    del membership, columns, rows

    best_score = np.full(clusters, -np.inf, dtype=np.float32)
    for start in range(0, n, args.batch):
        stop = min(start + args.batch, n)
        batch_labels = np.asarray(labels[start:stop], dtype=np.int32)
        score = np.einsum(
            "ij,ij->i",
            np.asarray(values[start:stop], dtype=np.float32),
            centroids[batch_labels],
            optimize=True,
        )
        np.maximum.at(best_score, batch_labels, score)

    representative = np.full(clusters, np.iinfo(np.int32).max, dtype=np.int32)
    tolerance = 1e-6
    for start in range(0, n, args.batch):
        stop = min(start + args.batch, n)
        batch_labels = np.asarray(labels[start:stop], dtype=np.int32)
        score = np.einsum(
            "ij,ij->i",
            np.asarray(values[start:stop], dtype=np.float32),
            centroids[batch_labels],
            optimize=True,
        )
        keep = score >= best_score[batch_labels] - tolerance
        indices = np.arange(start, stop, dtype=np.int32)[keep]
        np.minimum.at(representative, batch_labels[keep], indices)
    if np.any(representative == np.iinfo(np.int32).max):
        raise RuntimeError("at least one cluster has no representative")
    np.save(root / "cluster_representative.npy", representative)
    np.save(root / "cluster_centroid_pca.npy", centroids.astype(np.float16))
    report = {
        "method": "member with maximum cosine similarity to the normalized PCA-64 cluster centroid",
        "clusters": clusters,
        "mean_representative_similarity": float(best_score.mean()),
        "minimum_representative_similarity": float(best_score.min()),
        "elapsed_seconds": time.time() - started,
    }
    (root / "representative_manifest.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
