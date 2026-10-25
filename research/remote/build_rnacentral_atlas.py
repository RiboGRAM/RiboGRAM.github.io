#!/usr/bin/env python3
"""Build the RNAcentral–RiboGRAM pilot atlas in resumable phases.

The pipeline deliberately separates raw model outputs, direct database
annotations, graph clustering, and propagated labels.  A propagated label is
never written as if it were an RNAcentral/Rfam annotation.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import pickle
import shutil
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Sequence, Tuple

import numpy as np
from numpy.lib.format import open_memmap


SEED = 20261025
RECORD_ARRAYS = (
    "embeddings",
    "ids",
    "rna_types",
    "lengths",
    "taxa_counts",
    "ambiguous_counts",
    "ncrna13_id",
    "ncrna13_confidence",
    "ncrna88_id",
    "ncrna88_confidence",
)
NT_ARRAYS = ("attention_tracks", "pair_partner")


def log(message: str) -> None:
    print(f"[{time.strftime('%F %T')}] {message}", flush=True)


def atomic_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def copy_memmap(source: np.ndarray, target: np.ndarray, chunk: int) -> None:
    if source.shape != target.shape:
        raise RuntimeError(f"shape mismatch {source.shape} != {target.shape}")
    for start in range(0, source.shape[0], chunk):
        stop = min(start + chunk, source.shape[0])
        target[start:stop] = source[start:stop]


def load_manifest(pilot_dir: Path) -> dict:
    return json.loads((pilot_dir / "pilot_manifest.json").read_text(encoding="utf-8"))


def phase_merge(args: argparse.Namespace) -> None:
    output = args.output
    marker = output / "MERGE_SUCCESS"
    if marker.exists() and not args.force:
        log(f"merge already complete: {marker}")
        return
    output.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(args.pilot_dir)
    shards = sorted(manifest["shards"], key=lambda item: item["shard"])
    if len(shards) != 16:
        raise RuntimeError(f"expected 16 shards, observed {len(shards)}")
    for item in shards:
        shard = int(item["shard"])
        success = args.shards_dir / f"shard_{shard:02d}" / "SUCCESS"
        if not success.exists():
            raise RuntimeError(f"incomplete shard: {success}")

    total_records = int(sum(int(item["records"]) for item in shards))
    total_nt = int(sum(int(item["nt"]) for item in shards))
    if total_records != int(manifest["selected_records"]):
        raise RuntimeError("record total disagrees with pilot manifest")
    if total_nt != int(manifest["selected_nt"]):
        raise RuntimeError("nucleotide total disagrees with pilot manifest")

    first_dir = args.shards_dir / "shard_00"
    schemas = []
    for item in shards:
        shard_dir = args.shards_dir / f"shard_{int(item['shard']):02d}"
        schemas.append(json.loads((shard_dir / "schema.json").read_text(encoding="utf-8")))
    if any(schema != schemas[0] for schema in schemas[1:]):
        raise RuntimeError("shard schemas differ")
    atomic_json(output / "schema.json", schemas[0])

    specifications: Dict[str, Tuple[Tuple[int, ...], np.dtype]] = {}
    for name in RECORD_ARRAYS:
        source = np.load(first_dir / f"{name}.npy", mmap_mode="r")
        specifications[name] = ((total_records,) + source.shape[1:], source.dtype)
    for name in NT_ARRAYS:
        source = np.load(first_dir / f"{name}.npy", mmap_mode="r")
        specifications[name] = ((total_nt,) + source.shape[1:], source.dtype)
    specifications["offsets"] = ((total_records + 1,), np.dtype(np.uint64))

    arrays = {
        name: open_memmap(output / f"{name}.npy", mode="w+", dtype=dtype, shape=shape)
        for name, (shape, dtype) in specifications.items()
    }
    arrays["offsets"][0] = 0
    row_cursor = 0
    nt_cursor = 0
    shard_bounds = []
    for item in shards:
        shard = int(item["shard"])
        records = int(item["records"])
        nt = int(item["nt"])
        shard_dir = args.shards_dir / f"shard_{shard:02d}"
        log(f"merging shard {shard:02d}: records={records:,} nt={nt:,}")
        row_slice = slice(row_cursor, row_cursor + records)
        nt_slice = slice(nt_cursor, nt_cursor + nt)
        for name in RECORD_ARRAYS:
            source = np.load(shard_dir / f"{name}.npy", mmap_mode="r")
            copy_memmap(source, arrays[name][row_slice], 8192)
        for name in NT_ARRAYS:
            source = np.load(shard_dir / f"{name}.npy", mmap_mode="r")
            copy_memmap(source, arrays[name][nt_slice], 1_000_000)
        local_offsets = np.load(shard_dir / "offsets.npy", mmap_mode="r")
        arrays["offsets"][row_cursor + 1 : row_cursor + records + 1] = (
            local_offsets[1:] + np.uint64(nt_cursor)
        )
        shard_bounds.append(
            {
                "shard": shard,
                "record_start": row_cursor,
                "record_stop": row_cursor + records,
                "nt_start": nt_cursor,
                "nt_stop": nt_cursor + nt,
            }
        )
        row_cursor += records
        nt_cursor += nt
    for array in arrays.values():
        array.flush()
    if row_cursor != total_records or nt_cursor != total_nt:
        raise RuntimeError("merge cursor mismatch")

    ids = np.load(output / "ids.npy", mmap_mode="r")
    ids_path = output / "pilot_ids_with_index.tsv"
    with ids_path.open("wt", encoding="ascii", newline="\n", buffering=1 << 20) as handle:
        for start in range(0, total_records, 100_000):
            stop = min(start + 100_000, total_records)
            lines = [
                f"{index}\t{value.decode('ascii')}\n"
                for index, value in zip(range(start, stop), ids[start:stop])
            ]
            handle.writelines(lines)
    payload = {
        "schema_version": 1,
        "records": total_records,
        "nucleotides": total_nt,
        "source_manifest": manifest,
        "shard_bounds": shard_bounds,
    }
    atomic_json(output / "merged_manifest.json", payload)
    marker.write_text(json.dumps(payload) + "\n", encoding="utf-8")
    log(f"merge complete: {total_records:,} records; {total_nt:,} nt")


def create_array(path: Path, shape: Tuple[int, ...], dtype, fill=None) -> np.memmap:
    array = open_memmap(path, mode="w+", dtype=dtype, shape=shape)
    if fill is not None:
        array[...] = fill
    return array


def phase_annotations(args: argparse.Namespace) -> None:
    output = args.output
    marker = output / "ANNOTATION_SUCCESS"
    if marker.exists() and not args.force:
        log(f"annotations already complete: {marker}")
        return
    metadata_dir = args.metadata_dir
    if not (metadata_dir / "METADATA_SUCCESS").exists():
        raise RuntimeError(f"metadata filter is incomplete: {metadata_dir}")
    merged = json.loads((output / "merged_manifest.json").read_text(encoding="utf-8"))
    n = int(merged["records"])
    lengths = np.load(output / "lengths.npy", mmap_mode="r")

    permissive = create_array(output / "rfam_permissive.npy", (n,), "S8", b"")
    strict = create_array(output / "rfam_strict.npy", (n,), "S8", b"")
    best_evalue = create_array(output / "rfam_evalue.npy", (n,), np.float32, np.inf)
    best_bitscore = create_array(output / "rfam_bitscore.npy", (n,), np.float32, 0)
    best_coverage = create_array(output / "rfam_coverage.npy", (n,), np.float16, 0)
    family_count = create_array(output / "rfam_family_count.npy", (n,), np.uint16, 0)
    strict_rank: Dict[int, Tuple[float, float, float, bytes]] = {}
    descriptions: Dict[str, str] = {}

    current_index = -1
    observed_families: set = set()

    def finish_group() -> None:
        if current_index >= 0:
            family_count[current_index] = min(len(observed_families), 65535)

    kept_rows = 0
    with (metadata_dir / "rfam_filtered.tsv").open(
        "rt", encoding="utf-8", errors="replace", buffering=1 << 20
    ) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 9)
            if len(fields) < 10:
                continue
            index = int(fields[0])
            if index != current_index:
                finish_group()
                current_index = index
                observed_families = set()
            family = fields[2]
            observed_families.add(family)
            try:
                bitscore = float(fields[3])
                evalue = float(fields[4])
                seq_begin = int(fields[5])
                seq_end = int(fields[6])
            except ValueError:
                continue
            coverage = min(1.0, (abs(seq_end - seq_begin) + 1) / max(int(lengths[index]), 1))
            descriptions.setdefault(family, fields[9])
            rank = (evalue, -bitscore, -coverage)
            previous = (float(best_evalue[index]), -float(best_bitscore[index]), -float(best_coverage[index]))
            if rank < previous:
                permissive[index] = family.encode("ascii")
                best_evalue[index] = evalue
                best_bitscore[index] = bitscore
                best_coverage[index] = coverage
            is_strict = evalue <= 1e-10 and bitscore >= 20.0 and coverage >= 0.50
            if is_strict:
                previous_strict = strict_rank.get(index)
                candidate = (evalue, -bitscore, -coverage, family.encode("ascii"))
                if previous_strict is None or candidate[:3] < previous_strict[:3]:
                    strict_rank[index] = candidate
                    strict[index] = candidate[3]
            kept_rows += 1
            if kept_rows % 1_000_000 == 0:
                log(f"parsed {kept_rows:,} filtered Rfam rows")
    finish_group()
    strict_rank.clear()

    mapping_rows = create_array(output / "mapping_rows.npy", (n,), np.uint32, 0)
    mapping_taxids = create_array(output / "mapping_unique_taxids.npy", (n,), np.uint32, 0)
    taxon_mask = create_array(output / "taxon_mask.npy", (n,), np.uint8, 0)
    mapping_type = create_array(output / "mapping_primary_type.npy", (n,), "S48", b"")
    source_count = create_array(output / "mapping_source_count.npy", (n,), np.uint8, 0)
    mapping_ids = 0
    with (metadata_dir / "id_mapping_summary.tsv").open(
        "rt", encoding="utf-8", errors="replace", buffering=1 << 20
    ) as handle:
        next(handle, None)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            index = int(fields[0])
            mapping_rows[index] = min(int(fields[2]), np.iinfo(np.uint32).max)
            mapping_taxids[index] = min(int(fields[3]), np.iinfo(np.uint32).max)
            taxon_mask[index] = int(fields[4])
            mapping_type[index] = fields[5].encode("utf-8")[:48]
            source_count[index] = min(len([x for x in fields[6].split(";") if x]), 255)
            mapping_ids += 1

    for array in (
        permissive,
        strict,
        best_evalue,
        best_bitscore,
        best_coverage,
        family_count,
        mapping_rows,
        mapping_taxids,
        taxon_mask,
        mapping_type,
        source_count,
    ):
        array.flush()
    atomic_json(output / "rfam_descriptions.json", descriptions)
    payload = {
        "rfam_filtered_rows": kept_rows,
        "rfam_permissive_sequences": int(np.count_nonzero(permissive != b"")),
        "rfam_strict_sequences": int(np.count_nonzero(strict != b"")),
        "id_mapping_sequences": mapping_ids,
        "strict_definition": "E-value <= 1e-10, bitscore >= 20, sequence coverage >= 0.50",
        "permissive_definition": "Any match present in RNAcentral release 28 rfam_annotations.tsv.gz",
    }
    atomic_json(output / "annotation_manifest.json", payload)
    marker.write_text(json.dumps(payload) + "\n", encoding="utf-8")
    log(f"annotations complete: {payload}")


def normalized_rows(values: np.ndarray) -> np.ndarray:
    result = np.asarray(values, dtype=np.float32)
    norms = np.linalg.norm(result, axis=1, keepdims=True)
    result /= np.maximum(norms, 1e-12)
    return result


def phase_pca_knn(args: argparse.Namespace) -> None:
    output = args.output
    marker = output / "KNN_SUCCESS"
    if marker.exists() and not args.force:
        log(f"PCA/kNN already complete: {marker}")
        return
    from sklearn.decomposition import PCA
    import faiss

    embeddings = np.load(output / "embeddings.npy", mmap_mode="r")
    n, input_dim = embeddings.shape
    rng = np.random.default_rng(SEED)
    pca_sample_size = min(args.pca_sample, n)
    sample_indices = np.sort(rng.choice(n, size=pca_sample_size, replace=False))
    log(f"loading {pca_sample_size:,} embeddings for PCA")
    sample = normalized_rows(embeddings[sample_indices])
    pca = PCA(
        n_components=args.pca_dim,
        svd_solver="randomized",
        iterated_power=3,
        random_state=SEED,
    )
    log(f"fitting PCA {input_dim} -> {args.pca_dim}")
    pca.fit(sample)
    with (output / "pca_model.pkl").open("wb") as handle:
        pickle.dump(pca, handle, protocol=pickle.HIGHEST_PROTOCOL)
    del sample

    reduced = open_memmap(
        output / f"pca{args.pca_dim}.npy",
        mode="w+",
        dtype=np.float32,
        shape=(n, args.pca_dim),
    )
    for start in range(0, n, args.transform_batch):
        stop = min(start + args.transform_batch, n)
        batch = normalized_rows(embeddings[start:stop])
        batch = np.asarray(pca.transform(batch), dtype=np.float32)
        batch = normalized_rows(batch)
        reduced[start:stop] = batch
        if start % (args.transform_batch * 8) == 0:
            log(f"PCA transformed {stop:,}/{n:,}")
    reduced.flush()

    faiss.omp_set_num_threads(args.faiss_threads)
    d = args.pca_dim
    desired = max(1024, int(4 * math.sqrt(n)))
    nlist = 1 << int(round(math.log2(desired)))
    nlist = min(max(nlist, 1024), 16384)
    pq_m = args.pq_m
    if d % pq_m != 0:
        raise RuntimeError(f"PCA dimension {d} must be divisible by PQ m={pq_m}")
    quantizer = faiss.IndexFlatIP(d)
    index = faiss.IndexIVFPQ(
        quantizer, d, nlist, pq_m, 8, faiss.METRIC_INNER_PRODUCT
    )
    train_size = min(max(nlist * 48, args.knn_train), n)
    train_indices = np.sort(rng.choice(n, size=train_size, replace=False))
    train = np.ascontiguousarray(reduced[train_indices], dtype=np.float32)
    log(f"training FAISS IVFPQ: nlist={nlist} m={pq_m} train={train_size:,}")
    index.train(train)
    del train
    for start in range(0, n, args.faiss_batch):
        stop = min(start + args.faiss_batch, n)
        index.add(np.ascontiguousarray(reduced[start:stop], dtype=np.float32))
        if start % (args.faiss_batch * 8) == 0:
            log(f"FAISS added {stop:,}/{n:,}")
    index.nprobe = args.nprobe
    faiss.write_index(index, str(output / "pca_ivfpq.faiss"))

    search_k = args.neighbors + 1
    knn_indices = open_memmap(
        output / "knn_indices.npy", mode="w+", dtype=np.int32, shape=(n, search_k)
    )
    knn_similarity = open_memmap(
        output / "knn_similarity.npy", mode="w+", dtype=np.float16, shape=(n, search_k)
    )
    for start in range(0, n, args.faiss_batch):
        stop = min(start + args.faiss_batch, n)
        distance, neighbor = index.search(
            np.ascontiguousarray(reduced[start:stop], dtype=np.float32), search_k
        )
        knn_indices[start:stop] = neighbor.astype(np.int32)
        knn_similarity[start:stop] = distance.astype(np.float16)
        if start % (args.faiss_batch * 8) == 0:
            log(f"FAISS searched {stop:,}/{n:,}")
    knn_indices.flush()
    knn_similarity.flush()
    payload = {
        "records": n,
        "input_dimension": input_dim,
        "pca_dimension": args.pca_dim,
        "pca_sample": pca_sample_size,
        "pca_explained_variance": float(pca.explained_variance_ratio_.sum()),
        "faiss_index": "IVFPQ",
        "nlist": nlist,
        "nprobe": args.nprobe,
        "pq_m": pq_m,
        "neighbors": args.neighbors,
    }
    atomic_json(output / "knn_manifest.json", payload)
    marker.write_text(json.dumps(payload) + "\n", encoding="utf-8")
    log(f"PCA/kNN complete: {payload}")


def encode_fixed_strings(values: np.ndarray, dictionary: Mapping[str, int]) -> np.ndarray:
    unique, inverse = np.unique(values, return_inverse=True)
    mapped = np.full(len(unique), -1, dtype=np.int32)
    for index, value in enumerate(unique):
        text = value.decode("ascii") if isinstance(value, (bytes, np.bytes_)) else str(value)
        if text:
            mapped[index] = dictionary.get(text, -1)
    return mapped[inverse]


def wilson_lower(successes: int, total: int, z: float = 1.96) -> float:
    if total <= 0:
        return 0.0
    p = successes / total
    z2 = z * z
    return (
        p
        + z2 / (2 * total)
        - z * math.sqrt((p * (1 - p) + z2 / (4 * total)) / total)
    ) / (1 + z2 / total)


def phase_cluster(args: argparse.Namespace) -> None:
    output = args.output
    marker = output / "CLUSTER_SUCCESS"
    if marker.exists() and not args.force:
        log(f"clustering already complete: {marker}")
        return
    import igraph as ig
    import leidenalg

    knn = np.load(output / "knn_indices.npy", mmap_mode="r")
    similarity = np.load(output / "knn_similarity.npy", mmap_mode="r")
    n, search_k = knn.shape
    maximum_edges = n * search_k
    edges = np.empty((maximum_edges, 2), dtype=np.int32)
    weights = np.empty(maximum_edges, dtype=np.float32)
    cursor = 0
    chunk = 100_000
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        rows = np.repeat(np.arange(start, stop, dtype=np.int32), search_k)
        neighbors = np.asarray(knn[start:stop]).reshape(-1)
        scores = np.asarray(similarity[start:stop], dtype=np.float32).reshape(-1)
        valid = (neighbors >= 0) & (neighbors != rows) & np.isfinite(scores)
        count = int(valid.sum())
        edges[cursor : cursor + count, 0] = rows[valid]
        edges[cursor : cursor + count, 1] = neighbors[valid]
        weights[cursor : cursor + count] = np.maximum(scores[valid], 1e-3) ** 2
        cursor += count
    edges = edges[:cursor]
    weights = weights[:cursor]
    log(f"constructing igraph: vertices={n:,} directed-kNN-edges={cursor:,}")
    graph = ig.Graph(n=n, edges=edges, directed=False)
    del edges

    resolution = args.resolution
    log(f"running Leiden resolution={resolution}")
    partition = leidenalg.find_partition(
        graph,
        leidenalg.RBConfigurationVertexPartition,
        weights=weights,
        resolution_parameter=resolution,
        n_iterations=args.leiden_iterations,
        seed=SEED,
    )
    labels = np.asarray(partition.membership, dtype=np.int32)
    cluster_count = int(labels.max()) + 1
    if cluster_count < args.min_clusters:
        resolution *= 2.0
        log(f"only {cluster_count:,} clusters; retrying resolution={resolution}")
        partition = leidenalg.find_partition(
            graph,
            leidenalg.RBConfigurationVertexPartition,
            weights=weights,
            resolution_parameter=resolution,
            n_iterations=args.leiden_iterations,
            seed=SEED,
        )
        labels = np.asarray(partition.membership, dtype=np.int32)
        cluster_count = int(labels.max()) + 1
    elif cluster_count > args.max_clusters:
        resolution *= 0.5
        log(f"{cluster_count:,} clusters; retrying resolution={resolution}")
        partition = leidenalg.find_partition(
            graph,
            leidenalg.RBConfigurationVertexPartition,
            weights=weights,
            resolution_parameter=resolution,
            n_iterations=args.leiden_iterations,
            seed=SEED,
        )
        labels = np.asarray(partition.membership, dtype=np.int32)
        cluster_count = int(labels.max()) + 1
    np.save(output / "cluster_id.npy", labels)
    del graph, partition, weights

    schema = json.loads((output / "schema.json").read_text(encoding="utf-8"))
    descriptions = json.loads((output / "rfam_descriptions.json").read_text(encoding="utf-8"))
    all_families = sorted(set(descriptions) | set(schema["ncrna88"]))
    family_to_index = {family: index for index, family in enumerate(all_families)}
    atomic_json(
        output / "family_dictionary.json",
        {
            "families": [
                {"accession": family, "description": descriptions.get(family, "")}
                for family in all_families
            ]
        },
    )
    strict = np.load(output / "rfam_strict.npy", mmap_mode="r")
    permissive = np.load(output / "rfam_permissive.npy", mmap_mode="r")
    strict_code = encode_fixed_strings(strict, family_to_index)
    permissive_code = encode_fixed_strings(permissive, family_to_index)
    strict_indices = np.flatnonzero(strict_code >= 0)
    strict_total = np.bincount(labels[strict_indices], minlength=cluster_count).astype(np.int32)

    top_family = np.full(cluster_count, -1, dtype=np.int32)
    top_count = np.zeros(cluster_count, dtype=np.int32)
    if len(strict_indices):
        pairs = (
            labels[strict_indices].astype(np.int64) * len(all_families)
            + strict_code[strict_indices].astype(np.int64)
        )
        unique_pairs, counts = np.unique(pairs, return_counts=True)
        for pair, count in zip(unique_pairs, counts):
            cluster = int(pair // len(all_families))
            family = int(pair % len(all_families))
            if count > top_count[cluster]:
                top_count[cluster] = int(count)
                top_family[cluster] = family

    purity = np.divide(
        top_count,
        strict_total,
        out=np.zeros(cluster_count, dtype=np.float32),
        where=strict_total > 0,
    )
    lower = np.asarray(
        [wilson_lower(int(a), int(b)) for a, b in zip(top_count, strict_total)],
        dtype=np.float32,
    )
    consensus_family = np.where(
        (strict_total >= args.min_anchors)
        & (purity >= args.min_purity)
        & (lower >= args.min_wilson),
        top_family,
        -1,
    ).astype(np.int32)

    cluster_size = np.bincount(labels, minlength=cluster_count).astype(np.int32)
    first_order = np.argsort(labels, kind="stable")
    first_positions = np.r_[0, np.flatnonzero(np.diff(labels[first_order])) + 1]
    representative = first_order[first_positions].astype(np.int32)
    if len(representative) != cluster_count:
        raise RuntimeError("failed to select one representative per cluster")

    model88_id = np.load(output / "ncrna88_id.npy", mmap_mode="r")
    model88_conf = np.asarray(
        np.load(output / "ncrna88_confidence.npy", mmap_mode="r"), dtype=np.float32
    )
    model_family_lookup = np.asarray(
        [family_to_index[family] for family in schema["ncrna88"]], dtype=np.int32
    )
    model_family = model_family_lookup[np.asarray(model88_id, dtype=np.int64)]

    final_family = np.full(n, -1, dtype=np.int32)
    status = np.zeros(n, dtype=np.uint8)
    has_strict = strict_code >= 0
    has_permissive = permissive_code >= 0
    final_family[has_strict] = strict_code[has_strict]
    status[has_strict] = 1
    permissive_only = ~has_strict & has_permissive
    final_family[permissive_only] = permissive_code[permissive_only]
    status[permissive_only] = 2
    model_high = model88_conf >= args.model_confidence
    model_only = (final_family < 0) & model_high
    final_family[model_only] = model_family[model_only]
    status[model_only] = 3
    propagated = consensus_family[labels]
    propagated_only = (final_family < 0) & (propagated >= 0)
    final_family[propagated_only] = propagated[propagated_only]
    status[propagated_only] = 4
    conflict = has_permissive & model_high & (permissive_code != model_family)
    status[conflict] = 6

    anchor_fraction = strict_total / np.maximum(cluster_size, 1)
    novel_cluster = (
        (cluster_size >= args.novel_min_size)
        & (anchor_fraction < args.novel_max_anchor_fraction)
        & (consensus_family < 0)
    )
    novel_point = (
        (final_family < 0)
        & (model88_conf < args.novel_max_model_confidence)
        & novel_cluster[labels]
    )
    status[novel_point] = 5

    np.save(output / "cluster_size.npy", cluster_size)
    np.save(output / "cluster_representative.npy", representative)
    np.save(output / "cluster_strict_anchors.npy", strict_total)
    np.save(output / "cluster_consensus_family.npy", consensus_family)
    np.save(output / "cluster_purity.npy", purity)
    np.save(output / "cluster_wilson_lower.npy", lower)
    np.save(output / "cluster_novel_candidate.npy", novel_cluster)
    np.save(output / "final_family_id.npy", final_family)
    np.save(output / "annotation_status.npy", status)

    payload = {
        "records": n,
        "edges": cursor,
        "clusters": cluster_count,
        "resolution": resolution,
        "leiden_iterations": args.leiden_iterations,
        "strict_annotated": int(has_strict.sum()),
        "permissive_annotated": int(has_permissive.sum()),
        "cluster_consensus_count": int((consensus_family >= 0).sum()),
        "propagated_sequences": int(propagated_only.sum()),
        "novel_candidate_clusters": int(novel_cluster.sum()),
        "novel_candidate_sequences": int(novel_point.sum()),
        "conflicts": int(conflict.sum()),
        "annotation_status_codes": {
            "0": "unresolved",
            "1": "direct_rfam_strict",
            "2": "direct_rfam_permissive",
            "3": "ribogram_88_head_high_confidence",
            "4": "leiden_cluster_consensus_propagated",
            "5": "computational_novel_candidate",
            "6": "direct_rfam_vs_model_conflict",
        },
    }
    atomic_json(output / "cluster_manifest.json", payload)
    marker.write_text(json.dumps(payload) + "\n", encoding="utf-8")
    log(f"clustering complete: {payload}")


def phase_display(args: argparse.Namespace) -> None:
    output = args.output
    marker = output / "DISPLAY_SUCCESS"
    if marker.exists() and not args.force:
        log(f"display projection already complete: {marker}")
        return
    import faiss
    import umap

    labels = np.load(output / "cluster_id.npy", mmap_mode="r")
    representatives = np.load(output / "cluster_representative.npy")
    status = np.load(output / "annotation_status.npy", mmap_mode="r")
    n = len(labels)
    rng = np.random.default_rng(SEED)
    chosen = set(int(x) for x in representatives)

    novel = np.flatnonzero(status == 5)
    if len(novel):
        count = min(args.display_novel, len(novel))
        chosen.update(int(x) for x in rng.choice(novel, size=count, replace=False))
    direct = np.flatnonzero((status == 1) | (status == 2) | (status == 6))
    if len(direct):
        count = min(args.display_direct, len(direct))
        chosen.update(int(x) for x in rng.choice(direct, size=count, replace=False))
    remaining = args.display_points - len(chosen)
    if remaining > 0:
        candidates = np.setdiff1d(np.arange(n, dtype=np.int32), np.fromiter(chosen, dtype=np.int32))
        count = min(remaining, len(candidates))
        chosen.update(int(x) for x in rng.choice(candidates, size=count, replace=False))
    display_indices = np.asarray(sorted(chosen), dtype=np.int32)
    if len(display_indices) > args.display_points:
        mandatory = np.asarray(sorted(set(int(x) for x in representatives)), dtype=np.int32)
        keep = max(0, args.display_points - len(mandatory))
        optional = np.setdiff1d(display_indices, mandatory)
        optional = rng.choice(optional, size=keep, replace=False) if keep else np.empty(0, np.int32)
        display_indices = np.sort(np.r_[mandatory, optional]).astype(np.int32)
    np.save(output / "display_indices.npy", display_indices)
    log(f"display subset selected: {len(display_indices):,}/{n:,}")

    knn_manifest = json.loads((output / "knn_manifest.json").read_text(encoding="utf-8"))
    dimension = int(knn_manifest["pca_dimension"])
    reduced = np.load(output / f"pca{dimension}.npy", mmap_mode="r")
    values = np.ascontiguousarray(reduced[display_indices], dtype=np.float32)
    count = len(values)
    faiss.omp_set_num_threads(args.faiss_threads)
    nlist = min(4096, max(256, 1 << int(round(math.log2(math.sqrt(count) * 3)))))
    quantizer = faiss.IndexFlatIP(dimension)
    index = faiss.IndexIVFFlat(quantizer, dimension, nlist, faiss.METRIC_INNER_PRODUCT)
    train_count = min(150_000, count)
    train_index = rng.choice(count, size=train_count, replace=False)
    index.train(np.ascontiguousarray(values[train_index]))
    index.add(values)
    index.nprobe = min(64, nlist)
    distances, neighbors = index.search(values, args.umap_neighbors + 1)
    umap_indices = np.empty((count, args.umap_neighbors), dtype=np.int32)
    umap_distances = np.empty((count, args.umap_neighbors), dtype=np.float32)
    for row in range(count):
        candidates = neighbors[row]
        candidate_scores = distances[row]
        mask = candidates != row
        take_i = candidates[mask][: args.umap_neighbors - 1]
        take_s = candidate_scores[mask][: args.umap_neighbors - 1]
        if len(take_i) < args.umap_neighbors - 1:
            pad = args.umap_neighbors - 1 - len(take_i)
            take_i = np.r_[take_i, np.full(pad, row, dtype=np.int64)]
            take_s = np.r_[take_s, np.ones(pad, dtype=np.float32)]
        umap_indices[row, 0] = row
        umap_indices[row, 1:] = take_i
        umap_distances[row, 0] = 0.0
        umap_distances[row, 1:] = np.sqrt(np.maximum(0.0, 2.0 - 2.0 * take_s))
    del index, values, distances, neighbors

    reducer = umap.UMAP(
        n_neighbors=args.umap_neighbors,
        n_components=2,
        min_dist=args.umap_min_dist,
        spread=1.0,
        metric="euclidean",
        init="random",
        n_epochs=args.umap_epochs,
        n_jobs=args.umap_jobs,
        random_state=None,
        low_memory=False,
        verbose=True,
        precomputed_knn=(umap_indices, umap_distances, None),
    )
    log("optimizing display UMAP")
    coordinates = reducer.fit_transform(np.zeros((count, 1), dtype=np.float32))
    coordinates = np.asarray(coordinates, dtype=np.float32)
    center = np.median(coordinates, axis=0)
    scale = np.quantile(np.abs(coordinates - center), 0.995, axis=0)
    coordinates = (coordinates - center) / np.maximum(scale, 1e-6)
    np.save(output / "display_xy.npy", coordinates)
    with (output / "display_umap.pkl").open("wb") as handle:
        pickle.dump(reducer, handle, protocol=pickle.HIGHEST_PROTOCOL)
    payload = {
        "points": count,
        "mandatory_cluster_representatives": int(len(representatives)),
        "umap_neighbors": args.umap_neighbors,
        "umap_min_dist": args.umap_min_dist,
        "umap_epochs": args.umap_epochs,
        "projection_note": "UMAP is a stochastic display projection of a stratified pilot subset; distances are not a phylogeny.",
    }
    atomic_json(output / "display_manifest.json", payload)
    marker.write_text(json.dumps(payload) + "\n", encoding="utf-8")
    log(f"display projection complete: {payload}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    subparsers = parser.add_subparsers(dest="phase", required=True)

    merge = subparsers.add_parser("merge")
    merge.add_argument("--pilot-dir", type=Path, required=True)
    merge.add_argument("--shards-dir", type=Path, required=True)
    merge.set_defaults(function=phase_merge)

    annotations = subparsers.add_parser("annotations")
    annotations.add_argument("--metadata-dir", type=Path, required=True)
    annotations.set_defaults(function=phase_annotations)

    knn = subparsers.add_parser("pca-knn")
    knn.add_argument("--pca-dim", type=int, default=64)
    knn.add_argument("--pca-sample", type=int, default=120000)
    knn.add_argument("--transform-batch", type=int, default=32768)
    knn.add_argument("--faiss-batch", type=int, default=32768)
    knn.add_argument("--faiss-threads", type=int, default=128)
    knn.add_argument("--knn-train", type=int, default=400000)
    knn.add_argument("--pq-m", type=int, default=16)
    knn.add_argument("--nprobe", type=int, default=32)
    knn.add_argument("--neighbors", type=int, default=16)
    knn.set_defaults(function=phase_pca_knn)

    cluster = subparsers.add_parser("cluster")
    cluster.add_argument("--resolution", type=float, default=1.0)
    cluster.add_argument("--leiden-iterations", type=int, default=2)
    cluster.add_argument("--min-clusters", type=int, default=2000)
    cluster.add_argument("--max-clusters", type=int, default=80000)
    cluster.add_argument("--min-anchors", type=int, default=5)
    cluster.add_argument("--min-purity", type=float, default=0.80)
    cluster.add_argument("--min-wilson", type=float, default=0.60)
    cluster.add_argument("--model-confidence", type=float, default=0.85)
    cluster.add_argument("--novel-min-size", type=int, default=20)
    cluster.add_argument("--novel-max-anchor-fraction", type=float, default=0.05)
    cluster.add_argument("--novel-max-model-confidence", type=float, default=0.50)
    cluster.set_defaults(function=phase_cluster)

    display = subparsers.add_parser("display")
    display.add_argument("--display-points", type=int, default=300000)
    display.add_argument("--display-novel", type=int, default=30000)
    display.add_argument("--display-direct", type=int, default=60000)
    display.add_argument("--faiss-threads", type=int, default=128)
    display.add_argument("--umap-neighbors", type=int, default=30)
    display.add_argument("--umap-min-dist", type=float, default=0.08)
    display.add_argument("--umap-epochs", type=int, default=180)
    display.add_argument("--umap-jobs", type=int, default=128)
    display.set_defaults(function=phase_display)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    args.function(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
