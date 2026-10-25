#!/usr/bin/env python3
"""Export a compact, browser-oriented RNAcentral atlas bundle."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import time
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Sequence, Tuple

import numpy as np


def log(message: str) -> None:
    print(f"[{time.strftime('%F %T')}] {message}", flush=True)


def json_gzip(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0) as zipped:
            text = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
            zipped.write(text.encode("utf-8"))
    os.replace(temporary, path)


def fasta_records(paths: Sequence[Path]) -> Iterator[Tuple[str, str]]:
    identifier = ""
    sequence: List[str] = []
    for path in paths:
        with path.open("rt", encoding="ascii", errors="replace") as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if identifier:
                        yield identifier, "".join(sequence).upper().replace("T", "U")
                    identifier = line[1:].split()[0]
                    sequence = []
                else:
                    sequence.append(line)
            if identifier:
                yield identifier, "".join(sequence).upper().replace("T", "U")
                identifier = ""
                sequence = []


def decode_ids(values: np.ndarray) -> List[str]:
    return [value.decode("ascii") for value in values]


def rounded(values: np.ndarray, digits: int = 5) -> List[float]:
    return np.round(np.asarray(values, dtype=np.float64), digits).tolist()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--pilot-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--detail-records", type=int, default=12000)
    parser.add_argument("--detail-chunk", type=int, default=128)
    args = parser.parse_args()

    results = args.results
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    details_dir = output / "details"
    details_dir.mkdir(parents=True, exist_ok=True)

    display_indices = np.load(results / "display_indices.npy")
    xy = np.load(results / "display_xy.npy")
    ids = np.load(results / "ids.npy", mmap_mode="r")
    lengths = np.load(results / "lengths.npy", mmap_mode="r")
    cluster_id = np.load(results / "cluster_id.npy", mmap_mode="r")
    cluster_size = np.load(results / "cluster_size.npy")
    cluster_rep = np.load(results / "cluster_representative.npy")
    cluster_anchors = np.load(results / "cluster_strict_anchors.npy")
    cluster_consensus = np.load(results / "cluster_consensus_family.npy")
    cluster_purity = np.load(results / "cluster_purity.npy")
    cluster_lower = np.load(results / "cluster_wilson_lower.npy")
    cluster_novel = np.load(results / "cluster_novel_candidate.npy")
    final_family = np.load(results / "final_family_id.npy", mmap_mode="r")
    status = np.load(results / "annotation_status.npy", mmap_mode="r")
    type13 = np.load(results / "ncrna13_id.npy", mmap_mode="r")
    type13_conf = np.load(results / "ncrna13_confidence.npy", mmap_mode="r")
    model88_conf = np.load(results / "ncrna88_confidence.npy", mmap_mode="r")
    taxon_mask = np.load(results / "taxon_mask.npy", mmap_mode="r")
    taxa_count = np.load(results / "taxa_counts.npy", mmap_mode="r")
    source_count = np.load(results / "mapping_source_count.npy", mmap_mode="r")
    rfam_bitscore = np.load(results / "rfam_bitscore.npy", mmap_mode="r")
    rfam_evalue = np.load(results / "rfam_evalue.npy", mmap_mode="r")

    schema = json.loads((results / "schema.json").read_text(encoding="utf-8"))
    families = json.loads((results / "family_dictionary.json").read_text(encoding="utf-8"))["families"]
    cluster_manifest = json.loads((results / "cluster_manifest.json").read_text(encoding="utf-8"))
    annotation_manifest = json.loads((results / "annotation_manifest.json").read_text(encoding="utf-8"))
    display_manifest = json.loads((results / "display_manifest.json").read_text(encoding="utf-8"))
    merged_manifest = json.loads((results / "merged_manifest.json").read_text(encoding="utf-8"))

    cluster_count = len(cluster_size)
    displayed_clusters = np.asarray(cluster_id[display_indices], dtype=np.int32)
    centroid_x = np.bincount(displayed_clusters, weights=xy[:, 0], minlength=cluster_count)
    centroid_y = np.bincount(displayed_clusters, weights=xy[:, 1], minlength=cluster_count)
    displayed_per_cluster = np.bincount(displayed_clusters, minlength=cluster_count)
    centroid_x /= np.maximum(displayed_per_cluster, 1)
    centroid_y /= np.maximum(displayed_per_cluster, 1)

    cluster_rows = []
    for cluster in range(cluster_count):
        cluster_rows.append(
            [
                int(cluster_size[cluster]),
                int(cluster_rep[cluster]),
                int(cluster_anchors[cluster]),
                int(cluster_consensus[cluster]),
                round(float(cluster_purity[cluster]), 4),
                round(float(cluster_lower[cluster]), 4),
                1 if bool(cluster_novel[cluster]) else 0,
                round(float(centroid_x[cluster]), 5),
                round(float(centroid_y[cluster]), 5),
            ]
        )

    novel_clusters = np.flatnonzero(cluster_novel)
    consensus_clusters = np.flatnonzero(cluster_consensus >= 0)
    ranked_clusters: List[int] = []
    for candidates in (novel_clusters, consensus_clusters, np.arange(cluster_count)):
        order = candidates[np.argsort(-cluster_size[candidates], kind="stable")]
        ranked_clusters.extend(int(value) for value in order)
    seen_clusters = set()
    detail_global: List[int] = []
    for cluster in ranked_clusters:
        if cluster in seen_clusters:
            continue
        seen_clusters.add(cluster)
        detail_global.append(int(cluster_rep[cluster]))
        if len(detail_global) >= args.detail_records:
            break
    detail_set = set(detail_global)
    detail_chunk_by_global = {
        global_index: position // args.detail_chunk
        for position, global_index in enumerate(detail_global)
    }

    display_ids = decode_ids(ids[display_indices])
    point_columns = {
        "id": display_ids,
        "x": rounded(xy[:, 0]),
        "y": rounded(xy[:, 1]),
        "length": np.asarray(lengths[display_indices], dtype=np.int32).tolist(),
        "cluster": displayed_clusters.tolist(),
        "family": np.asarray(final_family[display_indices], dtype=np.int32).tolist(),
        "status": np.asarray(status[display_indices], dtype=np.int16).tolist(),
        "type13": np.asarray(type13[display_indices], dtype=np.int16).tolist(),
        "type13_conf": np.rint(np.asarray(type13_conf[display_indices], dtype=np.float32) * 100).astype(np.int16).tolist(),
        "model88_conf": np.rint(np.asarray(model88_conf[display_indices], dtype=np.float32) * 100).astype(np.int16).tolist(),
        "taxon_mask": np.asarray(taxon_mask[display_indices], dtype=np.int16).tolist(),
        "taxa_count": np.asarray(taxa_count[display_indices], dtype=np.int64).tolist(),
        "source_count": np.asarray(source_count[display_indices], dtype=np.int16).tolist(),
        "detail_chunk": [detail_chunk_by_global.get(int(index), -1) for index in display_indices],
    }

    status_counts = np.bincount(np.asarray(status, dtype=np.int64), minlength=7)
    payload = {
        "schema_version": 1,
        "title": "RNAcentral × RiboGRAM Universe",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "summary": {
            "rnacentral_release": 28,
            "pilot_sequences": int(merged_manifest["records"]),
            "pilot_nucleotides": int(merged_manifest["nucleotides"]),
            "displayed_sequences": int(len(display_indices)),
            "leiden_clusters": int(cluster_count),
            "direct_rfam_permissive": int(annotation_manifest["rfam_permissive_sequences"]),
            "direct_rfam_strict": int(annotation_manifest["rfam_strict_sequences"]),
            "propagated_sequences": int(cluster_manifest["propagated_sequences"]),
            "novel_candidate_clusters": int(cluster_manifest["novel_candidate_clusters"]),
        },
        "method": {
            "embedding": "RiboGRAM Layer-12 mean representation (1280 dimensions)",
            "reduction": "L2 normalization, randomized PCA-64, approximate cosine kNN",
            "clustering": "Leiden on a 16-neighbor FAISS IVFPQ graph",
            "projection": "UMAP of a stratified 300k display subset",
            "attention": "Layer-7 and Layer-12 head-mean APC received attention",
            "pairing": "Layer-12 Head-20 symmetrized APC; manuscript threshold 0.4",
            "warning": "Cluster-propagated and novel-candidate labels are computational hypotheses, not curated RNAcentral annotations.",
        },
        "dictionaries": {
            "families": families,
            "ncrna13": schema["ncrna13"],
            "annotation_status": cluster_manifest["annotation_status_codes"],
            "taxon_mask_bits": {"1": "Bacteria", "2": "Archaea", "4": "Eukaryota", "8": "Viruses"},
        },
        "cluster_columns": [
            "size",
            "representative_global_index",
            "strict_rfam_anchors",
            "consensus_family",
            "purity",
            "wilson_lower_95",
            "novel_candidate",
            "x",
            "y",
        ],
        "clusters": cluster_rows,
        "points": point_columns,
        "details": {
            "records": len(detail_global),
            "chunk_size": args.detail_chunk,
            "url_pattern": "details/chunk_{chunk:04d}.json.gz",
        },
        "status_counts_all_pilot": status_counts.tolist(),
        "projection": display_manifest,
    }
    json_gzip(output / "universe.json.gz", payload)
    log(f"wrote universe.json.gz with {len(display_indices):,} points")

    target_ids = {ids[index].decode("ascii"): index for index in detail_global}
    sequences: Dict[int, str] = {}
    fasta_paths = sorted(args.pilot_dir.glob("pilot_shard_*.fa"))
    for identifier, sequence in fasta_records(fasta_paths):
        index = target_ids.get(identifier)
        if index is not None:
            sequences[index] = sequence
    if len(sequences) != len(detail_global):
        missing = len(detail_global) - len(sequences)
        raise RuntimeError(f"failed to recover {missing} representative sequences")

    offsets = np.load(results / "offsets.npy", mmap_mode="r")
    tracks = np.load(results / "attention_tracks.npy", mmap_mode="r")
    partners = np.load(results / "pair_partner.npy", mmap_mode="r")
    total_chunks = math.ceil(len(detail_global) / args.detail_chunk)
    for chunk in range(total_chunks):
        begin = chunk * args.detail_chunk
        stop = min(begin + args.detail_chunk, len(detail_global))
        records = []
        for global_index in detail_global[begin:stop]:
            nt_begin = int(offsets[global_index])
            nt_stop = int(offsets[global_index + 1])
            score = np.asarray(tracks[nt_begin:nt_stop], dtype=np.float32)
            partner = np.asarray(partners[nt_begin:nt_stop], dtype=np.int32)
            partner[(partner == 65535) | (score[:, 2] < 0.4)] = -1
            family_index = int(final_family[global_index])
            records.append(
                {
                    "id": ids[global_index].decode("ascii"),
                    "sequence": sequences[global_index],
                    "global_index": global_index,
                    "cluster": int(cluster_id[global_index]),
                    "family": family_index,
                    "status": int(status[global_index]),
                    "rfam_bitscore": round(float(rfam_bitscore[global_index]), 3),
                    "rfam_evalue": None if not np.isfinite(rfam_evalue[global_index]) else float(rfam_evalue[global_index]),
                    "attention7": np.round(score[:, 0], 3).tolist(),
                    "attention12": np.round(score[:, 1], 3).tolist(),
                    "pair_score": np.round(score[:, 2], 3).tolist(),
                    "pair_partner": partner.tolist(),
                }
            )
        json_gzip(details_dir / f"chunk_{chunk:04d}.json.gz", {"records": records})
        if chunk % 10 == 0:
            log(f"wrote detail chunk {chunk + 1}/{total_chunks}")

    manifest = {
        "schema_version": 1,
        "main": "universe.json.gz",
        "main_bytes": (output / "universe.json.gz").stat().st_size,
        "detail_chunks": total_chunks,
        "detail_records": len(detail_global),
        "detail_bytes": sum(path.stat().st_size for path in details_dir.glob("*.json.gz")),
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "WEB_SUCCESS").write_text(json.dumps(manifest) + "\n", encoding="utf-8")
    log(f"web export complete: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
