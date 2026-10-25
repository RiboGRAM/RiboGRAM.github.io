# RNAcentral × RiboGRAM atlas

This repository publishes a browser-oriented view of a reproducible pilot of the RNAcentral sequence landscape. It deliberately separates **database evidence** from **model predictions** and **cluster-level hypotheses**.

## Dataset

- RNAcentral release 28 active FASTA: 48,053,536 records / 41,338,415,722 nt
- Release-28 MD5 table: 61,026,481 rows
- Deterministic pilot: 1,995,878 sequences / 653,939,172 nt
- Eligibility: 20–1,024 nt and no more than 5% ambiguous symbols
- Sampling: deterministic hash sampling with seed `20261025`

The active FASTA and release-wide MD5 totals describe different release artifacts and must not be treated as interchangeable denominators. The public atlas reports the exact FASTA-derived pilot denominator.

## RiboGRAM features

For every pilot RNA, the pipeline stores:

1. Layer-12 mean sequence representation (1,280 dimensions, FP16)
2. Released 13-class and 88-family classifier predictions with confidence
3. Layer-7 head-mean APC received attention at every nucleotide
4. Layer-12 head-mean APC received attention at every nucleotide
5. Layer-12 Head-20 symmetrized/APC zero-shot top-pair score and partner at every nucleotide

Dense `L × L` attention matrices are never retained. The manuscript pairing threshold is `0.4`.

## Graph and projection

- L2 normalization
- Randomized PCA to 64 dimensions (trained on a deterministic 120,000-sequence sample)
- FAISS IVFPQ approximate cosine graph (`k=16`)
- Leiden community detection
- A stratified display subset containing every cluster representative plus direct-evidence, novel-candidate and random points
- UMAP is used only as a two-dimensional display projection; it is not a phylogeny

The selected Leiden resolution produced 296 neighborhoods. Against 1,117,305 strict Rfam anchors, weighted family purity was 0.920, macro purity across anchored clusters was 0.873 and median anchored-cluster purity was 1.000. Exact machine-readable manifests are published in [`data/universe/provenance/`](data/universe/provenance/).

## Annotation hierarchy

The website uses the following precedence without rewriting provenance:

1. **Direct Rfam, strict** — RNAcentral release-28 hit with E-value ≤ 1e-10, bitscore ≥ 20 and sequence coverage ≥ 0.50
2. **Direct Rfam, permissive** — any hit distributed in the release-28 Rfam annotation table
3. **RiboGRAM-88** — released classifier prediction with confidence ≥ 0.85
4. **Cluster consensus** — propagated only from sufficiently supported, high-purity strict Rfam anchors
5. **Unresolved**

Direct/model conflicts remain marked as conflicts. A **novel candidate** is a computationally unannotated neighborhood satisfying predefined graph/evidence rules. It is not a claim of a newly discovered biological family.

## Reproduction

The principal scripts are in `research/remote/`:

- `sample_rnacentral_pilot.cpp`
- `extract_rnacentral_atlas.py`
- `ribogram_attention_layers.patch` (limits retained attention to requested layers)
- `filter_rnacentral_metadata.cpp`
- `build_rnacentral_atlas.py`
- `export_rnacentral_web.py`

Run manifests store exact counts, thresholds, model schemas and graph parameters alongside the generated arrays.
