# RiboGRAM Atlas

Two linked, buildless atlases make the representation space learned by RiboGRAM explorable in a browser:

- **RNAcentral Universe** — a deterministic release-28 pilot of **1,995,878** RNAs, with a stratified 300,000-point WebGL map, **296** evidence-calibrated Leiden neighborhoods and provenance-aware annotations
- **Rfam Seeds** — **92,781** Seed sequences spanning **4,178** Rfam 15.0 families
- Search, evidence filtering and coloring by RNA class, annotation status, taxonomy, family or cluster
- Per-nucleotide Layer-7/12 attention and zero-shot pairing tracks for representative RNAs
- Pan/zoom, point inspection and shareable RNAcentral accession URLs
- Static GitHub Pages deployment with no browser-side model inference

## Explore the atlas

**Launch the live website: [https://ribogram.github.io/](https://ribogram.github.io/)**

No installation is required. Open the link in a modern browser to explore the interactive landscapes.

## RiboGRAM resources

- Source code and analysis scripts: [Zenodo record 23201637](https://zenodo.org/records/23201637)
- Executable reproduction capsule: [Code Ocean capsule 2602506](https://codeocean.com/capsule/2602506/tree/v1)
- RiboGRAM resource portal: [ultrarnalab.com/ribogram.php](https://www.ultrarnalab.com/ribogram.php)
- RiboGRAM web server: [ultrarnalab.com/ribogram_usage.php](https://www.ultrarnalab.com/ribogram_usage.php)

## Deploy

The included GitHub Actions workflow publishes the repository root to GitHub Pages whenever `main` is updated. For an organization site, use the repository name `RiboGRAM.github.io` and set **Settings → Pages → Source** to **GitHub Actions**.

## Data provenance

The RNAcentral view is explicitly a sampled release-28 pilot rather than the complete release. Direct RNAcentral/Rfam evidence, RiboGRAM predictions, cluster-propagated labels and computational novel candidates remain separate in the interface. See [RNACENTRAL_ATLAS.md](RNACENTRAL_ATLAS.md) for sampling, feature, clustering and annotation definitions.

The Rfam view contains public Rfam 15.0 Seed sequences and precomputed UMAP coordinates from RiboGRAM Layer-12 representations. Detail bundles expose position-resolved contextual-dependency scores and symmetrized, APC-corrected Layer-12 Head-20 pairing signals. The interface does not perform model inference in the browser.

## Citation

Wang H., Chen Z., Huang H. et al. *RiboGRAM learns multi-scale non-coding RNA grammar for structural and regulatory inference.* bioRxiv (2026). DOI: [10.64898/2026.10.12.758266](https://doi.org/10.64898/2026.10.12.758266)
