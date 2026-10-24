# RiboGRAM Atlas

An immersive, single-page visualization of the representation space learned by RiboGRAM across the complete Rfam 15.0 Seed collection.

- **92,781** Seed sequences
- **4,178** Rfam families
- Search, family filtering and RNA-type/family coloring
- Semantic labels, pan/zoom, point inspection and shareable sequence URLs
- Buildless static deployment for GitHub Pages

## Explore the atlas

**Launch the live website: [https://ribogram.github.io/](https://ribogram.github.io/)**

No installation or local server is required. Open the link in a modern browser to explore the complete interactive RNA landscape.

## Deploy

The included GitHub Actions workflow publishes the repository root to GitHub Pages whenever `main` is updated. For an organization site, use the repository name `RiboGRAM.github.io` and set **Settings → Pages → Source** to **GitHub Actions**.

## Data provenance

The bundled projection contains public Rfam 15.0 Seed sequences and precomputed two-dimensional UMAP coordinates from RiboGRAM Layer-12 representations. The interface does not perform model inference in the browser.

## Citation

Wang H., Chen Z., Huang H. et al. *RiboGRAM learns multi-scale non-coding RNA grammar for structural and regulatory inference.* bioRxiv (2026). DOI: [10.64898/2026.10.12.758266](https://doi.org/10.64898/2026.10.12.758266)
