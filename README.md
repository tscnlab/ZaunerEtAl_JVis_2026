<img src="https://github.com/tscnlab/Templates/blob/main/logo/logo_with_text-01.png" width="400" />

# ZaunerEtAl_JVis_2026

This repository contains the Quarto project for the manuscript **"Analysis of human visual experience data"**, including the main article, supplementary material, and supporting assets. The website is automatically built and published to GitHub Pages.

- 📄 Main manuscript: `index.qmd`
- 📑 Supplementary material: `supplement.qmd`
- 🗂️ Data resources: `data/` (see `_quarto.yml` for included files)
- 🎨 Styling and extensions: `styles.css`, `_extensions/`
- 📜 Citation information: `CITATION.cff`

## View the site
The live, GitHub Pages–hosted version is available at: https://tscnlab.github.io/ZaunerEtAl_JVis_2026/

## Local development
1. Install Quarto and the project R dependencies (managed via `renv`).
2. Render the site locally:
   ```
   quarto preview
   ```
3. Published content will appear in the `output/` directory when rendered for deployment.
