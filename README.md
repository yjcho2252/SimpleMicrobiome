# SimpleMicrobiome

[![R](https://img.shields.io/badge/R-%3E%3D4.2-276DC3)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-App-1F77B4)](https://shiny.posit.co/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

SimpleMicrobiome is an R Shiny application for end-to-end microbiome analysis with a GUI.  
It integrates data loading, preprocessing, taxa profiling, diversity, differential abundance, network analysis, and association analysis in a single workflow.
Public app: https://simplemicrobiome.mglab.org/

## Current App Modules

- **Top**
  - App overview, quick workflow, version history, and external QZA converter link
- **Data Loading**
  - Upload ASV/OTU table, taxonomy table, metadata
  - Load/download bundled example data
- **Preprocessing**
  - Sample/taxa filtering before downstream analyses
- **Taxa Profiles**
  - Taxa Barplot
  - Taxa Comparison (including optional longitudinal/paired overlay and paired test mode when enabled)
- **Diversity**
  - Alpha Diversity
  - Beta Diversity
- **Differential Abundance**
  - ANCOM-BC2
  - MaAsLin2
  - Random Forest
- **Network**
  - SparCC
  - SpiecEasi
- **Association**
  - Correlation Heatmap
  - Association Biplot (dbRDA/CAP-based)
- **Citation**
  - Method/package references

## Requirements

- R (recommended: 4.2+)
- Main R packages:
  - shiny
  - phyloseq
  - tidyverse
  - ANCOMBC
  - DT
  - bslib
  - ggplot2
  - ggpubr
  - shinyWidgets
  - Maaslin2
  - microbiome
  - cluster
  - vegan
  - dplyr
  - randomForest
  - ComplexHeatmap
  - NetCoMi
  - igraph
  - ggraph
  - ggrepel
  - shapr

## Input Data

Upload three files in **Data Loading**:

1. ASV/OTU abundance matrix (`.csv`, `.tsv`, `.txt`)
2. Taxonomy table (`.csv`, `.tsv`, `.txt`)
3. Metadata file (`.csv`, `.tsv`, `.txt`)

Expected format:

- Abundance table: rows = taxa/features, columns = sample IDs, first column = taxa ID
- Taxonomy table: rows aligned to taxa IDs, first column = taxa ID
- Metadata: first column = sample IDs (mapped to `SampleID`), remaining columns = metadata variables

The app automatically:

- Detects delimiter (comma/tab)
- Intersects shared sample IDs between abundance and metadata
- Intersects shared taxa IDs between abundance and taxonomy
- Removes zero-sum taxa/samples

## Example Data

Bundled files:

- `sample/1_ASV_table.txt`
- `sample/2_taxonomy_table.txt`
- `sample/3_metadata.txt`

QZA conversion link (table.qza/taxonomy.qza):

- `https://simplemicrobiome.mglab.org/convert/ui`

## Recommended Workflow

1. Load data in **Data Loading**.
2. Apply filters in **Preprocessing**.
3. Run analysis modules in each menu (Taxa Profiles, Diversity, Differential Abundance, Network, Association).
4. Check **Citation** for method references.

## Troubleshooting

- Verify delimiter consistency and header rows.
- Verify sample ID matching between abundance and metadata.
- Verify taxa ID matching between abundance and taxonomy.
- If a module output is empty, revisit preprocessing filters and group sizes.

## License

This project is licensed under the MIT License.  
See [LICENSE](LICENSE) for details.
