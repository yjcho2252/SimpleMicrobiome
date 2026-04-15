# SimpleMicrobiome

SimpleMicrobiome is an R Shiny application for end-to-end microbiome analysis with a GUI.
It supports data loading, preprocessing, visualization, diversity analysis, differential abundance testing, and network analysis in one workflow.

## Features

- Guided workflow from data upload to downstream analysis
- Interactive sample filtering and selection in preprocessing
- Taxa visualization:
  - Taxa Barplot
  - Taxa Comparison
- Diversity analysis:
  - Alpha Diversity
  - Beta Diversity
- Differential abundance and predictive modeling:
  - ANCOM-BC2
  - MaAsLin2
  - Random Forest
- Network analysis:
  - SparCC
  - SpiecEasi
- External QIIME2 converter service for `.qza` to TSV conversion (`/convert/ui`)
- Built-in example dataset loader/downloader
- Citation page with module-wise references

## Tech Stack

- Language: R
- App framework: Shiny + bslib
- Core microbiome object model: phyloseq
- Major analysis libraries: ANCOMBC, Maaslin2, vegan, randomForest, NetCoMi

## Requirements

- R (recommended: 4.2+)
- R packages:
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

## Installation

Clone the repository:

```bash
git clone https://github.com/<your-org-or-user>/SimpleMicrobiome.git
cd SimpleMicrobiome
```

Install required packages in R:

```r
install.packages(c(
  "shiny", "tidyverse", "DT", "bslib", "ggplot2", "ggpubr", "shinyWidgets",
  "microbiome", "cluster", "vegan", "dplyr", "randomForest", "igraph",
  "ggraph", "ggrepel", "shapr"
))

# Bioconductor / non-CRAN packages as needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "ANCOMBC", "ComplexHeatmap"))

# Maaslin2 and NetCoMi may require project-specific installation steps
```

## Run the App

From the project root:

```r
shiny::runApp()
```

Or:

```r
source("app.R")
```

## Input Data Format

Load three files in the **Data Loading** tab:

1. `ASV/OTU Abundance Matrix` (`.csv`, `.tsv`, `.txt`)
2. `Taxonomy Table` (`.csv`, `.tsv`, `.txt`)
3. `Metadata File` (`.csv`, `.tsv`, `.txt`)

Expected structure:

- ASV/OTU table:
  - Rows: taxa/features
  - Columns: sample IDs
  - First column: taxa ID
- Taxonomy table:
  - Rows aligned to taxa IDs in ASV/OTU table
  - First column: taxa ID
- Metadata:
  - First column must contain sample IDs (internally mapped to `SampleID`)
  - Remaining columns are metadata variables

The app automatically:

- Detects delimiter (tab or comma)
- Intersects shared sample IDs between abundance table and metadata
- Intersects shared taxa IDs between abundance table and taxonomy
- Removes zero-sum taxa and samples

## Example Data

The repository includes example input files in `sample/`:

- `sample/1_ASV_table.txt`
- `sample/2_taxonomy_table.txt`
- `sample/3_metadata.txt`

You can use **Load Example** in the app or **Download Example** to export a zip.

For QIIME2 `.qza` conversion, use the external converter page linked from **Top**:

- `https://simplemicrobiome.mglab.org/convert/ui`

## Recommended Workflow

1. Open **Data Loading** and upload the 3 files (or load the example).
2. Move to **Preprocessing** and filter/select samples.
3. Run one or more analysis modules:
   - Visualization
   - Diversity Analysis
   - Differential Abundance
   - Network Analysis
4. Check the **Citation** tab for method references.

## Project Structure

```text
SimpleMicrobiome/
  app.R
  modules/
    mod_top.R
    mod_fileload.R
    mod_preprocessing.R
    mod_barplot.R
    mod_taxa_comparison.R
    mod_alpha.R
    mod_beta.R
    mod_ancom.R
    mod_maaslin2.R
    mod_randomforest.R
    mod_sparcc.R
    mod_spieceasi.R
    mod_citation.R
  convert_service/
    app.py
    deploy/qiime-converter.service
  sample/
```

## Troubleshooting

- If loading fails, verify:
  - Delimiter (comma/tab) is consistent
  - Sample IDs match between abundance table and metadata
  - Taxa IDs match between abundance and taxonomy tables
  - Input files contain headers
- If a module returns empty results, check preprocessing filters and group sizes.

## Citation

Please cite this repository and the underlying method packages used in your analysis.
Module-wise references are available directly in the app’s **Citation** tab.

## License

This project is licensed under the MIT License.
See [LICENSE](LICENSE) for details.
