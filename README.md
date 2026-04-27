# SimpleMicrobiome

[![R](https://img.shields.io/badge/R-%3E%3D4.2-276DC3)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-App-1F77B4)](https://shiny.posit.co/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

SimpleMicrobiome is an R Shiny application for end-to-end microbiome analysis with a GUI.  
It integrates data loading, preprocessing, taxa profiling, diversity, differential abundance, network analysis, and association analysis in a single workflow.
Public app: https://simplemicrobiome.mglab.org/

## About the App

SimpleMicrobiome is developed and maintained by the Microbial Genomics Laboratory at Kangwon National University in Korea.

## Current App Modules

- **Top**
  - App overview, quick workflow and version history
- **Data Loading**
  - Upload ASV/OTU table, taxonomy table, metadata
  - Load/download bundled example data
  - external QZA converter link
- **Preprocessing**
  - Sample/taxa filtering before downstream analyses
- **Taxa Profiles**
  - Taxa Bar plot
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

## How to Use

SimpleMicrobiome supports two usage modes:

- Use the hosted web app (no local installation required)
- Run locally (requires dependency installation)

### Use on Web

- Open: https://simplemicrobiome.mglab.org/
- Upload your input files in **Data Loading**
- Continue analysis through each module

### Run Locally

#### Prerequisites

- R (recommended: 4.2+)
- Internet access for package installation
- System build tools may be required for some packages:
  - Windows: Rtools
  - macOS: Xcode Command Line Tools
  - Linux: compiler/system libraries

#### 1) Clone the repository

```bash
git clone https://github.com/yjcho2252/SimpleMicrobiome.git
cd SimpleMicrobiome
```

#### 2) Install dependencies

Run the installer script from the project root:

```r
source("install.R")
```

#### 3) Run the app

```r
shiny::runApp("app.R")
```

If startup fails, rerun `source("install.R")` and check package/system dependency errors in the console.

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
- For local mode, ensure `source("install.R")` completes without missing-package errors.

## License

This project is licensed under the MIT License.  
See [LICENSE](LICENSE) for details.
