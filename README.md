# SimpleMicrobiome

<p align="center">
  <img src="www/260617_SimpleMicrobiome_icon.png" alt="SimpleMicrobiome icon" width="200" />
</p>

[![R](https://img.shields.io/badge/R-4.5.1-276DC3)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-App-1F77B4)](https://shiny.posit.co/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

SimpleMicrobiome is an R Shiny application for downstream microbiome analysis with a GUI.  
It integrates data loading, preprocessing, taxa profiling, diversity, differential abundance, network analysis, and association analysis in a single workflow.
Public app: https://simplemicrobiome.mglab.org/

## About the App

SimpleMicrobiome is developed and maintained by the Microbial Genomics Laboratory at Kangwon National University, Korea.

## Current App Modules

- **Top**
  - App overview, quick workflow and version history
- **Data Loading**
  - Upload ASV/OTU table, taxonomy table, metadata
  - Select, load, and download bundled example datasets
  - external QZA converter link
- **Preprocessing**
  - Sample/taxa filtering before downstream analyses
- **Taxa Profiles**
  - Taxa Bar plot
  - Taxa Comparison
- **Diversity**
  - Alpha Diversity
  - Beta Diversity
- **Differential Abundance**
  - ANCOM-BC2
  - MaAsLin2
  - Random Forest
- **Network**
  - SparCC
  - SPIEC-EASI
- **Association**
  - Correlation Heatmap
  - Association Biplot (dbRDA/CCA)
- **Citation**
  - Method/package references

## How to Use

SimpleMicrobiome supports three usage modes:

- Use the hosted web app (no local installation required)
- Run locally with the Docker image (recommended for local use)
- Install from source code (mainly for development or code modification)

### Use on Web

- Open: https://simplemicrobiome.mglab.org/
- Upload your input files in **Data Loading**
- Continue analysis through each module

### Run Locally with Docker (Recommended)

For local use, the Docker image is recommended because it avoids most R package
dependency conflicts that can occur when installing the latest CRAN,
Bioconductor, and GitHub packages directly on a local machine.

Docker image: https://hub.docker.com/r/mglab2252/simplemicrobiome

#### Prerequisites

- Docker Desktop for Windows or macOS, or Docker Engine for Linux
- Docker must be installed and running before executing the commands below.

The current Docker image targets the `linux/amd64` architecture. On Apple
Silicon Macs, Docker Desktop may run the image through amd64 emulation, which
can result in slower startup and execution.

#### 1) Pull the image

```bash
docker pull mglab2252/simplemicrobiome:v1.0
```

#### 2) Run the app

```bash
docker run --rm -p 3838:3838 mglab2252/simplemicrobiome:v1.0
```

#### 3) Open the app

Open your browser at:

```text
http://localhost:3838
```

If port 3838 is already in use, map another local port:

```bash
docker run --rm -p 3839:3838 mglab2252/simplemicrobiome:v1.0
```

Then open:

```text
http://localhost:3839
```

#### Build the image from GitHub

To build the same application image directly from the source repository:

```bash
git clone https://github.com/yjcho2252/SimpleMicrobiome.git
cd SimpleMicrobiome
docker build --pull -t simplemicrobiome:local .
docker run --rm -p 3838:3838 --name simplemicrobiome simplemicrobiome:local
```

Then open `http://localhost:3838`. The first build can take a long time because
the required CRAN, Bioconductor, and GitHub packages are compiled and installed.
Later builds reuse Docker layers while `install.R`, `renv.lock`, and the
dependency-install section of the Dockerfile are unchanged.

### Install from Source Code

#### Prerequisites

- R 4.5.1
- Internet access for package installation
- System build tools may be required for some packages:
  - Windows: Rtools
  - macOS: Xcode Command Line Tools
  - Linux: compiler/system libraries
- The standard `source("install.R")` workflow restores the pinned dependency
  set and does not require `clarabel` or a Rust toolchain. Rust (`rustc` and
  `cargo`) may be required only when packages are installed or updated outside
  the pinned `renv.lock` and a Rust-based package such as `clarabel` must be
  compiled from source.

#### 1) Clone the repository

```bash
git clone https://github.com/yjcho2252/SimpleMicrobiome.git
cd SimpleMicrobiome
```

#### 2) Install dependencies

Run the installer script from the project root. The script restores the tested
package versions and GitHub commits recorded in `renv.lock`:

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
- For Docker mode, check that Docker Desktop or Docker Engine is running and that the selected local port is not already in use.
- For source-code installation, ensure `source("install.R")` completes without missing-package errors.

## License

This project is licensed under the MIT License.  
See [LICENSE](LICENSE) for details.
