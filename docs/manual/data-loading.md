# Data Loading Manual

## 1. What This Module Does
Loads all base inputs required by downstream modules.

## 2. Left Panel Parameters

### 2.1 Feature table upload
- Expected role: sample-by-feature abundance/count matrix.
- If orientation/header is wrong, most downstream modules fail silently or show empty selectors.

### 2.2 Taxonomy table upload
- Maps feature IDs to taxonomic ranks.
- Determines whether modules can display/test at Genus/Species/etc.
- Current modules support taxonomy selection up to `Strain` where available.

### 2.3 Metadata upload
- Defines grouping/covariate columns used throughout app.
- Column types (categorical vs numeric) affect available analyses and model behavior.

### 2.4 Example data load/download
- Selects and loads a known-valid demo set to validate app state quickly.
- Available example datasets include the HMP V3-V5 body-site cohort and the SprockettTH multi-age fecal cohort.
- `Load Example` loads the currently selected example dataset into the app.
- `Download Example` exports the currently selected example dataset as a zip file.
- Useful to isolate whether issues are data-specific.
- While the data are being parsed and validated, the file-load status area shows a waiting message and a spinner only during active processing.

### 2.5 QZA converter link/entry
- Path for converting QIIME2 artifact formats before import.
- Required when raw input is not app-native tabular format.

## 3. Validation Checklist
- Sample IDs overlap across feature table and metadata.
- Feature IDs overlap between feature and taxonomy tables.
- Metadata grouping columns have expected levels and limited missing values.
