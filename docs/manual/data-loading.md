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

### 2.5 Bundled example datasets

The bundled example datasets are compact, app-ready subsets. They are intended for testing the workflow and learning the modules, not as complete reproductions of the original studies.

#### HMP V3-V5 body-site cohort
- Source: Human Microbiome Project V3-V5 16S rRNA gene data from Qiita study ID 1928.
- Current bundled subset: 80 samples and 2,895 features.
- Biological scope: oral, skin, and stool body habitats, with body-site labels including buccal mucosa, supragingival plaque, tongue dorsum, anterior nares, right antecubital fossa, and stool.
- Metadata columns: `SampleID`, `body_habitat`, `bodysite`, and `host_subject_id`.
- Useful grouping variables: `body_habitat` for broad oral/skin/stool comparisons and `bodysite` for finer body-site comparisons.
- Preparation note: the subset was processed outside SimpleMicrobiome from QIIME 2/DADA2-derived outputs and converted into app-ready abundance, taxonomy, and metadata tables.
- References: Human Microbiome Project Consortium, 2012; Gonzalez et al., 2018.

#### SprockettTH multi-age fecal cohort
- Source: `microbiomeDataSets::SprockettTHData()`.
- Original dataset scope: 16S rRNA gene profiling data from Bolivia, Finland, and Bangladesh, including adults, children, and infants with longitudinal sampling.
- Current bundled subset: Bolivia, Finland, and Bangladesh; 20 subjects per country; 2 fecal/stool timepoints per subject; 120 samples and 1,276 retained features.
- Metadata columns include `Country`, `Subject_ID`, `Sex`, `Age_Days`, `Age_Months`, `Age_Years`, `Age_Group`, `Age_Class`, `Feeding_Status_Corrected`, `Delivery_Mode`, `Community`, `Cohort`, and `Sample_Collection_Date`.
- Useful grouping variables: `Country` for between-country contrasts, `Age_Class` or `Age_Group` for age-related analyses, and `Subject_ID` for paired or repeated-measures workflows where supported.
- References: Sprockett et al., 2020; Subramanian et al., 2014; Vatanen et al., 2016; `microbiomeDataSets` SprockettTHData reference page.

### 2.6 QZA converter link/entry
- Path for converting QIIME2 artifact formats before import.
- Required when raw input is not app-native tabular format.

## 3. Validation Checklist
- Sample IDs overlap across feature table and metadata.
- Feature IDs overlap between feature and taxonomy tables.
- Metadata grouping columns have expected levels and limited missing values.

## 4. References
- Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, et al. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology. 37: 852-857.
- Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, et al. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods. 13: 581-583.
- Gonzalez A, Navas-Molina JA, Kosciolek T, McDonald D, Vazquez-Baeza Y, et al. 2018. Qiita: rapid, web-enabled microbiome meta-analysis. Nature Methods. 15: 796-798.
- Human Microbiome Project Consortium. 2012. Structure, function and diversity of the healthy human microbiome. Nature. 486: 207-214.
- Sprockett DD, Martin M, Costello EK, et al. 2020. Microbiota assembly, structure, and dynamics among Tsimane horticulturalists of the Bolivian Amazon. Nature Communications. 11: 3772.
- Subramanian S, Huq S, Yatsunenko T, et al. 2014. Persistent gut microbiota immaturity in malnourished Bangladeshi children. Nature. 510: 417-421.
- Vatanen T, Kostic AD, d'Hennezel E, et al. 2016. Variation in microbiome LPS immunogenicity contributes to autoimmunity in humans. Cell. 165: 842-853.
- microbiomeDataSets SprockettTHData reference: https://microbiome.github.io/microbiomeDataSets/reference/SprockettTHData.html
