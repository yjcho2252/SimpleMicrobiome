# Correlation Heatmap Manual

## 1. What This Module Does
Correlation Heatmap visualizes associations between taxa and metadata groups as a matrix.
The module supports abundance-based associations and prevalence-based binary associations.

## 2. Left Panel Parameters

### 2.1 Primary grouping variable
- Defines the metadata variable used for group-aware association summaries.
- In subgroup mode, the analysis can be restricted to one primary level before using the secondary grouping variable.

### 2.2 Select subgroup / primary level / secondary grouping variable
- `Select subgroup` restricts samples to one primary level.
- When enabled, the secondary grouping variable is used as the group variable for the heatmap context.
- Small subgroup sample sizes can make associations unstable.

### 2.3 Group type
- Options: `Auto`, `Continuous`, `Categorical`.
- `Auto` infers the group type from metadata values.
- Continuous and categorical groups are interpreted differently in association calculations and legends.

### 2.4 Taxonomic level
- Options: `ASV`, `Genus`, `Species`, `Strain`.
- Lower levels are more specific and sparse.
- Higher levels are usually more stable and readable.

### 2.5 Association mode
- `Abundance-based`: uses abundance values.
- `Prevalence-based (binary)`: binarizes taxa to presence/absence before computing associations.
- Prevalence mode is useful when detection patterns are more relevant than abundance magnitude.

### 2.6 Heatmap value scale
- `Raw`: shows the association values on the original module scale.
- `Z-score (by taxa)`: standardizes values by taxa to emphasize relative patterns.
- Do not compare color intensity across heatmaps unless value scale and settings are the same.

### 2.7 Advanced options
- `Data transform`: `TSS` or `CLR` for abundance-based calculations.
- `Correlation method`: `Pearson` or `Spearman` for abundance-based calculations.
- `Prevalence filter (%)`: removes low-prevalence taxa before heatmap construction.
- `Max taxa`: limits the number of taxa shown.
- `Selected taxa (optional)`: if selected, only those taxa are shown.
- `Show row dendrogram` and `Show column dendrogram`: controls clustering display.
- `Color scale abs limit (0 = auto)`: fixes the symmetric color scale limit when nonzero.
- `Mask non-significant cells (FDR)`: hides cells that do not pass the selected FDR cutoff.
- `FDR cutoff (q)`: q-value threshold used when masking is enabled.

### 2.8 Plot dimensions
- `Cell width`, `Cell height`, and `Base Font Size` control readability.
- Increase cell size for long labels or dense matrices.

## 3. Result Interpretation

### 3.1 Abundance-based mode
- Values summarize abundance-based association direction and magnitude.
- Positive values indicate positive association.
- Negative values indicate inverse association.
- The selected transform and correlation method affect the scale and interpretation.

### 3.2 Prevalence-based mode
- Taxa are converted to binary present/absent values.
- For categorical groups, values reflect directional prevalence differences between group contexts.
- This mode answers detection-pattern questions rather than abundance-magnitude questions.

### 3.3 FDR masking
- When enabled, non-significant cells are masked based on the selected q-value cutoff.
- Masking improves readability but can hide borderline cells.
- Report the q cutoff when exporting or presenting masked heatmaps.

### 3.4 Dendrograms
- Row/column dendrograms cluster similar association profiles.
- Clustering is descriptive and depends on the selected values and scale.

## 4. Practical Presets

### Preset A: Abundance association scan
- Association mode: `Abundance-based`
- Transform: `CLR`
- Correlation method: `Pearson` or `Spearman`
- Value scale: `Raw`
- FDR masking: off for first-pass exploration

### Preset B: Detection-pattern scan
- Association mode: `Prevalence-based (binary)`
- Moderate prevalence filter
- Value scale: `Raw`
- FDR masking: optional for report figures

### Preset C: Report-ready matrix
- Curated selected taxa
- Fixed color scale abs limit
- FDR masking on with reported cutoff
- Larger cell dimensions for label readability

## 5. Common Pitfalls
- Comparing colors across heatmaps with different value scales or transforms.
- Treating masked cells as zero values; masked means not displayed under the selected FDR rule.
- Overloading too many taxa into one matrix.
- Interpreting association as causation.

## 6. Recommended Reporting Items
- Primary/subgroup setup
- Group type
- Taxonomic level
- Association mode
- Transform and correlation method when abundance-based
- Prevalence filter and max taxa
- Value scale and color scale limit
- Whether FDR masking was used and the q cutoff
