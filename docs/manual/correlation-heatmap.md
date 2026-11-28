# Association Heatmap Manual

## 1. What This Module Does
Association Heatmap visualizes associations between taxa and metadata groups as a matrix.
The module supports abundance-based associations and prevalence-based binary associations.
Use the sample-level view to inspect individual samples, or the grouped view to compare groups and display significance marks.

## 2. Left Panel Parameters

### 2.1 Primary grouping variable
- Defines the metadata variable used for group-aware association summaries.
- All primary levels are included. In Individual Samples, each primary level defines a panel.

### 2.2 Secondary grouping variable
- Choose `None` to use the primary variable alone. In Individual Samples, a secondary variable divides each primary panel into labeled subpanels.
- In Group Summary, categorical secondary groups are compared within each primary panel. Panels share a color scale and secondary-column order.
- For a continuous secondary variable in Group Summary, the correlation pools samples with a recorded primary level.
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
- `Color palette`: selects a diverging palette, including a colorblind-friendly option.
- **Signed Color Scale Abs Limit (0 = Auto)** sets a symmetric color range from the negative limit to the positive limit. By default, Auto uses the 98th percentile of the absolute values of the displayed cells, with a fallback when that percentile cannot provide a positive limit. This setting applies to signed color scales.
- `Mask non-significant cells (FDR)`: hides cells that do not pass the selected FDR cutoff.
- An asterisk marks cells at or below the selected FDR cutoff. With masking enabled, cells above the cutoff or without a valid q-value are shown in gray.

### 2.8 Display Options
- `Plot Size` defaults to `Auto`. Choose `Manual` to set the width and height, and use `Base Font Size` to adjust text.

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
