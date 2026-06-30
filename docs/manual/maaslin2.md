# MaAsLin2 Manual

## 1. What This Module Does
MaAsLin2 fits multivariable association models between microbial features and metadata variables.
This module is designed for selected group-level contrasts with optional fixed covariates, interaction terms, and random effects.

## 2. Left Panel Parameters

### 2.1 Primary grouping variable
- Defines the main metadata variable used for analysis.
- If subgroup mode is off, this variable is the comparison variable.

### 2.2 Select subgroup / primary level to include
- Restricts the analysis to a selected primary level before model fitting.
- When subgroup mode is enabled, the secondary grouping variable becomes the comparison variable.

### 2.3 Comparison groups (levels)
- Selects levels included in the model.
- The reference level must be one of the selected levels.

### 2.4 Reference level
- Defines the baseline level for model contrasts.
- Positive log2FC means increased in the comparison level versus the reference level.
- Negative log2FC means decreased in the comparison level versus the reference level.

### 2.5 Taxonomic level
- Options: `ASV`, `Genus`, `Species`, `Strain`.
- Higher ranks are usually more stable; lower ranks are more specific and sparse.

### 2.6 Statistical metric
- Options: `q-value (FDR)` and `p-value`.
- This controls the volcano plot y-axis and the numeric labels shown on the bar plot.

### 2.7 Prevalence filter
- Default: `5%`.
- Allowed UI range: `0-20%`.
- Taxa are retained only when prevalence is greater than the selected threshold in the included sample subset.
- Zero-variance taxa are removed after prevalence filtering.

### 2.8 Analysis method and automatic normalization/transform
- Options: `LM` and `ZINB`.
- The module does not expose normalization and transform as separate user controls.
- They are selected automatically from the analysis method:

| Analysis method | MaAsLin2 normalization | MaAsLin2 transform |
|---|---|---|
| `LM` | `TSS` | `LOG` |
| `ZINB` | `NONE` | `NONE` |

- These settings affect coefficient scale and should be reported with MaAsLin2 results.

### 2.9 Fixed effects, interactions, and random effects
- `Additional covariates for fixed effects`: metadata variables added to the fixed-effect model.
- `Interaction terms for fixed effects`: selected interaction terms are expanded into derived metadata columns before fitting.
- `Random effects`: metadata grouping variables passed to MaAsLin2 as random effects.
- The model formula preview shows the applied fixed effects, selected interactions, and random effects.

## 3. Result Tabs

### 3.1 Volcano plot
- x-axis: log2FC for the selected comparison level against the reference level.
- y-axis: `-log10(q-value)` or `-log10(p-value)`, depending on the selected statistical metric.
- Positive log2FC means increased in the comparison level versus the reference level.

### 3.2 Bar plot
- Bars show the top 10 significant taxa ranked by absolute log2FC.
- If no significant rows are available, the plot falls back to the top 10 finite log2FC rows.
- The number printed beside each bar is the selected statistical metric (`q-value` or `p-value`).

### 3.3 Table
- Reports MaAsLin2 coefficient estimates, transformed effect sizes, significance values, selected contrast labels, and taxonomy fields.

## 4. Result Columns

### 4.1 Feature / taxon
- Microbial feature tested in the model.

### 4.2 Metadata / variable
- Predictor associated with the reported coefficient.
- With covariates, the coefficient is conditional on the other model terms.

### 4.3 Coefficient / log2FC-style effect
- Direction is relative to the reference level for group contrasts.
- Magnitude depends on the automatic normalization/transform and selected analysis method.

### 4.4 p-value
- Unadjusted evidence level.
- Useful for screening but not sufficient alone for high-dimensional reporting.

### 4.5 q-value
- Multiple-testing adjusted p-value.
- Preferred significance value for reporting.

## 5. Interpretation Workflow
1. Confirm comparison variable, selected levels, and reference level.
2. Confirm analysis method and the automatic normalization/transform it implies.
3. Check sample counts and prevalence filtering.
4. Interpret coefficient direction relative to the reference level.
5. Prioritize q-value significant associations.
6. Review covariates, interactions, and random effects for biological and statistical plausibility.

## 6. Recommended Reporting Items
- Primary/subgroup setup and comparison variable
- Included levels and reference level
- Taxonomic level
- Analysis method
- Automatic normalization and transform
- Prevalence cutoff
- Fixed effects, interactions, and random effects
- Significant associations with coefficient direction and q-values
