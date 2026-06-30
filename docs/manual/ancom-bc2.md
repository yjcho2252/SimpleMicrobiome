# ANCOM-BC2 Manual

## 1. What This Module Does
ANCOM-BC2 performs differential abundance (DA) testing while accounting for compositional bias.
This module compares selected metadata levels and reports taxon-level log2 fold-change, p-values, q-values, and ANCOM-BC2 differential-abundance flags.

## 2. Left Panel Parameters

### 2.1 Primary grouping variable
- Defines the main metadata variable used for comparison.
- If subgroup mode is off, this variable is also the comparison variable.

### 2.2 Select subgroup / primary level to include
- Restricts the analysis to one selected primary level before fitting ANCOM-BC2.
- When subgroup mode is enabled, the selected secondary grouping variable becomes the comparison variable.
- Narrow subsets reduce sample count and can destabilize model estimates.

### 2.3 Comparison groups (levels)
- Selects the group levels included in the ANCOM-BC2 model.
- At least one non-reference level is required.
- Up to five levels can be selected in the current UI.

### 2.4 Reference level
- Defines the baseline group for log2FC contrasts.
- Positive log2FC means the taxon is increased in the comparison level relative to the reference level.
- Negative log2FC means the taxon is decreased in the comparison level relative to the reference level.

### 2.5 Taxonomic level
- Options: `ASV`, `Genus`, `Species`, `Strain`.
- Higher ranks are usually more stable but less specific.
- Lower ranks are more specific but often sparser.

### 2.6 Statistical metric
- Options: `q-value (FDR)` and `p-value`.
- This controls the volcano plot y-axis and the numeric labels shown on the bar plot.
- It does not change ANCOM-BC2 fitting. ANCOM-BC2 is run with `alpha = 0.05`.

### 2.7 Prevalence filter
- Default: `5%`.
- Allowed UI range: `0-20%`.
- Taxa are retained only when prevalence is greater than the selected threshold in the currently included sample subset.
- Prevalence filtering is applied before ANCOM-BC2 fitting.

### 2.8 Zero-variance taxa removal
- Applied after prevalence filtering and before model fitting.
- Taxa with no remaining variation in the selected subset are removed.

### 2.9 Structural zero handling
- Enabled in ANCOM-BC2 through `struc_zero = TRUE`.
- The module also uses `neg_lb = TRUE`.
- Structural zero handling addresses group-specific all-zero patterns during ANCOM-BC2 fitting.

### 2.10 Advanced options
- `Additional covariates for fix_formula`: metadata variables added as fixed effects.
- `Interaction terms for fix_formula`: selected interaction terms among fixed effects.
- `Random-effect grouping variables`: random intercept terms written as `(1|variable)`.
- `Enable ordered trend test`: runs ANCOM-BC2 trend testing only when enabled and at least three levels are selected.
- The formula preview shows the fixed and random model terms that will be used.

## 3. Result Tabs

### 3.1 Volcano plot
- x-axis: log2FC for the selected contrast against the reference level.
- y-axis: `-log10(q-value)` or `-log10(p-value)`, depending on the selected statistical metric.
- Dashed vertical lines mark `|log2FC| = 0.5`.
- Dashed horizontal line marks `p/q = 0.05`.

### 3.2 Bar plot
- Bars show the top 10 taxa ranked by absolute log2FC.
- The primary selection is ANCOM-BC2 `diff = TRUE`, which corresponds to taxa significant at `q < alpha` from ANCOM-BC2. In this module `alpha = 0.05`.
- If no `diff = TRUE` taxa are available, the plot falls back to the top 10 taxa by absolute finite log2FC.
- The number printed beside each bar is the selected statistical metric (`q-value` or `p-value`).

### 3.3 Table
- Reports selected ANCOM-BC2 coefficient estimates, standard errors, W statistics, p-values, q-values, differential flags, contrast labels, and taxonomy fields.

## 4. Result Columns

### 4.1 feature_id
- The tested taxon or feature ID at the selected taxonomic level.

### 4.2 contrast
- The non-reference group level being compared against the reference level.

### 4.3 lfc
- Log2 fold-change for the contrast.
- Positive values indicate increase in the comparison group versus the reference group.
- Negative values indicate decrease in the comparison group versus the reference group.

### 4.4 se
- Standard error of the log2FC estimate.

### 4.5 W
- ANCOM-BC2 test statistic.

### 4.6 p_val
- Raw p-value.

### 4.7 q_val
- Multiple-testing adjusted p-value.
- This is the preferred significance value for reporting.

### 4.8 diff
- ANCOM-BC2 differential-abundance flag.
- In ANCOM-BC2 output this is based on whether the adjusted value is below `alpha`.

## 5. Interpretation Workflow
1. Confirm the comparison variable, selected levels, and reference level.
2. Check sample counts in the status box.
3. Prioritize `diff = TRUE` and q-value significant taxa.
4. Interpret log2FC direction relative to the reference level.
5. Check whether results are robust to prevalence cutoff, covariates, and taxonomic level.

## 6. Recommended Reporting Items
- Primary/subgroup setup and comparison variable
- Included levels and reference level
- Taxonomic level
- Prevalence cutoff
- Fixed effects, interactions, and random effects
- Whether ordered trend testing was enabled
- ANCOM-BC2 `alpha`
- Significant taxa with log2FC direction, q-values, and contrast labels
