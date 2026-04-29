# MaAsLin2 Manual

## 1. What This Module Does
MaAsLin2 fits multivariable association models between microbial features and metadata variables.
It is suitable when adjustment covariates are required.

## 2. Left Panel Parameters (Detailed)

### 2.1 Primary target / grouping variable
- Defines main association focus.
- Must match the hypothesis being tested.

### 2.2 Select subgroup / primary level
- Restricts analysis to a specific stratum when needed.
- Helps reduce heterogeneity but can reduce power.

### 2.3 Secondary grouping / comparison levels
- Focuses analysis on selected level combinations.
- Ensure adequate sample size per included level.

### 2.4 Taxonomic level
- Chooses feature aggregation rank.
- Lower rank gives higher specificity but more sparsity.

### 2.5 Transform/normalization options
- Changes input scale before fitting.
- Coefficients and significance can vary by transform.

### 2.6 Covariate selection
- Adds adjustment variables for confounding control.
- Over-adjustment can reduce interpretability and power.

### 2.7 Significance/output filters
- Controls display/reporting thresholds.
- Prefer adjusted significance for final reporting.

## 3. Result Columns: How to Read
(Exact labels may vary by module version.)

### 3.1 Feature / Taxon
- Microbial feature tested in the model.

### 3.2 Metadata / Variable
- Predictor associated with the reported coefficient.
- In multi-covariate models, coefficient is conditional on other variables.

### 3.3 Coefficient
- Direction:
  - Positive: feature tends to increase with predictor (or compared level).
  - Negative: feature tends to decrease with predictor (or compared level).
- Magnitude is model-scale dependent (transform aware).

### 3.4 Standard error
- Precision of coefficient estimate.
- Large SE suggests unstable estimate or weak support.

### 3.5 p-value
- Unadjusted evidence level.
- Use primarily for screening, not final claim.

### 3.6 q-value / adjusted p-value
- Multiple-testing corrected significance.
- Recommended main significance criterion.

### 3.7 Model fit/status fields (if shown)
- Help identify convergence or estimation problems.
- Flagged rows should be interpreted cautiously.

## 4. Interpretation Workflow
1. Prioritize q-value significant associations.
2. Read coefficient direction and relative magnitude.
3. Check whether direction is consistent across related taxa levels.
4. Confirm robustness under a second transform setting.

## 5. Common Pitfalls
- Treating coefficients as unadjusted bivariate effects.
- Ignoring collinearity among covariates.
- Reporting only significance without effect direction.

## 6. Recommended Reporting Items
- Outcome/predictor/covariate definitions
- Transform and taxonomic level
- Correction method
- Significant associations with coefficient direction and q-values
- Any subgroup restrictions
