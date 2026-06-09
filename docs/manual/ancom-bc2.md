# ANCOM-BC2 Manual

## 1. What This Module Does
ANCOM-BC2 performs differential abundance (DA) testing while accounting for compositional bias.
It is used to identify taxa whose abundance differs between defined groups.

## 2. Left Panel Parameters (Detailed)

### 2.1 Primary grouping variable
- Defines the main biological comparison axis.
- All DA interpretation depends on this grouping context.

### 2.2 Select subgroup
- Enables stratified analysis mode.
- Restricts analysis to a selected primary level before testing.

### 2.3 Primary level to include
- Active in subgroup mode.
- Too strict subset can reduce power and destabilize estimates.

### 2.4 Comparison groups (levels)
- Chooses which levels are included in pairwise contrasts.
- Include only scientifically relevant levels to reduce noise/multiplicity.

### 2.5 Taxonomic level
- Sets aggregation rank for testing.
- Higher rank: more stable, less specific.
- Lower rank: more specific, more sparse.

### 2.6 Prevalence filter
- Default: `5%`.
- Allowed UI range: `0-20%`.
- What it does:
  - Keeps taxa whose prevalence is greater than the selected threshold within the currently included sample subset.
  - Prevalence is computed before model fitting on the relative-abundance table.
- Practical note:
  - Lower values keep more rare taxa but increase sparsity.
  - Higher values simplify the model but may remove biologically relevant low-prevalence taxa.

### 2.7 Zero-variance taxa removal
- Applied after prevalence filtering.
- What it does:
  - Removes taxa with no remaining variation in the selected subset before ANCOM-BC2 fitting.
- Why it matters:
  - Prevents degenerate taxa from entering the model.
  - Can reduce the number of tested features substantially in sparse subsets.

### 2.8 Structural zero handling
- Enabled in the module by default through ANCOM-BC2 fitting.
- What it does:
  - Treats taxa that are completely absent in a group as structural zeros rather than ordinary missing detections.
- Interpretation note:
  - This is different from prevalence filtering.
  - Prevalence filtering removes globally rare taxa first; structural zero handling addresses group-specific all-zero patterns during model fitting.

### 2.9 Statistical/model options
- Controls module-specific fitting and display behavior.
- Use consistent settings when comparing runs across groups.

### 2.10 Multiple testing correction
- Controls adjusted significance in multi-taxa testing context.
- Prefer corrected significance for reporting.

## 3. Result Columns: How to Read
(Exact labels can vary by module version.)

### 3.1 Taxon / Feature ID
- The tested taxonomic unit at current selected rank.

### 3.2 Effect size (e.g., coefficient / log-fold-like estimate)
- Direction:
  - Positive: relatively higher in reference contrast direction.
  - Negative: relatively lower in reference contrast direction.
- Magnitude reflects modeled difference size, not absolute abundance.

### 3.3 Standard error (SE)
- Uncertainty of effect estimate.
- Large SE indicates unstable estimate (often small group/sample support).

### 3.4 Test statistic / W statistic
- Strength of evidence under model assumptions.
- Use together with p-value/q-value, not in isolation.

### 3.5 p-value
- Unadjusted significance.
- Not sufficient alone for high-dimensional DA claims.

### 3.6 Adjusted p-value / q-value
- Multiplicity-corrected significance (recommended primary criterion).

### 3.7 Detection/significance flags (if shown)
- Convenience indicators for thresholded significance.
- Always verify corresponding adjusted value.

## 4. Interpretation Workflow
1. Filter by adjusted significance threshold.
2. Confirm the taxa passed prevalence filtering and zero-variance removal.
3. Check effect direction and magnitude.
4. Confirm biological plausibility at neighboring taxonomic levels.
5. Validate robustness with alternative filtering/transform context upstream.

## 5. Common Pitfalls
- Interpreting unadjusted p-values as final findings.
- Ignoring group imbalance after subgroup filtering.
- Overstating tiny effect sizes that are only marginally significant.

## 6. Recommended Reporting Items
- Grouping setup and included levels
- Taxonomic level
- Correction method
- Significant taxa list with effect direction and adjusted p-values
- Any subgroup restriction used
