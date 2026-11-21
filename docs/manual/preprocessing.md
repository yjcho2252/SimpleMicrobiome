# Preprocessing Manual

## 1. What This Module Does
Applies sample/taxa filtering and selection before all downstream analyses.

## 2. Left Panel Parameters

### 2.1 Sample filtering controls
- Removes samples by quality/coverage criteria.
- Directly changes group balance and statistical power.

### 2.2 Rarefaction
- Turn on **Enable rarefaction** to use the lowest read count among the selected samples as the target depth. The option is off by default.
- Reads are randomly drawn without replacement to the displayed target depth.
- Deselect samples with zero or invalid read counts before applying rarefaction.
- Taxa with zero total reads are removed from the resulting dataset.

### 2.3 Selection state controls
- Determines what subset is carried into downstream modules.
- Sample selection and rarefaction settings apply when you switch to an analysis tab.
- The selection summary box reports counts as `Total: N samples` and `Selected: M samples`.

## 3. Recommended Workflow
1. Set the minimum read count and select the samples to analyze.
2. If needed, enable rarefaction and check the automatically calculated depth.
3. Confirm the sample and group counts, then open an analysis tab.
