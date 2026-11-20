# Preprocessing Manual

## 1. What This Module Does
Applies sample/taxa filtering and selection before all downstream analyses.

## 2. Left Panel Parameters

### 2.1 Sample filtering controls
- Removes samples by quality/coverage criteria.
- Directly changes group balance and statistical power.

### 2.2 Subsampling
- Turn on **Enable subsampling** and enter a positive whole number in **Reads per sample**. The option is off by default, and the initial depth is 10,000.
- Each selected sample is randomly subsampled without replacement to the target depth. For example, samples with 15,000, 30,000, and 40,000 reads each yield 10,000 reads at a target of 10,000. The app uses a fixed seed of 42 for reproducibility.
- If any selected sample has fewer reads than the target, lower the depth or deselect those samples before continuing.
- Taxa with zero total reads are removed from the resulting dataset.

### 2.3 Selection state controls
- Determines what subset is carried into downstream modules.
- Sample selection and subsampling settings apply when you switch to an analysis tab.
- The selection summary box reports counts as `Total: N samples` and `Selected: M samples`.

## 3. Recommended Workflow
1. Set the minimum read count and select the samples to analyze.
2. If needed, enable subsampling and choose a depth supported by all selected samples.
3. Confirm the sample and group counts, then open an analysis tab.
