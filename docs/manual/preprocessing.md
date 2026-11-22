# Preprocessing Manual

## 1. What This Module Does
Applies sample/taxa filtering and selection before all downstream analyses.

## 2. Left Panel Parameters

### 2.1 Sample filtering controls
- Removes samples by quality/coverage criteria.
- Directly changes group balance and statistical power.

### 2.2 Selection state controls
- Sample selection applies when you switch to an analysis tab.
- The selected read counts are passed to downstream modules. Alpha Diversity, Jaccard beta diversity, and Taxa Bar Plot handle subsampling within their own modules.
- The selection summary box reports counts as `Total: N samples` and `Selected: M samples`.

## 3. Recommended Workflow
1. Set the minimum read count and select the samples to analyze.
2. Confirm the sample and group counts, then open an analysis tab.
