# Preprocessing Manual

## 1. What This Module Does
Applies sample/taxa filtering and selection before all downstream analyses.

## 2. Left Panel Parameters

### 2.1 Sample filtering controls
- Removes samples by quality/coverage criteria.
- Directly changes group balance and statistical power.

### 2.2 Taxa filtering controls
- Removes sparse or low-abundance taxa.
- Improves runtime/stability but can remove biologically meaningful rare taxa.

### 2.3 Selection state controls
- Determines what subset is carried into downstream modules.
- Changing selection here immediately changes results elsewhere.
- The selection summary box reports counts as `Total: N samples` and `Selected: M samples`.

## 3. Recommended Workflow
1. Apply moderate sample filter first.
2. Apply taxa prevalence/abundance filters conservatively.
3. Confirm remaining sample/group counts before DA/modeling.
