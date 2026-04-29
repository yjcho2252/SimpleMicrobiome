# Random Forest Manual

## 1. What This Module Does
Builds predictive models from microbiome features for:
- **Classification** (categorical outcome)
- **Regression** (continuous outcome)

Main outputs include:
- Model status/summary
- Performance metrics (classification or regression)
- SHAP importance views
- Permutation importance views

## 2. Left Panel Parameters (Detailed)

### 2.1 Primary grouping variable
- Defines the main metadata context for filtering and subgroup workflow.
- Use a biologically meaningful grouping column.

### 2.2 Select subgroup
- Enables nested analysis mode.
- When enabled, model is trained on a restricted subset (primary level + secondary grouping context).

### 2.3 Primary level to include
- Active only when subgroup mode is on.
- Restricts samples to one level before modeling.
- Too narrow subset can cause unstable metrics.

### 2.4 Secondary grouping variable
- Often used as target context in subgroup scenarios.
- Check class distribution before training.

### 2.5 Outcome type
- `Auto`: inferred from target format.
- `Classification`: force class prediction.
- `Regression`: force continuous-value prediction.
- Recommendation: set explicitly unless outcome type is unambiguous.

### 2.6 Levels to include (optional)
- Restricts included outcome classes.
- Useful for focused binary tasks (e.g., Control vs Case).
- Excluding classes changes scientific scope.

### 2.7 SHAP target class
- Classification only.
- Chooses which class-specific SHAP explanation is shown.
- Important: SHAP ranking can differ by target class.

### 2.8 Taxonomic level
- Controls feature granularity (`ASV`, `Genus`, `Species`, etc. depending on data/module state).
- Lower rank increases dimensionality/sparsity.
- `Genus` is often a stable default for first-pass modeling.

### 2.9 Feature transform
- `TSS`: relative abundance style; intuitive baseline.
- `CLR`: compositional-aware transform; often preferred for robust inference.
- `Presence/Absence`: uses occurrence signal only.
- `log`: compresses extreme abundance and reduces skew.
- Recommendation: compare at least `TSS` and `CLR` for sensitivity.

### 2.10 Feature prevalence cutoff (0-20%)
- Removes sparse features before model fitting.
- Higher cutoff reduces noise and runtime.
- Too high can remove biologically useful rare taxa.

### 2.11 Train ratio
- Used in holdout mode.
- Higher value gives more training data but less reliable test estimate.
- Typical operating range: 0.75-0.85.

### 2.12 Top N features by mean abundance
- Preselects candidate feature space.
- Lower N reduces runtime and overfitting risk.
- Too small N may miss predictive taxa.

### 2.13 Number of trees (ntree)
- More trees typically stabilize model estimates.
- Runtime increases with ntree.
- Practical baseline: 500.

### 2.14 mtry (0 = auto)
- Number of features considered per split.
- `0` uses module default heuristic.
- Manual tuning can improve performance when feature dimension is high.

### 2.15 Validation mode
- `Holdout split`: faster, single-split estimate.
- `K-fold CV`: slower, more stable estimate of generalization.

### 2.16 K for K-fold CV
- Larger K can improve estimate stability but increases runtime.
- Common baseline: K=5.

### 2.17 Random seed
- Ensures reproducibility of splits and model randomness.
- Keep fixed for fair parameter comparisons.

### 2.18 Plot Width / Height / Base Font Size
- Controls readability of SHAP/importance figures and labels.

## 3. Validation and Metrics Interpretation

### 3.1 Classification metrics (typical)
- **Accuracy**: overall correct classification rate.
- **AUC/ROC** (if available): discrimination quality across thresholds.
- **Confusion pattern**: identifies which classes are commonly confused.

Interpretation notes:
- High accuracy with imbalanced classes can be misleading.
- Always inspect class-level behavior, not only one aggregate metric.

### 3.2 Regression metrics (typical)
- **R²**: explained variance proportion.
- **RMSE/MAE**: prediction error magnitude.

Interpretation notes:
- Compare errors to practical/biological scale of target variable.
- R² alone is not enough when target range is narrow.

### 3.3 Holdout vs K-fold
- Holdout may vary substantially by split.
- K-fold is preferred when sample size is limited or class balance is fragile.

## 4. SHAP Interpretation Guide

### 4.1 What SHAP shows
- Feature contribution magnitude and direction for model predictions.
- Helps explain model behavior, not causality.

### 4.2 SHAP target class (classification)
- SHAP values are class-conditional in multi-class settings.
- Always report which class SHAP is based on.

### 4.3 Reading SHAP bar importance
- Higher absolute SHAP importance = stronger average model contribution.
- Use as ranking signal, then cross-check biological plausibility.

### 4.4 Common SHAP pitfalls
- Treating SHAP ranking as causal effect ranking.
- Over-interpreting very small differences among top features.
- Ignoring class-conditional differences.

## 5. Permutation Importance Interpretation

### 5.1 What it measures
- Performance drop when a feature is permuted.
- Larger drop indicates stronger predictive dependency.

### 5.2 Practical use
- Compare with SHAP top features.
- Stable overlap between SHAP and permutation importance increases confidence.

### 5.3 Pitfalls
- Correlated features can dilute individual permutation importance.
- Low importance does not always mean biological irrelevance.

## 6. Recommended Workflows

### Preset A: Fast baseline
- Tax level: `Genus`
- Transform: `TSS`
- Prevalence cutoff: `5%`
- Top N: `100`
- ntree: `500`
- Validation: `Holdout`, train ratio `0.8`
- Seed fixed

### Preset B: Robust report run
- Run both `TSS` and `CLR`
- Validation: `K-fold CV` (K=5)
- Keep seed fixed
- Compare metric consistency + SHAP overlap

### Preset C: Sparse data stabilization
- Increase prevalence cutoff moderately
- Reduce Top N feature space
- Use K-fold for more stable estimate

## 7. Common Failure Patterns and Fixes
- **Overfitting (train high, test low)**:
  - Reduce feature space (Top N), increase prevalence cutoff, use CV.
- **Unstable runs across repeats**:
  - Fix seed, move from holdout to K-fold.
- **Poor minority-class behavior**:
  - Restrict levels for focused comparison or rebalance design.
- **Slow runtime**:
  - Lower Top N, reduce ntree moderately, simplify subgroup filters.

## 8. Reporting Checklist
Include:
- Outcome type and included levels
- Transform method and taxonomic level
- Prevalence cutoff and Top N features
- Validation mode and seed
- Key metrics (with class context if classification)
- SHAP target class and top feature summary
