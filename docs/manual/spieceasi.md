# SpiecEasi Manual

## 1. What This Module Does
SpiecEasi infers sparse microbial association networks using graphical-model approaches.
Compared with correlation-style methods, it focuses on conditional dependency structure.

## 2. Left Panel Parameters (Detailed)

### 2.1 Group / Subgroup Filters
- Defines sample subset used for model fitting.
- Smaller subsets reduce stability of sparse graph estimation.

### 2.2 Taxonomic Level / Feature Scope
- Determines node count and dimensionality.
- High-dimensional low-sample settings can destabilize model selection.

### 2.3 Feature Filtering Controls
- Removes sparse/noisy taxa before model fit.
- Stronger filtering improves convergence/runtime.
- Over-filtering may remove biologically meaningful taxa.

### 2.4 Method / Model Controls (module-dependent labels)
- Controls graphical model estimation behavior and sparsity path.
- Conservative settings -> fewer edges, often more stable.
- Aggressive settings -> more edges, higher false-link risk.

### 2.5 Stability / Tuning Controls
- Parameters that affect edge selection robustness across tuning path.
- Stricter stability criteria typically reduce edge count but improve confidence.

### 2.6 Edge Threshold / Selection Rule (if exposed)
- Determines which inferred associations are retained in final network.
- Record this rule when exporting/reporting.

### 2.7 Visualization Controls
- Layout, label, and style parameters affect readability only.
- Do not interpret visual proximity as model strength unless encoded.

### 2.8 Export Controls
- Export graph/table outputs together with parameter settings for reproducibility.

## 3. How to Interpret SpiecEasi Results

### 3.1 Edge meaning
- Edges represent conditional dependency structure under model assumptions.
- This is stronger than simple correlation but still not causal proof.

### 3.2 Graph sparsity
- Very dense graph often indicates weak regularization/overfitting.
- Very sparse graph may indicate over-regularization or over-filtering.

### 3.3 Hub interpretation
- High-degree nodes can be key network elements, but check robustness across settings.

## 4. Practical Parameter Presets

### Preset A: Stable baseline
- Genus-level input
- Moderate prevalence filter
- Conservative sparsity/tuning setting
- Goal: convergent and interpretable baseline graph

### Preset B: Sensitivity analysis
- Re-run with slightly relaxed and stricter sparsity settings
- Keep all other settings fixed
- Compare edge overlap to assess robustness

## 5. Common Pitfalls
- Treating one parameter setting as final truth without sensitivity checks.
- Ignoring convergence/runtime warnings.
- Comparing graphs across runs with different filtering and interpreting as biological change.
- Over-interpreting weakly stable edges.

## 6. Recommended Reporting Items
- Sample subset and taxonomic level
- Filtering criteria
- Core model/tuning settings
- Final edge-selection rule
- Node/edge counts and robustness summary across sensitivity runs
