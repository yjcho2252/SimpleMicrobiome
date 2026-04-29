# SparCC Manual

## 1. What This Module Does
SparCC estimates correlation-like associations between taxa under compositional constraints.
It is typically used to build co-occurrence/co-exclusion style microbial networks.

## 2. Left Panel Parameters (Detailed)

### 2.1 Group / Subgroup Filters
- Defines which samples are included before network inference.
- Different subsets can produce very different edge structures.
- Recommendation: keep subset definition biologically coherent and sample size adequate.

### 2.2 Taxonomic Level
- Controls node granularity (e.g., ASV/Genus/Species depending on module settings).
- Lower level -> more nodes, more sparsity, noisier correlations.
- Higher level -> more stable, less specific network.

### 2.3 Pre-filter / Feature Inclusion Controls
- Removes low-prevalence or low-abundance taxa before SparCC.
- Higher filtering stringency:
  - Pros: less noise, faster runtime, cleaner network.
  - Cons: can remove weak but real biological interactions.

### 2.4 Correlation Threshold (|r| cutoff)
- Keeps only edges above absolute correlation strength.
- Lower cutoff -> dense network (“hairball” risk).
- Higher cutoff -> sparse network (clearer, but may miss moderate links).

### 2.5 Significance Threshold (if enabled)
- Filters edges by p-value support.
- Lower p-cutoff (stricter) improves confidence but reduces edge count.

### 2.6 Positive/Negative Edge Display (if available)
- Controls whether co-occurrence (+), co-exclusion (-), or both are shown.
- Important for biological interpretation of interaction directionality.

### 2.7 Network Visualization Controls
- Node/edge size, labels, and layout controls mainly affect readability.
- These do not change inferred associations.

### 2.8 Plot/Export Controls
- Export options should be used with threshold settings recorded.
- Always save threshold values with exported network for reproducibility.

## 3. How to Interpret SparCC Results

### 3.1 Edge sign
- Positive edge: taxa tend to vary together.
- Negative edge: taxa tend to vary inversely.
- Note: correlation-like association is not causal interaction.

### 3.2 Edge strength
- Larger |r| suggests stronger association pattern under current preprocessing/subset.
- Compare strengths within the same run settings only.

### 3.3 Node connectivity
- Highly connected taxa can be ecological hubs, but may also reflect compositional or prevalence artifacts.

## 4. Practical Parameter Presets

### Preset A: Exploratory network
- Moderate prevalence filter
- Moderate |r| threshold
- No overly strict p-cutoff initially
- Goal: map overall structure

### Preset B: Report-ready network
- Stronger prevalence filter
- Higher |r| threshold
- Significant edges only
- Goal: robust, interpretable core network

## 5. Common Pitfalls
- Treating network edges as direct biological interactions.
- Comparing networks built with different thresholds as if directly equivalent.
- Over-interpreting hubs from very sparse/imbalanced datasets.
- Ignoring sample-size reduction from subgroup filtering.

## 6. Recommended Reporting Items
- Sample subset definition
- Taxonomic level
- Pre-filter criteria
- Correlation and significance thresholds
- Number of nodes/edges and positive/negative edge counts
