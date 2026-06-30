# SpiecEasi Manual

## 1. What This Module Does
SpiecEasi infers sparse microbial association networks using graphical-model approaches.
This module uses NetCoMi with `measure = "spieceasi"` to estimate taxon association networks.

Compared with SparCC, SpiecEasi focuses on sparse conditional-dependency structure rather than correlation-like association alone.

The module supports:
- `Single network`: estimate one or more group-specific networks.
- `Compare two groups`: estimate two group-specific networks and compare their edges.

## 2. Left Panel Parameters

### 2.1 Primary grouping variable
- Defines the metadata variable used for grouping samples.
- If subgroup mode is off, selected levels from this variable determine which group networks are estimated.

### 2.2 Select subgroup / primary level to include
- Restricts samples to a selected primary level before network inference.
- When subgroup mode is enabled, the secondary grouping variable becomes the grouping variable for network estimation.

### 2.3 Group levels
- Selects the group levels included in the run.
- In `Single network` mode, one or more levels can be selected.
- In `Compare two groups` mode, exactly two levels must be selected.

### 2.4 Analysis mode
- `Single network`: estimates networks for the selected eligible group levels.
- `Compare two groups`: estimates one SpiecEasi network per selected group, then compares shared and group-specific edges.
- Compare mode is descriptive. It does not run a separate edge-level statistical test between groups.

### 2.5 Taxonomic level
- Options: `ASV`, `Genus`, `Species`, `Strain`.
- Lower levels give more specific but sparser networks.
- Higher levels usually produce more stable and readable networks.

### 2.6 Prevalence filter
- Default: `10%`.
- Taxa below the selected prevalence threshold are removed before network inference.
- Stronger filtering can improve convergence and runtime but may remove rare taxa.

### 2.7 Node size and node color
- Node size options: `Connectivity`, `Abundance`.
- Node color is currently fixed to `None` in the UI.
- These options affect visualization, not the estimated graph.

### 2.8 Minimum absolute edge weight
- Default: `0.1`.
- Edges are shown in plots/tables only when `abs(weight) >= threshold`.
- This is an effect-size threshold for displayed edges, not a p-value threshold.

### 2.9 Max taxa for network
- Limits the number of taxa entering network estimation.
- Lower values reduce runtime and graph complexity.

### 2.10 Seed
- Sets the random seed used before network estimation for reproducibility.

## 3. Result Tabs

### 3.1 All
- Shows all nodes from the estimated network after edge filtering.
- Isolated nodes may remain visible.

### 3.2 Connected
- Shows the connected-node view after edge filtering.
- Isolated nodes are excluded.

### 3.3 Table
- Lists filtered edges from all estimated group networks.
- Edge weight sign indicates positive or negative association.

### 3.4 Summary
- Reports sample count, taxa count, group variable, nodes, edges, and connected components for each estimated group network.

### 3.5 Comparison Network
- Available when `Analysis mode = Compare two groups`.
- Visualizes edges classified as `Shared`, `Only <group 1>`, or `Only <group 2>`.

### 3.6 Differential Edges
- Available in compare mode.
- Reports edge-level differences between the two estimated networks:
  - edge endpoints
  - status (`Shared`, `Only group 1`, `Only group 2`)
  - group-specific weights
  - weight delta

### 3.7 Comparison Summary
- Summarizes node/edge counts, shared edges, unique edges, hub taxa, and modularity for the two groups.

### 3.8 Hub Table
- Lists hub candidates based on centrality-derived hub scores.
- The hub score combines standardized degree, betweenness, and eigenvector centrality.

## 4. Interpretation Guide

### 4.1 Edge meaning
- Edges represent sparse conditional-dependency structure under the SpiecEasi model.
- This is stronger than simple correlation but still not causal proof.

### 4.2 Graph sparsity
- Very dense graphs can indicate weak regularization or overfitting.
- Very sparse graphs can indicate over-filtering, insufficient sample size, or strong regularization.

### 4.3 Compare mode
- Compare mode compares two independently estimated group networks.
- `Only group` edges are edges retained in one group network but not the other after the current edge threshold.
- `weight delta` is descriptive and should not be interpreted as a formal differential-network p-value.

## 5. Common Pitfalls
- Treating one parameter setting as final truth without sensitivity checks.
- Comparing graphs across runs with different filtering or edge thresholds.
- Over-interpreting weak or unstable edges.
- Ignoring small group sample sizes. Groups with fewer than three samples are skipped.

## 6. Recommended Reporting Items
- Group/subgroup selection and included levels
- Analysis mode
- Taxonomic level
- Prevalence filter
- Max taxa for network
- Minimum absolute edge weight
- Number of samples, nodes, and edges by group
- For compare mode: shared edges, group-specific edges, and threshold used
