# Beta Diversity Manual

## 1. What This Module Does
Beta Diversity evaluates **between-sample dissimilarity** and visualizes sample structure in reduced dimensions (ordination).
It is used to assess whether groups differ in community composition.

## 2. Left Panel Parameters (Detailed)

### 2.1 Distance Method
The distance choice is the most important modeling decision in this module.

#### Bray-Curtis
- **What it captures**: abundance-based compositional difference.
- **Sensitive to**: medium/high abundance taxa.
- **Less sensitive to**: shared absences.
- **Use when**: abundance shifts are biologically central.

#### Jaccard
- **What it captures**: presence/absence turnover.
- **Ignores**: abundance magnitude.
- **Use when**: detection pattern matters more than abundance scale.

#### Jensen-Shannon
- **What it captures**: divergence between relative abundance distributions.
- **Sensitive to**: global distribution-shape differences across taxa.
- **Use when**: profile-level compositional pattern differences are the main focus.

#### Aitchison (CLR log(x+1) pseudocount)
- **What it captures**: log-ratio geometry for compositional data.
- **Requires/assumes**: CLR-type compositional treatment.
- **Use when**: compositional interpretation is prioritized.

### 2.2 Primary Variable (Color)
- **What it controls**: Main group color mapping and interpretation axis.
- **Effect of setting**:
  - Defines the primary visual separation narrative.
  - Too many levels reduce readability.
- **Current behavior**:
  - `SampleID` can be selected.
  - If non-`SampleID` variables exist, the default is a non-`SampleID` variable.

### 2.3 Secondary Variable (Shape)
- **What it controls**: Shape encoding for subgroup context.
- **Effect of setting**:
  - Adds another layer of interpretation without forcing facets.
  - Too many shape levels can become visually noisy.

### 2.4 Taxonomic Level
- **Options**: `ASV`, `Species`, `Genus`.
- **Effect of setting**:
  - `ASV`: highest resolution, often noisier/sparser.
  - `Genus`: more stable, often easier to interpret group trends.

### 2.5 Plot Settings
- **Dot Size**: larger points improve visibility but can hide overlap.
- **Show Dot Outline**: improves separation between overlapping groups.
- **Show Group Ellipses**: adds centroid/spread cues; helpful for group-level pattern reading.
- **Show Group Hulls**: draws a convex hull (`chull`) polygon per group in ordination space.
  - Hulls require at least 3 samples per group.
  - Hulls are useful for quick group footprint comparison, but are sensitive to outliers.
- **Show Sample Names**: useful for QC/outlier identification; can clutter dense plots.

### 2.6 Axis Limits
- Manual `x/y min/max` can standardize view across runs.
- Leave blank for auto scale.
- Use fixed limits when comparing multiple plots side by side.

## 3. Post-Ordination Analyses

### 3.1 Clustering
Used to segment samples in ordination-derived space.

#### Clustering mode
- `Manual k`: you set the number of clusters directly.
- `Auto k (Silhouette)`: module searches a range and picks k by silhouette quality.

#### Manual k
- Best when you have hypothesis-driven expected cluster count.
- Risk: forcing k can produce biologically weak partitions.

#### Auto k min / max
- Defines search range for optimal k.
- Too narrow range can miss better solutions.
- Very wide range can increase runtime and unstable small-cluster solutions.

#### Show cluster labels on plot
- Adds cluster IDs to ordination view.
- Helpful for communication and downstream cluster-specific review.

#### Overlay cluster colors on dots
- Overlays clustering-result colors directly on ordination dots.
- Adds a `Cluster` legend category for visual cluster membership tracing.
- Useful for checking how algorithmic clusters align with metadata-colored structure.

#### Run Clustering / Download Cluster Assignment
- Executes clustering with current settings.
- Export table should be used to trace sample-to-cluster mapping in downstream interpretation.

#### Interpretation notes for clustering
- Clusters are algorithmic partitions, not automatic biological classes.
- Always compare cluster assignment with known metadata.
- Validate whether cluster structure persists across distance methods.
- `Average silhouette width` is the mean silhouette value across all samples (not a mean of per-cluster means).
- In auto mode, optimal `k` is selected by maximizing average silhouette width.

### 3.2 EnvFit
Fits metadata/taxa vectors onto ordination to identify directional associations.

#### EnvFit variable selector
- Chooses metadata variables for vector fitting.
- Include only biologically interpretable candidates.

#### EnvFit taxa selector
- Chooses taxa to project as vectors.
- Large taxa lists can clutter figure; prefer focused selection.

#### Show only significant vectors
- Displays only vectors passing significance threshold.
- Improves readability and reduces false emphasis.

#### p-value cutoff
- Threshold for vector retention (default commonly 0.05).
- Lower cutoff = stricter, fewer vectors.

#### Overlay EnvFit vectors on plots
- Adds vectors to PCoA/NMDS panels.
- Useful for directionality interpretation, but can over-clutter dense plots.

#### Run EnvFit
- Recomputes vector fitting with current selections.

#### Interpretation notes for EnvFit
- Vector direction indicates association gradient.
- Vector length reflects relative strength in ordination space.
- Significant vectors do not imply causality.

## 4. How Parameter Choices Change Interpretation

### Scenario A: Separation in Bray, weak in Jaccard
- Suggests abundance magnitude drives group differences more than simple occurrence turnover.

### Scenario B: Strong clustering but weak metadata alignment
- Algorithmic partition exists, but biological meaning may be limited.

### Scenario C: EnvFit shows many vectors, most near threshold
- Prioritize robust vectors (lower p-cutoff or repeated confirmation) before strong claims.

## 5. Common Pitfalls
- Treating clustering output as definitive biological subgrouping.
- Over-interpreting EnvFit arrows without checking significance and consistency.
- Ignoring method sensitivity (distance choice can change conclusions).
- Interpreting hull area alone as statistical separation; hulls are descriptive overlays, not hypothesis tests.

## 6. Minimal Recommended Presets

### Preset A: Exploratory baseline
- Distance: `Bray-Curtis`
- Primary color variable only
- Ellipses on, hulls optional, sample names off
- No manual axis limits

### Preset B: Cluster-focused
- Clustering mode: `Auto k (Silhouette)`
- k range: moderate (e.g., 2-10)
- Export cluster assignment table

### Preset C: EnvFit-focused
- Select limited metadata + taxa candidates
- `Show only significant vectors` on
- Keep p-cutoff conservative for report figures
