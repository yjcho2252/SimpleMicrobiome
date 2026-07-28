# Association Biplot Manual

## 1. What This Module Does
Association Biplot displays multivariate structure, jointly showing sample configuration and variable directionality.
It is used to interpret how variables align with major ordination axes.

## 2. Left Panel Parameters (Detailed)

### 2.1 Group / Metadata Variable
- Defines the contextual grouping shown in biplot space.
- Strongly affects visual interpretation of separation patterns.

### 2.2 Feature / Taxa Scope
- Selects which taxa/features are included in ordination context.
- Too many features can make vectors and labels unreadable.

### 2.3 Taxonomic Level
- Sets granularity of taxon-level contributors.
- Lower rank increases detail and noise; higher rank improves stability.

### 2.4 Analysis Method
- `dbRDA (capscale)` is distance-based constrained ordination and uses the selected distance metric.
- `CCA` is constrained correspondence analysis and does not use the distance metric input.
- Choose one method and keep it fixed when comparing plots across groups.

### 2.5 Distance Metric
- Applies only when `dbRDA (capscale)` is selected.
- `Bray-Curtis` is commonly used for relative abundance dissimilarity.
- `Aitchison` is based on log-ratio geometry and can be more robust for compositional structure.

### 2.6 Automatic Transformation and Filtering
- Transformation is selected automatically according to the analysis method and distance:
  - `CCA`: TSS transformation (relative abundance).
  - `dbRDA + Bray-Curtis`: TSS transformation followed by Bray-Curtis distance.
  - `dbRDA + Aitchison`: CLR transformation using `log(count + 1)`, centered within each sample, followed by Euclidean distance.
- Prevalence filtering and the maximum-taxa selection are applied to aggregated counts before transformation.
- Taxa with zero or non-finite variance after transformation are removed.
- The selected transformation is reported in the Biplot Status panel and figure legend.
- These method-specific transformations are fixed to prevent unsupported combinations and ensure reproducible analysis settings.

### 2.7 Statistical Testing
- dbRDA reports a PERMANOVA result from `vegan::adonis2`.
- The constrained ordination model is also evaluated using `vegan::anova.cca`.
- Permutation-based tests use 999 permutations, with the random-number generator initialized using seed 42.
- Report the PERMANOVA R2 as an effect-size measure together with its permutation p-value.

### 2.8 Variable Arrow / Loading Display Options
- Controls whether and how variable vectors are drawn.
- Dense vector display can obscure sample pattern.
- Use filtered vectors for presentation-focused views.

### 2.9 Label Controls
- Controls sample/feature label visibility and density.
- Useful for balancing readability vs detail.

### 2.10 Scaling / Display Options
- Adjusts visual emphasis between sample distances and variable vectors.
- Keep scaling consistent when comparing multiple biplots.

### 2.11 Plot Width / Height / Base Font Size
- Increase dimensions for dense labels and multi-group overlays.

## 3. How to Interpret Biplot Elements

### 3.1 Sample points
- Relative distance suggests similarity under displayed ordination structure.
- Apparent clusters should be validated statistically elsewhere.

### 3.2 Arrows (variables)
- Direction: increasing gradient of the variable in ordination space.
- Length: relative contribution/association strength in displayed axes.

### 3.3 Angle between arrows
- Small angle: positive co-direction tendency.
- Opposite direction: inverse tendency.
- Near orthogonal: weak directional relation in shown space.

## 4. Practical Parameter Presets

### Preset A: Exploratory interpretation
- Moderate feature count
- Show key vectors only
- Minimal label clutter

### Preset B: Presentation figure
- Curated variables with clear biological relevance
- Consistent scaling and larger plot dimensions
- Reduced label density for readability

## 5. Common Pitfalls
- Over-interpreting exact distances/angles without considering dimensional reduction limits.
- Displaying too many vectors and losing interpretability.
- Comparing biplots across different scaling/settings without documentation.

## 6. Recommended Reporting Items
- Grouping variable and selected feature scope
- Taxonomic level
- Analysis method, distance metric, and automatically selected transformation
- Prevalence cutoff, maximum taxa, and permutation settings
- Vector/label/scaling settings
- Main directional patterns and biological interpretation
- Any complementary statistical validation used
