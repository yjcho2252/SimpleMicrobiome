# Constrained Ordination Manual

## 1. What This Module Does
Constrained Ordination uses dbRDA or CCA to model multivariate community structure while jointly showing sample configuration and variable directionality.
It is used to interpret how variables align with major ordination axes.

## 2. Left Panel Parameters (Detailed)

**Analysis Settings** includes the collapsible **Advanced Model Settings**, **Taxa Settings**, and **Explanatory Vectors** sections.

### 2.1 Primary and Additional Model Variables
- **Primary explanatory variable** defines the main constrained-model term, sample color, and primary centroid labels when categorical.
- **Additional explanatory variables** extend the constrained model with continuous environmental measurements or categorical metadata.
- **Conditioning variables** are included through `Condition()` to adjust for nuisance structure without drawing them as fitted vectors.
- Continuous model variables can be standardized before fitting, which is recommended when their measurement scales differ.
- The formula preview reports the selected model structure and the number of complete samples before the analysis is run.

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
- `Jaccard` compares presence/absence composition after equal-depth rarefaction.
- `Jensen-Shannon` compares relative-abundance probability profiles.
- `Aitchison` is based on log-ratio geometry and can be more robust for compositional structure.

### 2.6 Automatic Transformation and Filtering
- Transformation is selected automatically according to the analysis method and distance:
  - `CCA`: TSS transformation (relative abundance).
  - `dbRDA + Bray-Curtis`: TSS transformation followed by Bray-Curtis distance.
  - `dbRDA + Jaccard`: rarefaction to the minimum selected-sample depth before taxonomic aggregation, followed by presence/absence transformation and Jaccard distance.
  - `dbRDA + Jensen-Shannon`: TSS transformation followed by Jensen-Shannon distance.
  - `dbRDA + Aitchison`: CLR transformation using `log(count + 1)`, centered within each sample, followed by Euclidean distance.
- Prevalence filtering and the maximum-taxa selection are applied to aggregated counts before transformation.
- Taxa with zero or non-finite variance after transformation are removed.
- The selected transformation is reported in the Model Summary panel and figure legend.
- These method-specific transformations are fixed to prevent unsupported combinations and ensure reproducible analysis settings.

### 2.7 Statistical Testing
- The constrained ordination model is evaluated using `vegan::anova.cca`.
- Marginal permutation tests report each explanatory variable after accounting for the other selected variables.
- `vegan::vif.cca` values are displayed to help identify collinearity among model terms.
- Permutation-based tests use 999 permutations, with the random-number generator initialized using seed 42.

### 2.8 Explanatory Vectors and Taxa Scores
- Taxa arrows are taxa scores from the fitted constrained ordination; individual taxa p-values are not calculated.
- Continuous explanatory variables can be drawn as green dashed explanatory vectors.
- Environmental-vector labels can include marginal permutation p-values, and vectors can be filtered by significance or maximum count.
- Categorical variables are represented by centroid labels instead of directional arrows.
- Dense vector display can obscure sample pattern.
- Use filtered vectors for presentation-focused views.

### 2.9 Label Controls
- Controls sample/feature label visibility and density.
- Sample-to-centroid lines can connect each sample to its primary categorical-group centroid.
- Useful for balancing readability vs detail.

### 2.10 Scaling / Display Options
- Adjusts visual emphasis between sample distances and variable vectors.
- Keep scaling consistent when comparing multiple biplots.

### 2.11 Plot Width / Height / Base Font Size
- The base font defaults to 14 pt and dot size defaults to 5. **Show Dot Outline** adds a dark border to each sample marker.
- Increase dimensions for dense labels and multi-group overlays.
- **Connect Samples Within Subjects** draws a line for subjects with exactly two plotted samples. Choose the metadata column under **Subject ID Variable**.

## 3. How to Interpret Biplot Elements

### 3.1 Sample points
- Relative distance suggests similarity under displayed ordination structure.
- Apparent clusters should be validated statistically elsewhere.

### 3.2 Arrows (variables)
- Direction: increasing gradient of the variable in ordination space.
- Length: relative contribution/association strength in displayed axes.
- Taxa scores and explanatory vectors have different meanings: taxa use fitted feature scores, while explanatory arrows use continuous-variable biplot scores.

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
