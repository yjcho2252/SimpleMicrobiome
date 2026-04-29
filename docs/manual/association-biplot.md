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

### 2.4 Variable Arrow / Loading Display Options
- Controls whether and how variable vectors are drawn.
- Dense vector display can obscure sample pattern.
- Use filtered vectors for presentation-focused views.

### 2.5 Label Controls
- Controls sample/feature label visibility and density.
- Useful for balancing readability vs detail.

### 2.6 Scaling / Display Options
- Adjusts visual emphasis between sample distances and variable vectors.
- Keep scaling consistent when comparing multiple biplots.

### 2.7 Plot Width / Height / Base Font Size
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
- Vector/label/scaling settings
- Main directional patterns and biological interpretation
- Any complementary statistical validation used
