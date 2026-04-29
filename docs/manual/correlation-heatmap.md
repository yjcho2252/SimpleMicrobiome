# Correlation Heatmap Manual

## 1. What This Module Does
Correlation Heatmap visualizes association structure among taxa and/or metadata as a matrix.
It is useful for identifying co-varying patterns and candidate relationships.

## 2. Left Panel Parameters (Detailed)

### 2.1 Association Mode
- Typical options: `Abundance-based` and `Prevalence-based` (module-dependent).

#### Abundance-based
- Uses abundance magnitude information.
- Better for gradient-like quantitative relationships.

#### Prevalence-based
- Uses occurrence pattern logic (presence/absence context).
- Better when detection pattern differences are the focus.

### 2.2 Target Variable Selection
- Defines which taxa/metadata are included in rows/columns.
- Broad selection increases coverage but reduces readability.
- Curated selection improves interpretability.

### 2.3 Taxonomic Level
- Controls taxa aggregation rank used in matrix construction.
- Lower rank increases detail and sparsity.
- Higher rank improves stability and visual clarity.

### 2.4 Data Transform Options
- Controls feature scale before association calculation.
- Choice affects dynamic range and sign/magnitude distribution in heatmap.

### 2.5 Correlation / Association Method
- Method determines how pairwise relationships are computed.
- Use method consistent with variable type and distribution assumptions.

### 2.6 Value Scale / Color Mapping
- Defines mapping of association values to colors.
- Symmetric scales are preferred for signed associations.
- Narrow scale highlights subtle differences; wide scale emphasizes strong effects only.

### 2.7 Filtering / Threshold Options (if shown)
- Can hide weak or non-significant associations.
- Useful for decluttering, but may hide borderline patterns.

### 2.8 Plot Width / Height / Base Font Size
- Key for large matrices.
- Increase dimensions for many variables or long labels.

## 3. How to Interpret Heatmap Values

### 3.1 Positive vs negative values
- Positive: variables tend to increase together.
- Negative: inverse pattern tendency.
- Association does not imply direct interaction or causality.

### 3.2 Magnitude
- Larger absolute values indicate stronger modeled association.
- Compare magnitudes only within the same analysis settings.

### 3.3 Block patterns
- Grouped high/low blocks can indicate module-like structure.
- Validate biological plausibility with metadata and external knowledge.

## 4. Practical Parameter Presets

### Preset A: Exploratory scan
- Moderate variable count
- Default association mode/method
- Broad value scale

### Preset B: Report-ready matrix
- Curated variable subset
- Consistent taxonomic level with other analyses
- Symmetric value scale
- Optional significance/threshold filtering for clarity

## 5. Common Pitfalls
- Interpreting color intensity across heatmaps with different scales.
- Overloading too many variables into one matrix.
- Mixing transforms/methods and comparing results as directly equivalent.

## 6. Recommended Reporting Items
- Association mode and method
- Transform and taxonomic level
- Variable selection criteria
- Value scale/threshold settings
- Key positive/negative association clusters
