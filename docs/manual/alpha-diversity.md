# Alpha Diversity Manual

## 1. What This Module Does
Alpha Diversity quantifies **within-sample diversity** and compares those values across groups.
This module supports four indices:
- `Observed`
- `Chao1`
- `Shannon`
- `Simpson (1-D)`

It also provides group-wise statistical comparisons with optional multiple-testing correction.

## 2. Left Panel Parameters (Detailed)

### 2.1 Primary Group
- **What it controls**: Main metadata variable used for comparison and x-axis grouping.
- **Effect of setting**:
  - Few levels (2-4): easier pairwise interpretation.
  - Many levels: more pairwise tests and denser p-value annotations.
- **Recommendation**:
  - Use the biologically primary variable (e.g., treatment, disease class).
  - Avoid ID-like columns.

### 2.2 Secondary Group
- **What it controls**: Additional stratification/facet context.
- **Effect of setting**:
  - Helps separate subgroup patterns.
  - Reduces effective sample count per panel; significance can drop.
- **Recommendation**:
  - Enable only when subgroup interpretation is required.

### 2.3 Index Selection (`Observed`, `Chao1`, `Shannon`, `Simpson`)
- You can select one or multiple indices at once.
- **Effect of selecting multiple indices**:
  - Creates multi-index panel output.
  - Increases interpretation scope but also figure complexity.

#### Index meaning and when to use
- **Observed**
  - Definition: number of observed taxa/features (raw richness count).
  - Sensitive to sequencing depth and rare taxa detection.
  - Use when you want a direct richness count summary.
- **Chao1**
  - Definition: richness estimator that accounts for unseen taxa using rare-feature information.
  - Typically larger than Observed when many rare taxa exist.
  - Use when richness under-detection is a concern.
- **Shannon**
  - Definition: diversity metric combining richness and evenness.
  - Increases when community has both many taxa and balanced proportions.
  - Use as a balanced “overall diversity” index.
- **Simpson (1-D)**
  - Definition: dominance/evenness-sensitive diversity (higher = more even, less dominated).
  - More influenced by abundant taxa than very rare taxa.
  - Use when dominance structure is important.

### 2.4 Plot Type
- **Options**: `Box plot`, `Bar plot`
- **Effect of setting**:
  - `Box plot`: shows distribution, spread, and outliers per group.
  - `Bar plot`: emphasizes group summary (mean-like view), less detail on dispersion.
- **Recommendation**:
  - Statistical interpretation: `Box plot`
  - Presentation summary: `Bar plot`

### 2.5 Color Palette
- **What it controls**: Group color mapping.
- **Effect of setting**:
  - High-contrast palettes improve multi-group readability.
- **Recommendation**:
  - Keep a consistent palette across related figures.

### 2.6 Use ggpattern (Bar plot)
- **What it controls**: Pattern fill in bars (bar plot mode only).
- **Effect of setting**:
  - Improves grayscale/print accessibility.
  - Can increase rendering complexity.

### 2.7 Show p-value bars
- **What it controls**: Whether pairwise comparison bars are drawn.
- **Effect of setting**:
  - Enabled: direct statistical overlay on plot.
  - Disabled: cleaner figure, no inline significance context.

### 2.8 Show only p-value < 0.05
- Active when p-value bars are enabled.
- **Effect of setting**:
  - Enabled: only significant comparisons displayed.
  - Disabled: all pairwise comparisons shown (can clutter plot).

### 2.9 Show as significant marks
- Active when p-value bars are enabled.
- **Effect of setting**:
  - Enabled: symbolic marks (compact display).
  - Disabled: numeric p-values (more explicit reporting).

### 2.10 Statistical Method
- **Options**: `Wilcoxon`, `T-test`
- **When to use**:
  - `Wilcoxon`: safer for non-normal/skewed distributions.
  - `T-test`: use when normality/variance assumptions are acceptable.
- **Recommendation**:
  - Default to `Wilcoxon` for robust microbiome exploratory comparison.

### 2.11 Pairwise Multiple Testing Correction
- **Options**: `None`, `Holm`, `BH`, `Bonferroni`
- **Effect of setting**:
  - `None`: most sensitive, highest false-positive risk.
  - `BH`: balanced FDR control, good default for screening.
  - `Holm`: stricter than BH, less strict than Bonferroni.
  - `Bonferroni`: very strict, lowest false positives, higher false negatives.
- **Recommendation**:
  - Use `BH` for routine reporting.

### 2.12 Plot Width / Plot Height / Base Font Size
- **What it controls**: readability of multi-group/multi-index output.
- **Practical guidance**:
  - Many groups or long labels -> increase width.
  - Many selected indices -> increase height.
  - Presentation figure -> increase base font.

## 3. How Parameter Choices Change Interpretation

### Scenario A: Richness-focused question
- Use `Observed` + `Chao1`.
- Prefer `Box plot` to inspect group spread.
- Use `Wilcoxon` + `BH` for robust pairwise inference.

### Scenario B: Community balance question
- Use `Shannon` + `Simpson`.
- If one group is highly dominated by few taxa, Simpson often highlights it strongly.

### Scenario C: Clean publication figure
- Use limited index set (1-2 indices).
- Turn on p-value bars, keep corrected method (`BH` or `Holm`).
- Increase width/height and base font for readability.

## 4. Common Pitfalls
- Comparing p-values without checking group sample sizes.
- Treating different indices as interchangeable evidence.
- Reporting uncorrected pairwise significance when many comparisons exist.
- Over-interpreting bar plots without distribution context.

## 5. Minimal Recommended Preset
- Indices: `Chao1`, `Shannon`
- Plot type: `Box plot`
- Statistical method: `Wilcoxon`
- Multiple testing correction: `BH`
- p-value bars: On
- Show only p < 0.05: On
