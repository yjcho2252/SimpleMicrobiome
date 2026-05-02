# Taxa Profiles Manual

## 1. Scope
Taxa Profiles includes two menus:
- **Taxa Bar plot**: composition overview
- **Taxa Comparison**: group-wise statistical comparison

This document explains what each left-panel parameter does and how different settings change interpretation.

---

## 2. Taxa Bar plot: Left Panel Parameters

### 2.1 Primary Group (Facet)
- **What it controls**: Main metadata grouping used to split panels.
- **If changed**:
  - Group with many levels -> more facets, smaller panel area.
  - Group with 2-3 levels -> clearer side-by-side interpretation.
- **Current selector behavior**:
  - `SampleID` is selectable.
  - If variables other than `SampleID` exist, a non-`SampleID` variable is used as the default.
- **Recommendation**:
  - Use biologically primary factor (e.g., treatment, disease status).
  - Avoid IDs or high-cardinality columns.

### 2.2 Secondary Group
- **What it controls**: Additional stratification layer (typically color/order context).
- **If changed**:
  - Adds subgroup structure but can make legend crowded.
  - Useful for nested designs (e.g., Treatment within Timepoint).
- **Current selector behavior**:
  - `SampleID` is selectable unless already used as primary grouping.
  - `None` remains the default selected option.
- **Recommendation**:
  - Keep levels limited.
  - Use only when subgroup interpretation is needed.

### 2.3 Plot Display Mode
- **Options**:
  - `Sample Level`: each sample shown individually.
  - `Group Mean`: group-averaged composition.
- **Interpretation impact**:
  - `Sample Level` preserves within-group variability.
  - `Group Mean` suppresses variance and shows trend only.
- **When to use**:
  - Early QC/exploration: `Sample Level`
  - Summary figure for report: `Group Mean`

### 2.4 Primary Sort Criterion (within Group)
- **Options**:
  - `Top Taxa Abundance`
  - `Sample Name`
  - `Metadata Variable`
- **Interpretation impact**:
  - Sorting changes visual pattern readability, not underlying data.
  - `Top Taxa Abundance` emphasizes dominant-taxa gradient.
  - `Sample Name` is best for reproducible fixed ordering.
  - `Metadata Variable` aligns bars with a specific covariate trend.

### 2.5 Secondary Sort (after Metadata)
- **Options**:
  - `None`
  - `Top Taxa Abundance`
  - `Sample Name`
- **Use case**:
  - Resolves tie/within-bin order after metadata-based primary sorting.

### 2.6 Taxonomic Level
- **Common choices**: Phylum/Class/Order/Family/Genus/Species (module-dependent availability).
- **Interpretation impact**:
  - Higher rank (e.g., Phylum): stable, broad signal, less detail.
  - Lower rank (e.g., Genus/Species): detailed but sparse/noisier.
- **Recommendation**:
  - Start at `Genus` for balance.
  - Move to higher level if plot is too fragmented.

### 2.7 Taxa Numbers to Display (Top N)
- **What it controls**: Number of most abundant taxa shown explicitly.
- **Interpretation impact**:
  - Low N: clean chart but may hide meaningful minor taxa.
  - High N: more complete but harder to read.
- **Practical range**:
  - Exploratory view: 10-20
  - Supplementary detail: 20-30

### 2.8 Taxa Rank to Display
- **Options**:
  - `Full Hierarchy`
  - `Current Rank`
- **Interpretation impact**:
  - `Full Hierarchy`: less ambiguity, longer labels.
  - `Current Rank`: compact legend, possible duplicate names across lineages.

### 2.9 Color Palette
- **What it controls**: Taxa color mapping.
- **Interpretation impact**:
  - High-contrast palette improves multi-taxa discrimination.
  - Poor contrast can create false visual grouping.
- **Recommendation**:
  - Keep one palette fixed across related figures for consistency.

### 2.10 Plot Width / Height / Base Font Size
- **What it controls**: Figure readability.
- **Guideline**:
  - Many taxa or facets -> increase width/height.
  - Long taxon names -> increase width and/or reduce displayed taxa.

---

## 3. Taxa Comparison: Left Panel Parameters

### 3.1 Primary grouping variable
- **What it controls**: Main comparison axis for statistical tests.
- **Interpretation impact**:
  - Defines all p-value contexts.
  - Wrong grouping variable = invalid biological interpretation.

### 3.2 Secondary grouping variable
- **What it controls**: Additional subgroup context.
- **Interpretation impact**:
  - Helps stratified interpretation.
  - Can reduce per-cell sample size and lower power.

### 3.3 Within-Subject Pairing
- **Options**: Enabled/Disabled.
- **Interpretation impact**:
  - Enabled: paired logic (within-subject change focus).
  - Disabled: independent-group comparison.
- **Critical note**:
  - Enable only when repeated measures truly match by subject ID.

### 3.4 Subject Identifier
- **Used when pairing is enabled**.
- **Requirement**:
  - Must uniquely map repeated samples to the same subject.
- **Risk**:
  - Incorrect ID column causes invalid paired statistics.

### 3.5 Show Regression Line / Group variable type / Regression method
- **Show Regression Line**:
  - When enabled, the plot automatically switches to `Scatter plot`.
  - Pairwise p-value bars are turned off while trend-line mode is active.
  - Facet-level trend p-values are shown as `p-value = ...`.
- **Group variable type**:
  - `Auto`: if group values are numeric for enough samples, uses continuous mode; otherwise categorical encoding.
  - `Categorical`: group levels are encoded as ordered indices for trend fitting.
  - `Continuous (numeric)`: group values are converted to numeric directly.
- **Regression method**:
  - `Spearman` (default): uses a loess smoothing line for visualization.
  - `Pearson`: uses linear-model smoothing.
- **When disabled**:
  - Plot type returns to `Bar plot`.
  - `Show p-value bars` is re-enabled.

### 3.6 Taxonomic Level
- Same interpretation principles as Bar plot.
- Lower ranks increase sparsity and multiple-testing burden.

### 3.7 Taxa to Compare
- **What it controls**: Explicit taxon subset for testing/plotting.
- **Interpretation impact**:
  - Broad selection increases discovery chance but also testing burden.
  - Focused selection improves interpretability.

### 3.8 Show only taxa with p-value < 0.05
- **Options**: Enabled/Disabled.
- **Interpretation impact**:
  - Enabled: concise significant-only report, but hides near-threshold trends.
  - Disabled: full transparency, easier to assess effect landscape.
- **Calculation mode**:
  - When `Show Regression Line` is disabled, filtering uses group-comparison p-values (two-group test or multi-group Kruskal-Wallis).
  - When `Show Regression Line` is enabled, filtering is recalculated from the selected `Regression method` (`Spearman` or `Pearson`) and `Group variable type`.
  - When `Within-Subject Pairing` is enabled and group count is 3 or more, overall group comparison uses `Friedman test` (repeated measures) instead of Kruskal-Wallis.

### 3.9 Transformation
- **Options**:
  - `CLR Abundance`
  - `Relative Abundance (%)`
  - `Log TSS (log10(% + 1))`
- **Interpretation impact**:
  - `CLR`: compositional-aware, centered log-ratio scale.
  - `Relative`: intuitive percent scale, but compositional constraints remain.
  - `Log TSS`: compresses extreme values and improves visual stability.
- **Recommendation**:
  - Primary inference: `CLR`
  - Communication figure: `Relative` or `Log TSS` as supplementary.

### 3.10 Plot Type
- **Options**: `Box plot`, `Bar plot`, `Scatter plot`.
- **Interpretation impact**:
  - `Box plot`: shows distribution/spread/outliers.
  - `Bar plot`: emphasizes mean-level differences.
  - `Scatter plot`: emphasizes sample-level point structure; useful with trend lines.
- **Recommendation**:
  - Use box plot for statistical interpretation.
  - Use bar plot for summary presentation.

### 3.11 Color Palette / Use ggpattern
- **Color Palette**: group visual distinction.
- **Use ggpattern (bar plot)**:
  - Adds pattern cues for accessibility/print contexts.
  - Slightly heavier rendering.

### 3.12 Show p-value bars / significance display
- **Show p-value bars**: toggles bracket annotations.
- **Show only p-value < 0.05**: hides non-significant brackets.
- **Show as significant marks**: converts numeric labels to symbolic marks.
- **Interpretation note**:
  - Numeric p-values are better for technical reports.
  - Significance marks are better for compact slides.

### 3.13 Statistical Method
- **Options**: `Wilcoxon`, `T-test`.
- **Selection logic**:
  - `Wilcoxon`: robust for non-normal/heavy-tailed data.
  - `T-test`: suitable when approximate normality assumptions are acceptable.

### 3.14 Pairwise Multiple Testing Correction
- **Options**: `None`, `Holm`, `BH (FDR)`, `Bonferroni`.
- **Interpretation impact**:
  - `None`: highest sensitivity, highest false-positive risk.
  - `BH`: balanced default for microbiome screening.
  - `Bonferroni`: strict, low false positives, high false negatives.

### 3.15 Manual facet rows/cols (ncol/nrow)
- **What it controls**: explicit panel grid.
- **Use case**:
  - Standardize layout across figures for publication consistency.

### 3.16 Plot Width / Height / Base Font Size
- **Recommendation**:
  - Multi-taxa comparisons typically need larger width.
  - Increase base font for presentation exports.

---

## 4. Minimal Recommended Presets

### Preset A: Fast exploratory scan
- Taxa Bar plot:
  - Plot mode = `Sample Level`
  - Taxonomic level = `Genus`
  - Top N = `15`
- Taxa Comparison:
  - Transformation = `CLR`
  - Plot type = `Box plot`
  - Method = `Wilcoxon`
  - Correction = `BH`

### Preset B: Report-ready summary
- Taxa Bar plot:
  - Plot mode = `Group Mean`
  - Top N = `10-15`
  - Full hierarchy labels for disambiguation
- Taxa Comparison:
  - Keep p-value bars on
  - Use corrected p-values (`BH` or `Holm`)
  - Increase plot size for readable labels
