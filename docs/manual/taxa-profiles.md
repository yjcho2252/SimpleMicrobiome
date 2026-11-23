# Taxa Profiles Manual

## 1. Scope
Taxa Profiles includes two menus:
- **Taxa Bar Plot**: composition overview
- **Taxa Comparison**: group-wise statistical comparison

Taxa Bar Plot subsamples samples above 2,000 reads to that depth. Samples at or below 2,000 reads retain their counts.

This document explains what each left-panel parameter does and how different settings change interpretation.

---

## 2. Taxa Bar Plot: Left Panel Parameters

### 2.1 Primary Group (Facet)
- **What it controls**: Main metadata grouping used to split panels.
- **If changed**:
  - Group with many levels -> more facets, smaller panel area.
  - Group with 2-3 levels -> clearer side-by-side interpretation.
- **Current selector behavior**:
  - `SampleID` is selectable.
  - `None` displays all samples in a single panel.
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

### 2.7 Taxa Display Mode
- **Options**:
  - `Top N Taxa`
  - `Abundance Threshold`
- **What it controls**: Which taxa are shown explicitly instead of being aggregated into `Others`.
- **Interpretation impact**:
  - `Top N Taxa`: displays the highest-ranked taxa across all selected samples or within each group.
  - `Abundance Threshold`: adaptive view that keeps all taxa above a chosen mean-abundance cutoff.
- **Practical range**:
  - `Top N Taxa`: 10-20 for exploration, 20-30 for denser summaries.
  - `Abundance Threshold`: 0.5-2.0% for most summary plots, depending on sparsity.

### 2.8 Taxa Numbers to Display / Abundance Threshold
- **What it controls**:
  - `Top N Taxa` mode uses the selected count.
  - `Abundance Threshold` mode uses the selected mean-abundance percentage.
- **Interpretation note**:
  - Mean abundance is computed on the relative-abundance scale after sample normalization.
  - Under `Select top taxa from`, `All selected samples` ranks taxa by their mean relative abundance across all selected samples. `Each group` ranks taxa by their mean relative abundance within each group. The common legend includes taxa selected in any group and can contain more than N taxa.
  - Under `Threshold applied to`, `Mean in at least one group` displays a taxon if its mean relative abundance meets the threshold in any group. `Mean across selected samples` uses the overall mean.
  - With a secondary grouping variable, each primary and secondary group combination is evaluated separately for group-based selection.

### 2.9 Taxa Rank to Display
- **Options**:
  - `Full Hierarchy`
  - `Current Rank`
- **Interpretation impact**:
  - `Full Hierarchy`: less ambiguity, longer labels.
  - `Current Rank`: compact legend, possible duplicate names across lineages.

### 2.10 Color Palette
- **What it controls**: Taxa color mapping.
- Additional choices include Okabe–Ito, Pastel, Dynamic, Warm, and Cool.
- **Interpretation impact**:
  - High-contrast palette improves multi-taxa discrimination.
  - Poor contrast can create false visual grouping.
- **Recommendation**:
  - Keep one palette fixed across related figures for consistency.

### 2.11 Plot Width / Height / Base Font Size
- **What it controls**: Figure readability.
- `Fit plot to screen` adjusts the width to the available space. Turn it off to set the width manually.
- Hover over a bar segment to see its taxon and relative abundance. Click a legend item to highlight it across panels.
- Taxa Barplot rarefies to 3,000 reads per sample or the lowest selected read count, whichever is smaller.
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

Enable **Within-Subject Pairing** to compare repeated samples from the same subject across conditions such as visits or treatment time points. The comparison uses the primary grouping variable by default and the secondary grouping variable when one is selected.

- **Two groups** require at least two subjects with observations in both groups after sample filtering. The selected statistical method runs as a paired t-test or paired Wilcoxon signed-rank test.
- **Three or more groups** require at least two subjects with observations in every comparison group for the overall Friedman test. Each pairwise test uses the subjects with observations in both groups being compared.
- When a subject has multiple samples in the same group, the tests use the median abundance for that subject and group.

### 3.4 Subject Identifier

Use a separate metadata column such as `SubjectID`. Give every sample a unique `SampleID`, reuse the exact same `SubjectID` for all samples from one subject, and assign different IDs to different subjects. Rows with missing or empty subject IDs are excluded from paired tests.

| SampleID | SubjectID | Timepoint |
| --- | --- | --- |
| S01_Pre | P01 | Pre |
| S01_Post | P01 | Post |
| S02_Pre | P02 | Pre |
| S02_Post | P02 | Post |

For this example, set **Primary grouping variable** to `Timepoint`, leave **Secondary grouping variable** unset, enable **Within-Subject Pairing**, and select `SubjectID` under **Subject Identifier**.

### 3.5 Show Regression Line / Group variable type / Regression method
- **Show Regression Line**:
  - When enabled, the plot automatically switches to `Scatter plot`.
  - Pairwise p-value bars are turned off while trend-line mode is active.
  - Facet-level trend q-values are shown as `q-value = ...`.
- **Group variable type**:
  - `Categorical`: group levels are encoded as ordered indices for trend fitting.
  - `Continuous (numeric)`: group values are converted to numeric directly.
- **Regression method**:
  - `Spearman` (default): uses a loess smoothing line for visualization.
  - `Pearson`: uses linear-model smoothing.
- **When disabled**:
  - The previously selected plot type is restored.
  - `Show p-value bars` is re-enabled.
  - When pairing is enabled, the repeated-measures logic follows the same x-axis grouping used by the plot, and the overall test switches to `Friedman test` when there are 3 or more x-axis groups.

### 3.6 Taxonomic Level
- Same interpretation principles as Bar plot.
- Lower ranks increase sparsity and multiple-testing burden.

### 3.7 Taxa to Compare
- **What it controls**: Explicit taxon subset for testing/plotting.
- **Interpretation impact**:
  - Broad selection increases discovery chance but also testing burden.
  - Focused selection improves interpretability.

### 3.8 Show only q-value < 0.05
- **Options**: Enabled/Disabled.
- **Interpretation impact**:
  - Enabled: concise significant-only report, but hides near-threshold trends.
  - Disabled: full transparency, easier to assess effect landscape.
- **Calculation mode**:
  - Filtering always uses taxa-level multiple-testing-adjusted q-values (`q < 0.05`) from group-comparison statistics.
  - `Show Regression Line` does not change the taxa filtering rule.
  - When `Within-Subject Pairing` is enabled and group count is 3 or more, overall group comparison uses `Friedman test` (repeated measures) instead of Kruskal-Wallis.
- **Important distinction**:
  - This option filters which taxa are included in the taxa-selection workflow.
  - It is separate from the p-value bar option below.

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
- **Options**: `Violin plot` (default), `Box plot`, `Jitter plot`, `Bar plot`, and `Scatter plot`.
- Dot Size defaults to 4, Base Font Size defaults to 14, and categorical X-axis labels can be rotated to 45° or 90°.
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
- **Show only p-value < 0.05**: hides non-significant brackets after applying the selected pairwise multiple-testing correction.
- **Show as significant marks**: converts numeric labels to symbolic marks.
- **Interpretation note**:
  - Numeric p-values are better for technical reports.
  - Significance marks are better for compact slides.
  - If `Pairwise Multiple Testing Correction` is `None`, the displayed/filtering value is the raw pairwise p-value. Otherwise it is the adjusted pairwise p-value.

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
- **Exclude outliers from plot** applies the 1.5 × IQR rule within each group and panel. The plot is recalculated from retained observations, and statistical tests use all observations. The Comparison Status box describes the rule when it is active.
- **Recommendation**:
  - Multi-taxa comparisons typically need larger width.
  - Increase base font for presentation exports.

---

## 4. Minimal Recommended Presets

### Preset A: Fast exploratory scan
- Taxa Bar Plot:
  - Plot mode = `Sample Level`
  - Taxonomic level = `Genus`
  - Top N = `15`
- Taxa Comparison:
  - Transformation = `CLR`
  - Plot type = `Box plot`
  - Method = `Wilcoxon`
  - Correction = `BH`

### Preset B: Report-ready summary
- Taxa Bar Plot:
  - Plot mode = `Group Mean`
  - Top N = `10-15`
  - Full hierarchy labels for disambiguation
- Taxa Comparison:
  - Keep p-value bars on
  - Use corrected p-values (`BH` or `Holm`)
  - Increase plot size for readable labels
