# Taxa Comparison Manual

Use Taxa Comparison to compare selected taxa across groups or explore trends along a continuous metadata variable.

## Left Panel Parameters

### Primary grouping variable
- **What it controls**: Main comparison axis for statistical tests.
- **Interpretation impact**:
  - Defines all p-value contexts.
  - Wrong grouping variable = invalid biological interpretation.

### Secondary grouping variable
- **What it controls**: Additional subgroup context.
- **Interpretation impact**:
  - Helps stratified interpretation.
  - Can reduce per-cell sample size and lower power.

### Within-Subject Pairing

Enable **Within-Subject Pairing** to compare repeated samples from the same subject across conditions such as visits or treatment time points. The comparison uses the primary grouping variable by default and the secondary grouping variable when one is selected.

- **Two groups** require at least two subjects with observations in both groups after sample filtering. The selected statistical method runs as a paired t-test or paired Wilcoxon signed-rank test.
- **Three or more groups** require at least two subjects with observations in every comparison group for the overall Friedman test. Each pairwise test uses the subjects with observations in both groups being compared.
- When a subject has multiple samples in the same group, the tests use the median abundance for that subject and group.

### Subject Identifier

Use a separate metadata column such as `SubjectID`. Give every sample a unique `SampleID`, reuse the exact same `SubjectID` for all samples from one subject, and assign different IDs to different subjects. Rows with missing or empty subject IDs are excluded from paired tests.

| SampleID | SubjectID | Timepoint |
| --- | --- | --- |
| S01_Pre | P01 | Pre |
| S01_Post | P01 | Post |
| S02_Pre | P02 | Pre |
| S02_Post | P02 | Post |

For this example, set **Primary grouping variable** to `Timepoint`, leave **Secondary grouping variable** unset, enable **Within-Subject Pairing**, and select `SubjectID` under **Subject Identifier**.

### Show Regression Line / Group variable type / Regression method
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

### Taxonomic Level
- Same interpretation principles as Bar plot.
- Lower ranks increase sparsity and multiple-testing burden.

### Taxa to Compare
- **What it controls**: Explicit taxon subset for testing/plotting.
- **Interpretation impact**:
  - Broad selection increases discovery chance but also testing burden.
  - Focused selection improves interpretability.

### Show only q-value < 0.05
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

### Transformation
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

### Plot Type
- **Options**: `Violin plot` (default), `Box plot`, `Jitter plot`, `Bar plot`, and `Scatter plot`.
- Dot Size defaults to 4, Base Font Size defaults to 14, and categorical X-axis labels can be rotated to 45° or 90°.
- **Interpretation impact**:
  - `Box plot`: shows distribution/spread/outliers.
  - `Bar plot`: emphasizes mean-level differences.
  - `Scatter plot`: emphasizes sample-level point structure; useful with trend lines.
- **Recommendation**:
  - Use box plot for statistical interpretation.
  - Use bar plot for summary presentation.

### Color Palette / Use ggpattern
- **Color Palette**: group visual distinction.
- **Use ggpattern (bar plot)**:
  - Adds pattern cues for accessibility/print contexts.
  - Slightly heavier rendering.

### Show p-value bars / significance display
- **Show p-value bars**: toggles bracket annotations.
- **Show only p-value < 0.05**: hides non-significant brackets after applying the selected pairwise multiple-testing correction.
- **Show as significant marks**: converts numeric labels to symbolic marks.
- **Interpretation note**:
  - Numeric p-values are better for technical reports.
  - Significance marks are better for compact slides.
  - If `Pairwise Multiple Testing Correction` is `None`, the displayed/filtering value is the raw pairwise p-value. Otherwise it is the adjusted pairwise p-value.

### Statistical Method
- **Options**: `Wilcoxon`, `T-test`.
- **Selection logic**:
  - `Wilcoxon`: robust for non-normal/heavy-tailed data.
  - `T-test`: suitable when approximate normality assumptions are acceptable.

### Pairwise Multiple Testing Correction
- **Options**: `None`, `Holm`, `BH (FDR)`, `Bonferroni`.
- **Interpretation impact**:
  - `None`: highest sensitivity, highest false-positive risk.
  - `BH`: balanced default for microbiome screening.
  - `Bonferroni`: strict, low false positives, high false negatives.

### Manual facet rows/cols (ncol/nrow)
- **What it controls**: explicit panel grid.
- **Use case**:
  - Standardize layout across figures for publication consistency.

### Plot Width / Height / Base Font Size
- **Exclude outliers from plot** applies the 1.5 × IQR rule within each group and panel. The plot is recalculated from retained observations, and statistical tests use all observations. The Comparison Status box describes the rule when it is active.
- **Recommendation**:
  - Multi-taxa comparisons typically need larger width.
  - Increase base font for presentation exports.

## Practical Presets

- **Exploratory scan**. Use CLR, a box plot, Wilcoxon tests, and BH correction where appropriate.
- **Report summary**. Keep p-value bars visible, use BH or Holm correction, and increase plot size for readable labels.
