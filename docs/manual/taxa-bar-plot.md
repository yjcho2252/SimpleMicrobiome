# Taxa Bar Plot Manual

Use Taxa Bar Plot to compare the relative abundance and composition of taxa across samples or groups.

Taxa Bar Plot subsamples samples above 2,000 reads to that depth. Samples at or below 2,000 reads retain their counts.

## Left Panel Parameters

### Primary Group (Facet)
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

### Secondary Group
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

### Plot Display Mode
- **Options**:
  - `Sample Level`: each sample shown individually.
  - `Group Mean`: group-averaged composition.
- **Interpretation impact**:
  - `Sample Level` preserves within-group variability.
  - `Group Mean` suppresses variance and shows trend only.
- **When to use**:
  - Early QC/exploration: `Sample Level`
  - Summary figure for report: `Group Mean`

### Primary Sort Criterion (within Group)
- **Options**:
  - `Top Taxa Abundance`
  - `Sample Name`
  - `Metadata Variable`
- **Interpretation impact**:
  - Sorting changes visual pattern readability, not underlying data.
  - `Top Taxa Abundance` emphasizes dominant-taxa gradient.
  - `Sample Name` is best for reproducible fixed ordering.
  - `Metadata Variable` aligns bars with a specific covariate trend.

### Secondary Sort (after Metadata)
- **Options**:
  - `None`
  - `Top Taxa Abundance`
  - `Sample Name`
- **Use case**:
  - Resolves tie/within-bin order after metadata-based primary sorting.

### Taxonomic Level
- **Common choices**: Phylum/Class/Order/Family/Genus/Species (module-dependent availability).
- **Interpretation impact**:
  - Higher rank (e.g., Phylum): stable, broad signal, less detail.
  - Lower rank (e.g., Genus/Species): detailed but sparse/noisier.
- **Recommendation**:
  - Start at `Genus` for balance.
  - Move to higher level if plot is too fragmented.

### Taxa Display Mode
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

### Taxa Numbers to Display / Abundance Threshold
- **What it controls**:
  - `Top N Taxa` mode uses the selected count.
  - `Abundance Threshold` mode uses the selected mean-abundance percentage.
- **Interpretation note**:
  - Mean abundance is computed on the relative-abundance scale after sample normalization.
  - Under `Select top taxa from`, `All selected samples` ranks taxa by their mean relative abundance across all selected samples. `Each group` ranks taxa by their mean relative abundance within each group. The common legend includes taxa selected in any group and can contain more than N taxa.
  - Under `Threshold applied to`, `Mean in at least one group` displays a taxon if its mean relative abundance meets the threshold in any group. `Mean across selected samples` uses the overall mean.
  - With a secondary grouping variable, each primary and secondary group combination is evaluated separately for group-based selection.

### Taxa Rank to Display
- **Options**:
  - `Full Hierarchy`
  - `Current Rank`
- **Interpretation impact**:
  - `Full Hierarchy`: less ambiguity, longer labels.
  - `Current Rank`: compact legend, possible duplicate names across lineages.

### Color Palette
- **What it controls**: Taxa color mapping.
- Additional choices include Okabe–Ito, Pastel, Dynamic, Warm, and Cool.
- **Interpretation impact**:
  - High-contrast palette improves multi-taxa discrimination.
  - Poor contrast can create false visual grouping.
- **Recommendation**:
  - Keep one palette fixed across related figures for consistency.

### Plot Width / Height / Base Font Size
- **What it controls**: Figure readability.
- `Plot Size` defaults to `Auto`, which fits the width to the available space. Choose `Manual` to set the width and height. The base font defaults to 12 pt.
- Hover over a bar segment to see its taxon and relative abundance. Click a legend item to highlight it across panels.
- **Guideline**:
  - Many taxa or facets -> increase width/height.
  - Long taxon names -> increase width and/or reduce displayed taxa.

## Practical Presets

- **Exploratory scan**. Choose **Sample Level**, **Genus**, and a top-15 taxa limit.
- **Report summary**. Choose **Group Mean**, show the top 10–15 taxa, and use full hierarchy labels to disambiguate names.
