# Taxa Bubble Plot Manual

Use Taxa Bubble Plot to compare how often taxa are detected and how abundant they are across groups.

1. Choose a grouping variable. Choose **ASV** or an available taxonomic rank. The initial selection is **Genus** when available and **ASV** otherwise.
2. Set **Top Taxa by Mean Relative Abundance** to control how many taxa are included.
3. Under **Display Options**, choose what bubble size and color represent. Read both legends when comparing groups. Prevalence is the percentage of samples in a group in which the taxon was detected. Mean relative abundance summarizes the taxon's share of reads across samples in the group.

Use **Plot (PNG)** to save the figure and **Summary (TSV)** to export the underlying group summaries.

By default, the plot includes the top 20 taxa. Bubble size represents prevalence, bubble color uses **Auto**, and the palette is **Plasma**. X-axis labels are rotated 45°. Plot sizing defaults to **Auto**. Manual dimensions start at 500 × 700 px, and the base font size is 12. Bubble outlines are off, vertical grid lines are visible, and both the fixed and maximum bubble size settings start at 7.
