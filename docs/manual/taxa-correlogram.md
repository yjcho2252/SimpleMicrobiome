# Taxa Correlogram Manual

Use Taxa Correlogram to visualize pairwise correlations between taxa across samples.

1. Analyze all samples or use **Primary Grouping Variable** and **Group Level to Include** to focus on a subset.
2. Choose **ASV** or an available taxonomic rank. The initial selection is **Genus** when available and **ASV** otherwise. Use **Prevalence Filter (%)** and **Max Taxa** to limit the matrix. Eligible taxa are ranked by mean relative abundance.
3. Choose abundance-based or presence/absence association. For abundance, select CLR or TSS and Pearson or Spearman correlation. The presence/absence mode uses phi correlation.
4. Read the circle color for the sign and the circle size for the magnitude of each correlation. Use **Mask Non-Significant Pairs (FDR)** and **FDR Cutoff (q)** to hide pairs above the adjusted significance threshold.

Use **Plot (PNG)** to save the figure and **Matrix (TSV)** to export the pairwise results. Correlation describes an association across samples and does not establish a direct interaction.

The default prevalence filter is 10%, and **Max Taxa** is 20. FDR masking and significance marks are off by default, and the q cutoff is 0.05. **Cluster Taxa** is on, the default palette is blue–red, and the maximum circle size is 8. Plot sizing defaults to **Auto**. Manual dimensions start at 850 × 750 px, and the base font size is 11.
