# Core & Indicator Taxa

Use **Core Taxa** to see which taxa meet the detection and prevalence thresholds in each group. Use **Indicator Taxa (IndVal)** to test which taxa characterize a group or a group combination.

## Core Taxa UpSet

The results table and TSV export include group prevalence and mean relative abundance (%), averaged across all samples in the group used for the analysis.

1. Choose a **Grouping Variable** and select 2 to 12 groups under **Groups to Include**.
2. Choose **ASV** or an available taxonomic rank. The initial selection is **Genus** when available and **ASV** otherwise. Set **Prevalence Filter (%)**. A taxon belongs to a group when it is detected in at least one sample and meets the threshold within that group.
3. Read the filled dots to identify the groups in each intersection. The top bar counts taxa belonging to exactly those groups among the selected groups. The side bars show the total number of included taxa in each group.
4. Use **Groups per Pattern** and **Groups in Pattern** to filter intersections. Select **Top 15**, **Top 25**, or **All** under **Patterns to Show**. The Results Table includes taxa from all matching intersections, including those outside the plot limit. Select a plot column to view only the taxa in that intersection.

Use **Core PNG** to save the figure and **Core TSV** to export the taxa in the Results Table, including any active pattern selection.

## Indicator Taxa by Group Pattern

Select **Indicator Taxa (IndVal)** under **Analysis Mode**. Set **Maximum Groups per Pattern** to 2 or 3 to evaluate combinations alongside individual groups. Select at least one more group than the maximum pattern size so each candidate pattern has a comparison group. Use the same display filters as in Core Taxa and click a plot column to view its taxa in Results Table. Click the same column again to clear the selection. Use **IndVal PNG** to save the figure and **IndVal TSV** to export the current Results Table.

After data are loaded, the app keeps a valid grouping selection or selects the first available metadata column other than sample ID. It selects all available group levels when there are at most 12. When there are more than 12, choose the levels to include. Core Taxa defaults to 50% minimum group prevalence and 0.01% minimum relative abundance. The plot shows up to 25 patterns by default. Plot sizing defaults to **Auto**. Manual dimensions start at 900 × 550 px, and the base font size is 11.
