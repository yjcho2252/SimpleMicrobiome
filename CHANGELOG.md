# Changelog

All notable changes to this project are documented in this file.

## [2026-04-25]

### Added
- Added `Jaccard` distance option in beta diversity and wired it through ordination/permanova displays:
  - `modules/mod_beta.R`
- Added prevalence-focused association mode in heatmap:
  - `Association mode: Abundance-based / Prevalence-based (binary)`
  - Binary-aware testing and reporting updates
  - `modules/mod_heatmap.R`

### Changed
- Beta diversity (`modules/mod_beta.R`):
  - Clarified Aitchison labeling as `Aitchison (CLR log(x+1) pseudocount)`
  - EnvFit label cleanup:
    - Removed `Taxa::` prefix from taxa labels
    - Removed group-name prefixes from factor arrow labels
    - Applied bold font to EnvFit labels
  - Added origin marker (dot at `(0,0)`) for EnvFit arrow vectors
  - Ensured Jaccard uses binary presence/absence transformation
- Heatmap (`modules/mod_heatmap.R`):
  - Added prevalence mode behavior for taxa/taxa and taxa/group association workflows
  - For prevalence + categorical groups, heatmap value changed to delta prevalence (`in-group - out-group`)
  - Disabled controls not applicable in prevalence mode (`Data transform`, `Correlation method`, `Value scale`)
  - Updated legend title text and styling for clearer mode-specific interpretation
- Taxa comparison (`modules/mod_taxa_comparison.R`):
  - Aligned taxa filtering p-value logic (`show_p_lt_0_05_only`) with selected statistical method and paired mode
  - Updated CLR y-scale handling to avoid clipping negative values
  - Unified box/bar plot y-axis scale basis while allowing bar rendering from axis baseline
  - Refined y-axis break/label formatting and auto-width scaling factor (`*10 -> *20`)
  - Restored visible boxplot outlier points
- Bar plot (`modules/mod_barplot.R`):
  - Added zero-sum-safe relative abundance normalization guards
- Alpha diversity (`modules/mod_alpha.R`):
  - Updated legend significance text to reflect 2-group vs 3+-group comparisons and selected correction method

### Fixed
- Fixed taxa selection reset bug in taxa comparison when changing secondary grouping with multiple taxa preselected:
  - `modules/mod_taxa_comparison.R`

### Notes
- Detailed release notes: `docs/releases/2026-04-25.md`

## [2026-04-24]

### Added
- Introduced a two-stage grouping workflow (Primary + Secondary) in:
  - `modules/mod_ancom.R`
  - `modules/mod_maaslin2.R`
  - `modules/mod_randomforest.R`
  - `modules/mod_sparcc.R`
  - `modules/mod_spieceasi.R`
  - `modules/mod_heatmap.R`

### Changed
- Updated ANCOM/MaAsLin2 group selector wording for clarity:
  - `Secondary group levels to include` -> `Comparison groups (levels)`
  - Placeholder updated to `Select two or more groups to compare`
- Reduced placeholder/wait-message font size across major modules (`graphics::text(..., cex = 0.85)`).

### Fixed
- Fixed table output overlapping legend/status UI in:
  - `modules/mod_ancom.R`
  - `modules/mod_maaslin2.R`

### Notes
- Detailed release notes: `docs/releases/2026-04-24.md`
