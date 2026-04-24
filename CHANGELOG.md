# Changelog

All notable changes to this project are documented in this file.

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
