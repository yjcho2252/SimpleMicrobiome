# Changelog

All notable changes to this project are documented in this file.

## [2026-04-29]

### Added
- Added module guide documentation under `docs/manual/`:
  - `README.md` (guide index)
  - `data-loading.md`
  - `preprocessing.md`
  - `taxa-profiles.md`
  - `alpha-diversity.md`
  - `beta-diversity.md`
  - `ancom-bc2.md`
  - `maaslin2.md`
  - `random-forest.md`
  - `sparcc.md`
  - `spieceasi.md`
  - `correlation-heatmap.md`
  - `association-biplot.md`
  - `citation.md`
- Added browser-rendered guide pages in `docs/manual_html/`.
- Added transparent icon asset:
  - `www/icon_transparent.png`

### Changed
- Top page guide access (`modules/mod_top.R`):
  - Added Guide link on Top page.
  - Updated path to `/manual/README.html`.
- Static resource mapping (`app.R`):
  - Added `addResourcePath("docs", ...)`.
  - Added `addResourcePath("manual", ...)`.
- Manual index links (`docs/manual/README.md`):
  - Updated links to `/manual/*.html` endpoints.
- README branding (`README.md`):
  - Updated top icon to `www/icon_transparent.png`.
- Taxa comparison (`modules/mod_taxa_comparison.R`):
  - Added trend-line controls near Within-Subject Pairing (`Show trend line`, `Trend Line Method`).
  - Added `Scatter plot` mode and auto-switched to scatter when trend line is enabled.
  - Added facet-level trend p-value annotations (`p-value = ...`) during trend-line display.
  - Enabled smoothing interval display for trend lines (`se = TRUE`) and adjusted p-value text placement/margins to reduce overlap.
  - Included `SampleID` in Primary Group choices while defaulting to a non-`SampleID` variable when available.
- Beta diversity (`modules/mod_beta.R`):
  - Included `SampleID` in Primary Variable choices while defaulting to a non-`SampleID` variable when available.
  - Added optional cluster-color overlay on ordination dots with `Cluster` legend categories.
  - Clarified clustering result text: average silhouette width definition and optimal `k` selection rule (maximize average silhouette width).
- Association navigation (`app.R`):
  - Removed the Correlation Plot entry and server wiring after module removal (`mod_scatter`).

### Removed
- Removed previous icon file:
  - `www/SimpleMicrobiome_icon3.png`
- Removed module file:
  - `modules/mod_scatter.R`

### Notes
- Detailed release notes: `docs/releases/2026-04-29.md`

## [2026-04-28]

### Added
- Added a new icon asset for documentation branding updates:
  - `www/SimpleMicrobiome_icon3.png`

### Changed
- README branding (`README.md`):
  - Placed the project icon under the title section.
  - Tuned icon display width for GitHub readability (`width="220"`).

### Notes
- Detailed release notes: `docs/releases/2026-04-28.md`

## [2026-04-27]

### Added
- Replaced the workflow overview SVG with a full end-to-end app workflow diagram:
  - `www/workflow-overview.svg`
- Added explicit subgroup toggles to simplify primary/secondary grouping flows across major analysis modules:
  - `modules/mod_ancom.R`
  - `modules/mod_maaslin2.R`
  - `modules/mod_randomforest.R`
  - `modules/mod_sparcc.R`
  - `modules/mod_SpiecEasi.R`
  - `modules/mod_heatmap.R`
- Added a root-level local dependency installer script:
  - `install.R`

### Changed
- Alpha diversity (`modules/mod_alpha.R`):
  - Updated bar-plot rendering to use facet-aware baseline minima, matching the taxa comparison bar-plot scale behavior.
  - Applied the same baseline logic to regular and `ggpattern` bar plots.
- Heatmap (`modules/mod_heatmap.R`):
  - Renamed `Cell Size` to `Plot Dimensions`.
  - Updated the section icon and moved advanced transform/correlation controls under `Advanced options`.
- Preprocessing (`modules/mod_preprocessing.R`):
  - Added a small visual gap between the `Select All` and `Unselect All` buttons.
- Grouping-dependent modules:
  - Reset subgroup-specific inputs when subgroup selection is disabled.
  - Renumbered visible controls after hiding subgroup-only controls.
  - Improved metadata column resolution for subgroup-aware workflows.
- Association biplot/file loading UI:
  - Refined compact sidebar text/placeholder styling and upload guidance.
- License:
  - Updated the MIT copyright holder to `Kangwon National University`.
- Local setup workflow:
  - Aligned dependency installation logic with Docker build strategy, including Bioconductor and GitHub package installation flow for `SpiecEasi`, `SPRING`, and `NetCoMi`.
  - `install.R`
- Documentation:
  - Split README usage guidance into `Use on Web` and `Run Locally`.
  - Standardized local quick-start sequence to `source("install.R")` then `shiny::runApp("app.R")`.
  - `README.md`

### Notes
- Detailed release notes: `docs/releases/2026-04-27.md`

## [2026-04-26]

### Added
- Added `All` as a default/selectable option for `2. Primary level to include` to align filtering behavior across:
  - `modules/mod_ancom.R`
  - `modules/mod_maaslin2.R`
  - `modules/mod_randomforest.R`
- Added a dedicated `Result` tab in Random Forest for model status/metrics output:
  - `modules/mod_randomforest.R`

### Changed
- Random Forest (`modules/mod_randomforest.R`):
  - Split display responsibilities:
    - Bottom `Random Forest Status` now shows sample-selection context only
    - `Result` tab now shows model status + model metrics
  - Replaced verbose `SampleID preview` with a concise pattern-goal description in status text
  - Added training-target context in status text (classification group list or target type)
  - Increased result panel display area in the `Result` tab
  - Made SHAP bar-plot fonts follow `Base Font Size` for row labels, column labels, and axis text (render + download)
- Heatmap (`modules/mod_heatmap.R`):
  - Added normalized metadata-column resolution fallback
  - Aligned primary-level behavior with `All` option handling
  - Expanded status/legend text for analysis scope and z-score/prevalence context
- ANCOM/MaAsLin2 (`modules/mod_ancom.R`, `modules/mod_maaslin2.R`):
  - Updated primary-level subset logic so `All` behaves as no strict primary-level equality filter
- Taxa comparison (`modules/mod_taxa_comparison.R`):
  - Improved selected-taxa ordering stability and faceted bar-scale behavior
- Preprocessing (`modules/mod_preprocessing.R`):
  - Applied full-row default selection when switching away from preprocessing without manual selection

### Notes
- Detailed release notes: `docs/releases/2026-04-26.md`

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
