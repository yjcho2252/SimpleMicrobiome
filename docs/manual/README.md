# SimpleMicrobiome Guide

## v1.2.1 Highlights

- Alpha and Beta Diversity exports are in the left settings panel. Use **Analysis Method** in Beta Diversity to switch between PCoA and NMDS.
- The manual index links directly to Taxa Bubble Plot, Taxa Correlogram, and UpSet Plot.

## v1.2.0 Highlights

- Taxa Profiles includes Taxa Bubble Plot for comparing abundance and prevalence across groups.
- Association includes Taxa Correlogram for pairwise correlations and Upset Plot for group overlap, alongside Association Heatmap.
- On larger screens, scroll settings and results independently. Figure legends and status cards appear in the results panel.
- Data Loading accepts QZA feature and taxonomy files through the configured converter and Excel files for metadata.

## Module Guides
- [Data Loading](/manual/data-loading.html)
- [Preprocessing](/manual/preprocessing.html)
- [Taxa Profiles](/manual/taxa-profiles.html)
  - [Taxa Bubble Plot](/manual/taxa-profiles.html#taxa-bubble-plot)
- [Alpha Diversity](/manual/alpha-diversity.html)
- [Beta Diversity](/manual/beta-diversity.html)
- [Constrained Ordination](/manual/association-biplot.html)
- [ANCOM-BC2](/manual/ancom-bc2.html)
- [MaAsLin2](/manual/maaslin2.html)
- [MaAsLin3](/manual/maaslin3.html)
- [Random Forest](/manual/random-forest.html)
- [SparCC](/manual/sparcc.html)
- [SPIEC-EASI](/manual/spieceasi.html)
- [Association Heatmap](/manual/correlation-heatmap.html)
- [Taxa Correlogram](/manual/correlation-heatmap.html#taxa-correlogram)
- [UpSet Plot](/manual/correlation-heatmap.html#upset-plot)

## Quick Goal-Based Path
- **Start-to-finish first analysis**:
  - Data Loading -> Preprocessing -> Taxa Profiles -> Alpha/Beta Diversity
- **Differential abundance focus**:
  - ANCOM-BC2, MaAsLin2, or MaAsLin3
- **Prediction focus**:
  - Random Forest
- **Network focus**:
  - SparCC + SPIEC-EASI
- **Association visualization focus**:
  - Association Heatmap for taxa and metadata, Taxa Correlogram for taxon pairs, or Upset Plot for taxa shared across groups.
- **Constrained community analysis**:
  - Constrained Ordination under Diversity

## Notes
- Settings and default values may evolve with UI updates.
- When reporting results, include key parameter settings from each module.

## Cite SimpleMicrobiome

If you use SimpleMicrobiome in your research, please cite the following article.

Na, S.-I., Kim, J., Kim, S.-Y., Park, J., & Cho, Y.-J. (2026). SimpleMicrobiome: An integrated web-based platform for streamlined microbiome data analysis and visualization. *Journal of Microbiology, 64*(9), e2606011. [View article (DOI: 10.71150/jm.2606011)](https://www.jmicrobiol.or.kr/journal/view.php?number=3102)

The app's first screen and Citation tab both display this citation. The Citation tab also lists references for individual analysis methods under **Method References**.
