# SimpleMicrobiome Guide

## v1.1.4 Highlights
- Diversity now separates Alpha Diversity, Beta Diversity, and Constrained Ordination.
- Constrained Ordination supports expanded model variables, dbRDA distance choices, and optional taxa/environmental displays.
- Association Heatmap naming and 300 dpi PNG exports are used consistently across the app.

## Module Guides
- [Data Loading](/manual/data-loading.html)
- [Preprocessing](/manual/preprocessing.html)
- [Taxa Profiles](/manual/taxa-profiles.html)
- [Alpha Diversity](/manual/alpha-diversity.html)
- [Beta Diversity](/manual/beta-diversity.html)
- [Constrained Ordination](/manual/association-biplot.html)
- [ANCOM-BC2](/manual/ancom-bc2.html)
- [MaAsLin2](/manual/maaslin2.html)
- [MaAsLin3](/manual/maaslin3.html)
- [Random Forest](/manual/random-forest.html)
- [SparCC](/manual/sparcc.html)
- [SpiecEasi](/manual/spieceasi.html)
- [Association Heatmap](/manual/correlation-heatmap.html)

## Quick Goal-Based Path
- **Start-to-finish first analysis**:
  - Data Loading -> Preprocessing -> Taxa Profiles -> Alpha/Beta Diversity
- **Differential abundance focus**:
  - ANCOM-BC2, MaAsLin2, or MaAsLin3
- **Prediction focus**:
  - Random Forest
- **Network focus**:
  - SparCC + SpiecEasi
- **Association visualization focus**:
  - Association Heatmap
- **Constrained community analysis**:
  - Constrained Ordination under Diversity

## Notes
- Settings and default values may evolve with UI updates.
- When reporting results, include key parameter settings from each module.

## Cite SimpleMicrobiome

If you use SimpleMicrobiome in your research, please cite the following article.

Na, S.-I., Kim, J., Kim, S.-Y., Park, J., & Cho, Y.-J. (2026). SimpleMicrobiome: An integrated web-based platform for streamlined microbiome data analysis and visualization. *Journal of Microbiology, 64*(9), e2606011. [View article (DOI: 10.71150/jm.2606011)](https://www.jmicrobiol.or.kr/journal/view.php?number=3102)

The app's first screen and Citation tab both display this citation. The Citation tab also lists references for individual analysis methods under **Method References**.
