# MaAsLin3 Manual

## 1. What This Module Does

MaAsLin3 tests microbial abundance and prevalence associations with metadata. Abundance describes how much of a feature is observed among samples where it is detected. Prevalence describes the likelihood of detection. Both models are fitted when you click **Run MaAsLin3**.

## 2. Left Panel

1. **Primary variable**. Choose the metadata variable of interest. For a numeric variable, **Treat as** defaults to **Continuous**. Choose **Categorical** for numeric group codes such as 0 and 5.
2. **Comparison groups (levels)**. For categorical predictors, all observed levels start selected. Keep at least two. A warning appears above 10 observed levels, and you can still select and analyze the variable.
3. **Reference level**. Choose a baseline from the included levels for ordinary categorical contrasts. Group tests and ordered contrasts use their own comparison rules.
4. **Taxonomic level**. Choose a taxonomic rank or **Feature level** to analyze individual input features.
5. **Significance threshold (q)**. Sets the significance threshold for the summary plot and the table's significance filter. The default is 0.1.
6. **Prevalence filter (%)**. Sets the minimum detection frequency for features entering the models. The default is 5%.

### Advanced Options

7. **Additional covariates**. Add metadata variables to adjust for.
8. **Interaction terms**. Select pairs of included fixed effects. Variables used in group tests or ordered contrasts cannot also enter an interaction.
9. **Random intercepts**. Select one or more grouping variables with repeated observations.
10. **Adjust for input read depth**. Adds sample total reads as a covariate when depths differ. This is enabled by default. The covariate is skipped when depth is constant, including after fixed-depth subsampling.
11. **Group tests**. Test a categorical fixed effect as a whole. An overall test has no single directional effect size.
12. **Ordered contrasts**. Select a categorical fixed effect and supply its levels from lowest to highest. The model compares adjacent levels.
13. **Matched strata**. Select one matching variable. The primary predictor must vary within at least one retained stratum. This setting cannot be combined with random intercepts.

The preview lists selected terms, read depth adjustment, and the sample count used for analysis. The module uses TSS normalization, log transformation, and standardized numeric covariates. Abundance tests use a median coefficient comparison, and prevalence tests use zero as the null.

### Plot Dimensions

Set the width and height in pixels and the base font size in points. The controls apply to Summary Plot and Association Detail. For Association Detail, the two panels share the selected width. Plot settings can be changed after running the model.

## 3. Results

### Summary Plot

The package-generated summary shows coefficient plots for up to two leading comparison terms by default and a heatmap for remaining terms where available. Heatmap stars indicate small individual q-values. The preview is PNG, and exports are available as PNG or vector PDF using the current plot settings.

### Results Table

Switch among **Abundance**, **Prevalence**, and **All results**. Filter by the primary variable or by significance using individual or joint q-values. Individual q-values apply to a specific model, while joint q-values combine evidence across both models. The full package result is available as a TSV download.

The status panel is below the result tabs. Numeric table values are displayed to two decimal places.

Abundance effect sizes are relative abundance fold changes, and prevalence effect sizes are detection odds ratios. For standardized continuous predictors, effects are per standard deviation. Review model fitting warnings before interpreting an estimate.

### Association Detail

Select a table row for an ordinary main effect to see observed sample values alongside the adjusted model effect and individual q-value.

| Predictor | Abundance | Prevalence |
|---|---|---|
| Categorical | Boxplot with detected sample points | Detection percentage and detected/total count by group |
| Continuous | Scatter plot for detected samples | Predictor distribution for detected and undetected samples |

Interaction estimates, group tests, and ordered contrasts are available in the results table. Changing input data or model options requires a new run.

### Volcano Plot

Abundance and prevalence are shown in separate panels. The X-axis is the coefficient minus its tested null value, and the Y-axis is the negative base-10 logarithm of the individual q-value. The plot can be downloaded as PNG.

See the [MaAsLin3 1.2.0 documentation](https://bioconductor.posit.co/packages/3.22/bioc/vignettes/maaslin3/inst/doc/maaslin3_manual.html) for additional model and output details.
