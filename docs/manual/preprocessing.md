# Preprocessing Manual

## 1. What This Module Does
Applies sample/taxa filtering and selection before all downstream analyses.

## 2. Left Panel Parameters

### 2.1 Sample filtering controls
- Removes samples by quality/coverage criteria.
- Directly changes group balance and statistical power.

### 2.2 Selection state controls
- Sample selection applies when you switch to an analysis tab.
- The selected read counts are passed to downstream modules. Alpha Diversity, Jaccard beta diversity, and Taxa Bar Plot handle subsampling within their own modules.
- The selection summary box reports counts as `Total: N samples` and `Selected: M samples`.

## 3. Plot Order

Open the collapsed **Plot Order** panel in Preprocessing.

Use **Sample Order** to arrange sample axes, including samples within each Taxa Bar Plot panel.

- **Module default** restores each plot's own sorting options.
- **By metadata** sorts by the first selected column, then by the optional second column. For example, choose `SubjectID` and then `Visit` to keep each subject's visits together. Saved group level orders apply to these columns.

Check the sample order preview and click **Apply Sample Order**. The preview shows the first 12 entries for longer lists.

Use **Group Order** to arrange levels of a metadata variable on categorical axes and facets. Select a variable such as `Visit`, drag the levels into the desired order, and click **Apply Group Order**. You can also focus a list item and use the Up or Down key to move it. Each variable has its own setting. **Reset Group Order** restores the selected variable's default.

Numerical axes, categorical regression positions, and heatmap columns with clustering enabled follow their respective analysis settings. Plot order does not change model reference levels or comparison direction.

## 4. Recommended Workflow
1. Set the minimum read count and select the samples to analyze.
2. Confirm the sample and group counts, then open an analysis tab.
