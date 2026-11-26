(() => {
  const selector = ".heatmap-plot-stage > .shiny-plot-output";
  const pending = new WeakMap();

  function updateAfterDraw(output) {
    const image = output.querySelector("img");
    const src = pending.get(output);
    if (src && image?.getAttribute("src") === src && image.complete && image.naturalWidth > 0) {
      output.parentElement.dataset.heatmapLoading = "false";
      pending.delete(output);
    }
  }

  $(document)
    .on("shiny:outputinvalidated.heatmapLoading", selector, function () {
      pending.delete(this);
      this.parentElement.dataset.heatmapLoading = "true";
    })
    .on("shiny:value.heatmapLoading", selector, function (event) {
      if (!event.value?.src) return;
      pending.set(this, event.value.src);
      setTimeout(() => updateAfterDraw(this), 0);
    })
    .on("shiny:error.heatmapLoading", selector, function () {
      pending.delete(this);
      this.parentElement.dataset.heatmapLoading = "false";
    });

  document.addEventListener("load", (event) => {
    if (!event.target.matches?.("img")) return;
    const output = event.target.closest(selector);
    if (output) updateAfterDraw(output);
  }, true);
})();
