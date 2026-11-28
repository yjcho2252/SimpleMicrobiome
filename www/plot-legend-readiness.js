(() => {
  const legends = ".app-plot-legend[data-legend-results]";
  const imageOutputs = ".shiny-plot-output, .shiny-image-output";
  const observed = new WeakSet();
  const expectedImages = new WeakMap();
  const pendingWidgets = new WeakSet();
  const resultIds = new Set();
  let frame = null;

  const observer = new MutationObserver(schedule);

  function isReady(result) {
    if (!result.getClientRects().length ||
        result.matches(".recalculating, .shiny-output-error")) {
      return false;
    }
    if (result.matches(imageOutputs)) {
      const image = result.querySelector("img");
      const expected = expectedImages.get(result);
      return image?.complete && image.naturalWidth > 0 &&
        (!expected || image.getAttribute("src") === expected);
    }
    if (result.matches(".plotly.html-widget-output")) {
      return result.hasAttribute("data-taxa-ready") && !pendingWidgets.has(result);
    }
    return result.childElementCount > 0 || result.textContent.trim() !== "";
  }

  function schedule() {
    if (frame !== null) return;
    frame = requestAnimationFrame(() => {
      frame = null;
      resultIds.clear();
      document.querySelectorAll(legends).forEach((legend) => {
        const results = legend.dataset.legendResults.split(" ").map((id) => {
          resultIds.add(id);
          return document.getElementById(id);
        }).filter(Boolean);
        results.forEach((result) => {
          if (observed.has(result)) return;
          observed.add(result);
          observer.observe(result, {
            attributes: true,
            attributeFilter: ["class", "data-taxa-ready"],
            childList: true
          });
        });
        legend.toggleAttribute("data-legend-ready", results.some(isReady));
      });
    });
  }

  function relevant(target) {
    return resultIds.has(target.id) || target.matches?.(legends);
  }

  $(document)
    .on("shiny:value.plotLegend", (event) => {
      const result = event.target;
      if (!relevant(result)) return;
      if (result.matches?.(imageOutputs) && event.value?.src) {
        expectedImages.set(result, event.value.src);
      }
      if (result.matches?.(".plotly.html-widget-output")) {
        pendingWidgets.add(result);
      }
      schedule();
    })
    .on("shiny:bound.plotLegend shiny:unbound.plotLegend shiny:error.plotLegend shiny:outputinvalidated.plotLegend shiny:visualchange.plotLegend", (event) => {
      if (relevant(event.target)) schedule();
    });

  document.addEventListener("plotlegend:widgetready", (event) => {
    pendingWidgets.delete(event.target);
    if (relevant(event.target)) schedule();
  });
  document.addEventListener("shown.bs.tab", schedule);
  function onImageEvent(event) {
    if (event.target.matches?.("img") &&
        relevant(event.target.closest(imageOutputs) || event.target)) schedule();
  }
  document.addEventListener("load", onImageEvent, true);
  document.addEventListener("error", onImageEvent, true);
  schedule();
})();
