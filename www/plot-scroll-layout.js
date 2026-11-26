(() => {
  let frame;

  function updateHeight() {
    cancelAnimationFrame(frame);
    frame = requestAnimationFrame(() => {
      const viewportHeight = document.documentElement.clientHeight;
      const bottomGap = parseFloat(getComputedStyle(document.documentElement).fontSize);

      document.querySelectorAll(".plot-scroll-layout").forEach((layout) => {
        if (!layout.getClientRects().length) return;
        const top = layout.getBoundingClientRect().top + window.scrollY;
        const available = Math.max(0, Math.floor(viewportHeight - top - bottomGap));
        layout.style.setProperty("--plot-scroll-height", `${available}px`);
      });
    });
  }

  window.addEventListener("resize", updateHeight);
  document.addEventListener("shown.bs.tab", updateHeight);
  const navbar = document.querySelector(".navbar");
  if (navbar && typeof ResizeObserver !== "undefined") {
    new ResizeObserver(updateHeight).observe(navbar);
  }
  updateHeight();
})();
