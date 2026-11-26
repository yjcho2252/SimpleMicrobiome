(() => {
  const stateSelector = ".profile-group-selector-state";

  function applyState(state, value) {
    const fieldset = state.closest(".profile-group-selectors");
    if (!fieldset) return;
    fieldset.disabled = value !== "enabled";
    fieldset.querySelectorAll("select").forEach((select) => {
      if (!select.selectize) return;
      if (fieldset.disabled) select.selectize.disable();
      else select.selectize.enable();
    });
  }

  $(document).on("shiny:value.profileGroupSelectors", stateSelector, function (event) {
    applyState(this, event.value);
  });

  document.querySelectorAll(stateSelector).forEach((state) => {
    applyState(state, state.textContent);
  });
})();
