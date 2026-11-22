(() => {
  let draggedItem = null;
  let dropTarget = null;

  function renderItems(editor) {
    const items = editor.querySelector(".plot-order-items");
    const values = editor.querySelector("textarea").value.split("\n")
      .map(value => value.trim()).filter(Boolean);
    const fragment = document.createDocumentFragment();
    for (const value of values) {
      const item = document.createElement("li");
      item.textContent = value;
      item.draggable = true;
      item.tabIndex = 0;
      fragment.append(item);
    }
    items.replaceChildren(fragment);
  }

  function saveDraft(item) {
    const items = item.parentElement;
    const textarea = items.closest(".plot-order-editor").querySelector("textarea");
    const children = Array.from(items.children);
    const index = children.indexOf(item);
    textarea.value = children.map(child => child.textContent).join("\n");
    textarea.dispatchEvent(new Event("change", { bubbles: true }));
    items.children[index]?.focus();
  }

  function clearDropTarget() {
    if (dropTarget) delete dropTarget.dataset.drop;
    dropTarget = null;
  }

  // Shiny textarea updates dispatch a jQuery change event.
  window.jQuery(document).on("input change", ".plot-order-editor textarea", function () {
    renderItems(this.closest(".plot-order-editor"));
  });

  document.addEventListener("dragstart", event => {
    const item = event.target.closest(".plot-order-items > li");
    if (!item) return;
    draggedItem = item;
    event.dataTransfer.effectAllowed = "move";
    event.dataTransfer.setData("text/plain", item.textContent);
  });

  document.addEventListener("dragover", event => {
    clearDropTarget();
    const item = event.target.closest(".plot-order-items > li");
    if (!draggedItem || !item || item === draggedItem ||
        item.parentElement !== draggedItem.parentElement) return;
    event.preventDefault();
    event.dataTransfer.dropEffect = "move";
    const bounds = item.getBoundingClientRect();
    item.dataset.drop = event.clientY < bounds.top + bounds.height / 2 ? "before" : "after";
    dropTarget = item;
  });

  document.addEventListener("drop", event => {
    if (!draggedItem || !dropTarget || !dropTarget.contains(event.target)) return;
    event.preventDefault();
    dropTarget[dropTarget.dataset.drop](draggedItem);
    saveDraft(draggedItem);
    clearDropTarget();
    draggedItem = null;
  });

  document.addEventListener("dragend", () => {
    clearDropTarget();
    draggedItem = null;
  });

  document.addEventListener("keydown", event => {
    const item = event.target.closest(".plot-order-items > li");
    if (!item || !["ArrowUp", "ArrowDown"].includes(event.key)) return;
    event.preventDefault();
    const target = event.key === "ArrowUp" ? item.previousElementSibling : item.nextElementSibling;
    if (!target) return;
    if (event.key === "ArrowUp") target.before(item);
    else target.after(item);
    saveDraft(item);
  });
})();
