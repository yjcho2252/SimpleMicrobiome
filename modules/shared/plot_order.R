plot_order_levels <- function(default, preferred) {
  c(intersect(preferred, default), setdiff(default, preferred))
}

plot_order_factor <- function(values, variable, plot_order = NULL, default = levels(factor(values))) {
  preferred <- if (identical(variable, "SampleID")) {
    plot_order$samples
  } else if (length(variable) == 1L && !is.na(variable)) {
    groups <- plot_order$groups
    if (!variable %in% names(groups) && length(groups)) {
      matches <- names(groups)[make.names(names(groups)) == make.names(variable)]
      if (length(matches) == 1L) variable <- matches
    }
    groups[[variable]]
  }
  factor(values, levels = plot_order_levels(default, preferred))
}
