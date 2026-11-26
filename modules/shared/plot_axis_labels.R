plot_axis_label_size <- function(base_size, default_base_size) {
  if (length(base_size) != 1L || !is.finite(base_size)) base_size <- default_base_size
  10 * base_size / default_base_size
}

plot_axis_label_theme <- function(base_size, default_base_size) {
  ggplot2::theme(axis.text = ggplot2::element_text(
    size = plot_axis_label_size(base_size, default_base_size),
    colour = "black"
  ))
}
