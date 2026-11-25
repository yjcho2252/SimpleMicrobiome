add_subject_pair_lines <- function(plot, points, metadata, subject_col, x_col, y_col) {
  if (length(subject_col) != 1L || !subject_col %in% names(metadata) ||
      !all(c("SampleID", x_col, y_col) %in% names(points))) {
    return(plot)
  }

  subject <- trimws(as.character(metadata[[subject_col]][match(points$SampleID, rownames(metadata))]))
  keep <- !is.na(subject) & nzchar(subject) &
    is.finite(points[[x_col]]) & is.finite(points[[y_col]])
  rows_by_subject <- split(which(keep), subject[keep])
  pairs <- rows_by_subject[lengths(rows_by_subject) == 2L]
  if (!length(pairs)) return(plot)

  first <- vapply(pairs, function(rows) rows[1L], integer(1))
  second <- vapply(pairs, function(rows) rows[2L], integer(1))
  segments <- data.frame(
    x = points[[x_col]][first], y = points[[y_col]][first],
    xend = points[[x_col]][second], yend = points[[y_col]][second]
  )
  layer <- ggplot2::geom_segment(
    data = segments, mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE, colour = "grey55", linewidth = 0.45, alpha = 0.75
  )
  plot$layers <- c(list(layer), plot$layers)
  plot
}
