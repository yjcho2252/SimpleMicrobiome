additional_plot_palettes <- c(
  "Okabe–Ito" = "okabe_ito",
  "Pastel" = "pastel",
  "Dynamic" = "dynamic",
  "Warm" = "warm",
  "Cool" = "cool"
)

group_plot_palettes <- c("Set2" = "set2", "Dark2" = "dark2", "Paired" = "paired", "Grayscale" = "gray",
                         additional_plot_palettes)

group_plot_colors <- function(n, palette) {
  if (n == 0L) return(character())
  if (length(palette) != 1L || is.na(palette) || !palette %in% group_plot_palettes) palette <- "set2"
  if (palette == "gray") return(grDevices::gray.colors(n, start = 0.2, end = 0.85))
  if (palette == "okabe_ito") {
    colors <- grDevices::palette.colors(min(n, 8L), palette = "Okabe-Ito")
    if (n > 8L) colors <- c(colors, grDevices::hcl.colors(n - 8L, "Dynamic"))
    return(unname(colors))
  }
  hcl_palette <- c(set2 = "Set 2", dark2 = "Dark 3", paired = "Set 3",
                   pastel = "Pastel 1", dynamic = "Dynamic", warm = "Warm", cool = "Cold")
  grDevices::hcl.colors(n, palette = hcl_palette[[palette]])
}
