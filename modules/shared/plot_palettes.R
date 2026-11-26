additional_plot_palettes <- c(
  "Okabe–Ito" = "okabe_ito",
  "Pastel" = "pastel",
  "Dynamic" = "dynamic",
  "Warm" = "warm",
  "Cool" = "cool"
)

group_plot_palettes <- c("Set2" = "set2", "Dark2" = "dark2", "Paired" = "paired", "Grayscale" = "gray",
                         additional_plot_palettes)

brewer_plot_palettes <- c(set2 = "Set2", dark2 = "Dark2", paired = "Paired")

brewer_plot_colors <- function(n, palette) {
  if (n == 0L) return(character())
  limit <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
  colors <- RColorBrewer::brewer.pal(max(3L, min(n, limit)), palette)
  if (n <= limit) return(colors[seq_len(n)])
  grDevices::colorRampPalette(colors)(n)
}

group_plot_colors <- function(n, palette) {
  if (n == 0L) return(character())
  if (length(palette) != 1L || is.na(palette) || !palette %in% group_plot_palettes) palette <- "set2"
  if (palette == "gray") return(grDevices::gray.colors(n, start = 0.2, end = 0.85))
  if (palette %in% names(brewer_plot_palettes)) {
    return(brewer_plot_colors(n, brewer_plot_palettes[[palette]]))
  }
  if (palette == "okabe_ito") {
    colors <- grDevices::palette.colors(min(n, 8L), palette = "Okabe-Ito")
    if (n > 8L) colors <- c(colors, grDevices::hcl.colors(n - 8L, "Dynamic"))
    return(unname(colors))
  }
  hcl_palette <- c(pastel = "Pastel 1", dynamic = "Dynamic", warm = "Warm", cool = "Cold")
  grDevices::hcl.colors(n, palette = hcl_palette[[palette]])
}

diverging_plot_palettes <- c(
  "Blue - White - Red" = "blue_red",
  "Blue - White - Orange (colorblind friendly)" = "blue_orange",
  "Purple - White - Green" = "purple_green",
  "Teal - White - Brown" = "teal_brown"
)

diverging_plot_colors <- function(palette) {
  if (length(palette) != 1L || !palette %in% diverging_plot_palettes) palette <- "blue_red"
  switch(
    palette,
    blue_red = c("#2166AC", "#F7F7F7", "#B2182B"),
    blue_orange = c("#3B4CC0", "#F7F7F7", "#E66101"),
    purple_green = c("#762A83", "#F7F7F7", "#1B7837"),
    teal_brown = c("#018571", "#F5F5F5", "#A6611A")
  )
}

sequential_plot_palettes <- c(
  "Viridis" = "viridis", "Cividis" = "cividis", "Plasma" = "plasma",
  "Magma" = "magma", "Inferno" = "inferno",
  "Light to dark blue" = "blue", "Light to dark red" = "red",
  "Light to dark orange" = "orange"
)

sequential_plot_colors <- function(palette) {
  if (length(palette) != 1L || !palette %in% sequential_plot_palettes) palette <- "plasma"
  switch(
    palette,
    blue = c("#DEEBF7", "#2166AC"),
    red = c("#FDE0DD", "#B2182B"),
    orange = c("#FEE6CE", "#E66101"),
    viridisLite::viridis(256L, option = palette)
  )
}
