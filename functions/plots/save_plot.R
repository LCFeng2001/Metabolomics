## =========================================================
## save_plot.R
## Publication-level export: PDF, SVG, high-resolution PNG
## =========================================================

save_ggplot_multi <- function(
    p,
    filename_prefix,
    width = 6,
    height = 5,
    dpi = 600
) {
  suppressPackageStartupMessages({
    library(ggplot2)
  })

  if (is.null(p)) return(invisible(NULL))

  dir.create(dirname(filename_prefix), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = paste0(filename_prefix, ".pdf"),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    filename = paste0(filename_prefix, ".png"),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(
      filename = paste0(filename_prefix, ".svg"),
      plot = p,
      width = width,
      height = height,
      units = "in",
      device = svglite::svglite,
      bg = "white",
      limitsize = FALSE
    )
  }

  invisible(NULL)
}
