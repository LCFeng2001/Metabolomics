## =========================================================
## plot_pca.R
## Publication-level PCA
## =========================================================

run_pca_analysis <- function(expr_scaled, sample_info) {
  expr_scaled <- as_numeric_matrix(expr_scaled)
  sample_info <- as.data.frame(sample_info, check.names = FALSE)

  keep <- apply(expr_scaled, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    is.finite(s) && s > 0
  })

  expr_scaled <- expr_scaled[keep, , drop = FALSE]

  pca <- prcomp(
    t(expr_scaled),
    center = FALSE,
    scale. = FALSE
  )

  var_exp <- pca$sdev^2 / sum(pca$sdev^2)

  scores <- as.data.frame(pca$x, check.names = FALSE)
  scores$sample <- rownames(scores)
  scores <- dplyr::left_join(scores, sample_info, by = "sample")
  scores$time <- factor(scores$time, levels = c("0h", "24h", "48h"))
  scores$treatment <- factor(scores$treatment, levels = c("CK", "Fg"))

  list(
    pca = pca,
    scores = scores,
    var_exp = var_exp
  )
}

plot_pca <- function(
    pca_res,
    title = "PCA_CK_vs_Fg_all_genotypes",
    base_size = 12,
    base_family = "Arial"
) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
  })

  scores <- pca_res$scores
  var_exp <- pca_res$var_exp

  p <- ggplot(
    scores,
    aes(
      x = PC1,
      y = PC2,
      color = treatment,
      shape = time
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
    geom_point(size = 3.1, alpha = 0.92, stroke = 0.75) +
    facet_wrap(~ genotype, scales = "free") +
    scale_color_manual(values = metabo_palette("treatment"), drop = FALSE) +
    labs(
      x = paste0("PC1 (", round(var_exp[1] * 100, 2), "%)"),
      y = paste0("PC2 (", round(var_exp[2] * 100, 2), "%)"),
      color = "Treatment",
      shape = "Time",
      title = title
    ) +
    pub_theme(base_size = base_size, base_family = base_family) +
    theme(
      legend.position = "right",
      panel.grid.major = element_line(linewidth = 0.18, color = "grey90"),
      panel.grid.minor = element_blank()
    )

  p
}
