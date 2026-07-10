## =========================================================
## plot_volcano.R
## Volcano plot in the style requested by the user
## v24:
##   - style based on the provided ggplot code
##   - classic theme
##   - Up / NS / Down colors
##   - focus metabolites shown with a different shape
##   - VIP retained as point size
##   - geom_text_repel settings tuned to always show connector lines
## =========================================================

plot_volcano <- function(
    res,
    vip_df = NULL,
    logFC_cutoff = 1,
    pvalue_cutoff = 0.05,
    stat_col = "P.Value",
    y_label = "Pvalue",
    focus_classes = c(
      "Phenolamides",
      "Benzoic acid derivatives",
      "Hydroxycinnamoyl derivatives",
      "Flavone",
      "Flavonol",
      "Flavanone",
      "Flavone C-glycosides",
      "BX"
    ),
    focus_bx_patterns = c("^BX$", "Benzoxazinoid", "Benzoxazinoids"),
    label_focus_only = TRUE,
    label_sig_focus_only = TRUE,
    label_all_focus = TRUE,
    label_top = 5,
    title = "Volcano plot",
    base_size = 12,
    base_family = "Arial"
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
  })

  res <- as.data.frame(res, check.names = FALSE)

  if (!"metabolite_id" %in% colnames(res)) stop("res must contain metabolite_id.")
  if (!"logFC" %in% colnames(res)) stop("res must contain logFC.")
  if (!stat_col %in% colnames(res)) stop("res must contain ", stat_col, ".")

  res$logFC <- suppressWarnings(as.numeric(res$logFC))
  res$stat_p <- suppressWarnings(as.numeric(res[[stat_col]]))

  if (!is.null(vip_df)) {
    vip_df <- as.data.frame(vip_df, check.names = FALSE)
    if (all(c("metabolite_id", "VIP") %in% colnames(vip_df))) {
      vip_df$VIP <- suppressWarnings(as.numeric(vip_df$VIP))
      vip_df <- vip_df[, intersect(c("metabolite_id", "VIP", "VIP_method"), colnames(vip_df)), drop = FALSE]
      res <- dplyr::left_join(res, vip_df, by = "metabolite_id")
    }
  }
  if (!"VIP" %in% colnames(res)) res$VIP <- NA_real_

  res <- add_focus_class_columns(
    res,
    focus_classes = focus_classes,
    class_col = "Class",
    bx_patterns = focus_bx_patterns
  )

  min_p <- suppressWarnings(min(res$stat_p[is.finite(res$stat_p) & res$stat_p > 0], na.rm = TRUE))
  if (!is.finite(min_p)) min_p <- 1e-300

  res$stat_p_plot <- res$stat_p
  res$stat_p_plot[!is.finite(res$stat_p_plot)] <- NA_real_
  res$stat_p_plot[res$stat_p_plot <= 0] <- min_p / 10
  res$v <- -log10(res$stat_p_plot)

  res$group <- dplyr::case_when(
    is.finite(res$stat_p) & res$stat_p < pvalue_cutoff & res$logFC >= logFC_cutoff ~ "Up",
    is.finite(res$stat_p) & res$stat_p < pvalue_cutoff & res$logFC <= -logFC_cutoff ~ "Down",
    TRUE ~ "NS"
  )

  res$Shape <- ifelse(res$is_focus_class, "Focus metabolites", "Other metabolites")

  if ("Compound name" %in% colnames(res)) {
    res$ID <- ifelse(
      is.na(res[["Compound name"]]) | res[["Compound name"]] == "",
      res$metabolite_id,
      res[["Compound name"]]
    )
  } else {
    res$ID <- res$metabolite_id
  }

  ## VIP -> point size, but keep it in a reasonable range
  res$VIP_plot <- suppressWarnings(as.numeric(res$VIP))
  res$VIP_plot[!is.finite(res$VIP_plot)] <- 0
  if (max(res$VIP_plot, na.rm = TRUE) <= 0) {
    res$VIP_plot <- rep(1, nrow(res))
    vip_breaks <- c(1)
  } else {
    pos_vip <- res$VIP_plot[res$VIP_plot > 0]
    lo <- if (length(pos_vip) > 0) quantile(pos_vip, 0.05, na.rm = TRUE) else 0.2
    hi <- if (length(pos_vip) > 0) quantile(pos_vip, 0.95, na.rm = TRUE) else max(res$VIP_plot, na.rm = TRUE)
    lo <- as.numeric(lo); hi <- as.numeric(hi)
    if (!is.finite(lo) || lo <= 0) lo <- 0.2
    if (!is.finite(hi) || hi <= lo) hi <- max(res$VIP_plot, na.rm = TRUE)
    res$VIP_plot[res$VIP_plot <= 0] <- lo
    res$VIP_plot <- pmin(res$VIP_plot, hi)
    vip_breaks <- pretty(res$VIP_plot, n = 4)
    vip_breaks <- vip_breaks[vip_breaks > 0]
  }

  ## labels
  labels_df <- res
  if (isTRUE(label_focus_only)) {
    labels_df <- labels_df[labels_df$is_focus_class, , drop = FALSE]
  }
  if (isTRUE(label_sig_focus_only)) {
    labels_df <- labels_df[labels_df$group %in% c("Up", "Down"), , drop = FALSE]
  }
  labels_df <- labels_df[is.finite(labels_df$v), , drop = FALSE]
  labels_df <- labels_df[order(labels_df$stat_p, -abs(labels_df$logFC), -labels_df$VIP_plot), , drop = FALSE]
  if (is.finite(label_top)) {
    labels_df <- labels_df[seq_len(min(label_top, nrow(labels_df))), , drop = FALSE]
  }

  custom_colors <- c(
    "Up" = "#E64B35",
    "NS" = "#CCCCCC",
    "Down" = "#4DBBD5"
  )
  target_order <- c("Up", "NS", "Down")

  res$group <- factor(res$group, levels = target_order)
  res$Shape <- factor(res$Shape, levels = c("Other metabolites", "Focus metabolites"))

  axis_title_x <- "log2(Fold Change)"
  axis_title_y <- "-log10(Pvalue)"

  use_markdown <- requireNamespace("ggtext", quietly = TRUE)
  if (use_markdown) {
    axis_title_x <- "log<sub>2</sub>(Fold Change)"
    axis_title_y <- "-log<sub>10</sub>(Pvalue)"
  }

  volcano_plot_advanced <- ggplot(res, aes(x = logFC, y = v, color = group, shape = Shape, size = VIP_plot)) +
    geom_point(alpha = 0.8) +

    scale_color_manual(values = custom_colors, limits = target_order) +
    scale_shape_manual(values = c("Other metabolites" = 16, "Focus metabolites" = 17)) +
    scale_size_continuous(
      range = c(0.35, 3.2),
      breaks = vip_breaks,
      name = "VIP value"
    ) +

    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +

    geom_text_repel(
      data = labels_df,
      aes(label = ID),
      size = 3.5,
      family = base_family,
      box.padding = 0.6,
      point.padding = 0.3,
      min.segment.length = 0,
      segment.color = "grey40",
      segment.linewidth = 0.4,
      segment.alpha = 0.8,
      show.legend = FALSE
    ) +

    labs(
      title = title,
      x = axis_title_x,
      y = axis_title_y
    ) +

    theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = if (use_markdown) ggtext::element_markdown(size = 12, face = "bold") else element_text(size = 12, face = "bold"),
      axis.title.y = if (use_markdown) ggtext::element_markdown(size = 12, face = "bold") else element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      legend.text = element_text(size = 11),
      legend.title = element_blank(),
      legend.position = "right"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 3)),
      shape = guide_legend(order = 2, override.aes = list(size = 3, color = "black")),
      size = guide_legend(order = 3, title = "VIP value")
    )

  volcano_plot_advanced
}
