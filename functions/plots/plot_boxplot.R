## =========================================================
## plot_boxplot.R
## Publication-level box/jitter plots for significant focus-class metabolites
## v6:
##   - show differential metabolites from user-defined focus classes
##   - default rule uses nominal P.Value + |logFC| cutoff
##   - can display all significant focus metabolites instead of generic top_n
## =========================================================

plot_boxplot_top <- function(
    expr_log,
    res,
    sample_info,
    comparison_row,
    annotation = NULL,
    top_n = 12,
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
    stat_col = "P.Value",
    pvalue_cutoff = 0.05,
    logFC_cutoff = log2(1.5),
    focus_only = TRUE,
    show_all_sig_focus = TRUE,
    max_facets = Inf,
    facet_ncol = 4,
    title = "Boxplot",
    base_size = 12,
    base_family = "Arial"
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })

  expr_log <- as_numeric_matrix(expr_log)
  res <- as.data.frame(res, check.names = FALSE)
  sample_info <- as.data.frame(sample_info, check.names = FALSE)

  if (!"metabolite_id" %in% colnames(res)) stop("res must contain metabolite_id")
  if (!"logFC" %in% colnames(res)) stop("res must contain logFC")
  if (!stat_col %in% colnames(res)) stop("res must contain ", stat_col)

  res$logFC <- suppressWarnings(as.numeric(res$logFC))
  res[[stat_col]] <- suppressWarnings(as.numeric(res[[stat_col]]))

  res <- add_focus_class_columns(
    res,
    focus_classes = focus_classes,
    class_col = "Class",
    bx_patterns = focus_bx_patterns
  )

  stat_vec <- res[[stat_col]]
  res$selected_for_boxplot <- is.finite(res$logFC) &
    is.finite(stat_vec) &
    stat_vec < pvalue_cutoff &
    abs(res$logFC) >= logFC_cutoff

  res2 <- res
  if (isTRUE(focus_only)) {
    res2 <- res2[res2$is_focus_class, , drop = FALSE]
  }

  if (isTRUE(show_all_sig_focus)) {
    res2 <- res2[res2$selected_for_boxplot, , drop = FALSE]
  } else {
    res2 <- res2[order(res2[[stat_col]], -abs(res2$logFC)), , drop = FALSE]
    if (is.finite(max_facets)) {
      res2 <- res2[seq_len(min(max_facets, nrow(res2))), , drop = FALSE]
    } else {
      res2 <- res2[seq_len(min(top_n, nrow(res2))), , drop = FALSE]
    }
  }

  res2 <- res2[order(res2[[stat_col]], -abs(res2$logFC)), , drop = FALSE]

  if (nrow(res2) < 1) {
    warning("No significant focus-class metabolites for boxplot: ", comparison_row$comparison)
    return(NULL)
  }

  if (is.finite(max_facets)) {
    res2 <- res2[seq_len(min(max_facets, nrow(res2))), , drop = FALSE]
  }

  top_ids <- intersect(unique(res2$metabolite_id), rownames(expr_log))
  if (length(top_ids) < 1) {
    warning("Selected metabolites not found in expr_log: ", comparison_row$comparison)
    return(NULL)
  }

  selected_samples <- sample_info$sample[
    sample_info$group %in% c(comparison_row$control_group, comparison_row$case_group)
  ]
  selected_samples <- intersect(selected_samples, colnames(expr_log))

  if (length(selected_samples) < 2) {
    warning("Too few samples for boxplot: ", comparison_row$comparison)
    return(NULL)
  }

  mat <- expr_log[top_ids, selected_samples, drop = FALSE]
  mat <- as_numeric_matrix(mat)
  mat[!is.finite(mat)] <- NA_real_

  df <- as.data.frame(mat, check.names = FALSE) %>%
    tibble::rownames_to_column("metabolite_id") %>%
    tidyr::pivot_longer(
      cols = -metabolite_id,
      names_to = "sample",
      values_to = "value"
    ) %>%
    left_join(sample_info, by = "sample") %>%
    left_join(
      res2[, c("metabolite_id", "logFC", stat_col, "focus_class", "Class", "selected_for_boxplot"), drop = FALSE],
      by = "metabolite_id"
    )

  if (!is.null(annotation)) {
    annotation <- as.data.frame(annotation, check.names = FALSE)
    if (all(c("metabolite_id", "Compound name") %in% colnames(annotation))) {
      df <- df %>% left_join(
        annotation[, c("metabolite_id", "Compound name"), drop = FALSE],
        by = "metabolite_id"
      )
    }
  }

  if (!"Compound name" %in% colnames(df)) df[["Compound name"]] <- NA_character_

  df$compound_label <- ifelse(
    is.na(df[["Compound name"]]) | df[["Compound name"]] == "",
    df$metabolite_id,
    df[["Compound name"]]
  )
  df$focus_class[is.na(df$focus_class) | df$focus_class == ""] <- normalize_class_text(df$Class)

  meta_order <- unique(res2$metabolite_id)
  meta_table <- unique(df[, c("metabolite_id", "compound_label", "focus_class", stat_col, "logFC"), drop = FALSE])
  meta_table <- meta_table[match(meta_order, meta_table$metabolite_id), , drop = FALSE]
  meta_table$facet_label <- paste0(meta_table$compound_label, "\n[", meta_table$focus_class, "]")

  df <- df %>% left_join(meta_table[, c("metabolite_id", "facet_label"), drop = FALSE], by = "metabolite_id")
  df$treatment <- factor(
    df$treatment,
    levels = c(comparison_row$control_treatment, comparison_row$case_treatment)
  )
  df$facet_label <- factor(df$facet_label, levels = meta_table$facet_label)

  p <- ggplot(df, aes(x = treatment, y = value, fill = treatment, color = treatment)) +
    geom_boxplot(
      width = 0.56,
      outlier.shape = NA,
      alpha = 0.24,
      linewidth = 0.48
    ) +
    geom_point(
      position = position_jitter(width = 0.13, height = 0, seed = 1),
      size = 1.8,
      alpha = 0.88
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 2.5,
      fill = "white",
      color = "black",
      stroke = 0.45
    ) +
    facet_wrap(~ facet_label, scales = "free_y", ncol = facet_ncol) +
    scale_fill_manual(values = metabo_palette("treatment"), drop = FALSE) +
    scale_color_manual(values = metabo_palette("treatment"), drop = FALSE) +
    labs(
      x = NULL,
      y = expression(log[2]~normalized~Peak~Area),
      title = title
    ) +
    pub_theme(base_size = base_size, base_family = base_family) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = base_size - 2, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )

  attr(p, "n_facets") <- length(unique(df$metabolite_id))
  attr(p, "facet_ncol") <- facet_ncol
  p
}
