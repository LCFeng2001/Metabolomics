## =========================================================
## plot_vip.R
## Publication-level VIP ranking plot
## v6:
##   - can focus on user-defined classes
##   - can retain only metabolites significant in limma result
##   - used both as standalone figure and as volcano side panel
## =========================================================

plot_vip <- function(
    vip_df,
    annotation = NULL,
    res = NULL,
    top_n = 30,
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
    focus_only = FALSE,
    sig_only = FALSE,
    stat_col = "P.Value",
    pvalue_cutoff = 0.05,
    logFC_cutoff = log2(1.5),
    title = "VIP score",
    base_size = 12,
    base_family = "Arial"
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })

  if (is.null(vip_df) || nrow(vip_df) == 0) {
    warning("Empty VIP table.")
    return(NULL)
  }

  vip_df <- as.data.frame(vip_df, check.names = FALSE)
  vip_df$VIP <- suppressWarnings(as.numeric(vip_df$VIP))
  vip_df <- vip_df[is.finite(vip_df$VIP), , drop = FALSE]
  if (nrow(vip_df) == 0) {
    warning("No finite VIP values.")
    return(NULL)
  }

  if (!is.null(res)) {
    res <- as.data.frame(res, check.names = FALSE)
    if ("metabolite_id" %in% colnames(res)) {
      if (!stat_col %in% colnames(res)) res[[stat_col]] <- NA_real_
      if (!"logFC" %in% colnames(res)) res$logFC <- NA_real_
      res$logFC <- suppressWarnings(as.numeric(res$logFC))
      res[[stat_col]] <- suppressWarnings(as.numeric(res[[stat_col]]))
      res <- add_focus_class_columns(
        res,
        focus_classes = focus_classes,
        class_col = "Class",
        bx_patterns = focus_bx_patterns
      )
      res$vip_keep <- TRUE
      if (isTRUE(focus_only)) {
        res$vip_keep <- res$vip_keep & res$is_focus_class
      }
      if (isTRUE(sig_only)) {
        res$vip_keep <- res$vip_keep & is.finite(res[[stat_col]]) & res[[stat_col]] < pvalue_cutoff &
          is.finite(res$logFC) & abs(res$logFC) >= logFC_cutoff
      }
      keep_cols <- c("metabolite_id", "logFC", stat_col, "focus_class", "is_focus_class", "vip_keep", "Class")
      vip_df_all <- dplyr::left_join(vip_df, res[, keep_cols, drop = FALSE], by = "metabolite_id")
      if (isTRUE(focus_only) || isTRUE(sig_only)) {
        vip_df <- vip_df_all[!is.na(vip_df_all$vip_keep) & vip_df_all$vip_keep, , drop = FALSE]
        if (nrow(vip_df) == 0) {
          vip_df <- vip_df_all
        }
      } else {
        vip_df <- vip_df_all
      }
    }
  }

  ## fallback if filtering removed all rows
  if (nrow(vip_df) == 0) {
    warning("No VIP metabolites remained after filtering.")
    return(NULL)
  }

  if (!is.null(annotation)) {
    annotation <- as.data.frame(annotation, check.names = FALSE)
    keep_anno <- intersect(c("metabolite_id", "Compound name", "Class"), colnames(annotation))
    if (length(keep_anno) >= 2) {
      vip_df <- dplyr::left_join(vip_df, annotation[, keep_anno, drop = FALSE], by = "metabolite_id")
    }
  }

  if (!"Compound name" %in% colnames(vip_df)) vip_df[["Compound name"]] <- NA_character_
  if (!"focus_class" %in% colnames(vip_df)) {
    vip_df <- add_focus_class_columns(vip_df, focus_classes = focus_classes, class_col = "Class", bx_patterns = focus_bx_patterns)
  }

  vip_df$label <- ifelse(
    is.na(vip_df[["Compound name"]]) | vip_df[["Compound name"]] == "",
    vip_df$metabolite_id,
    vip_df[["Compound name"]]
  )

  vip_df$display_class <- vip_df$focus_class
  vip_df$display_class[is.na(vip_df$display_class) | vip_df$display_class == ""] <- "Other"

  vip_df <- vip_df[order(vip_df$VIP, decreasing = TRUE), , drop = FALSE]
  vip_df <- vip_df[seq_len(min(top_n, nrow(vip_df))), , drop = FALSE]

  vip_df$label <- paste0(vip_df$label, ifelse(vip_df$display_class != "Other", paste0(" [", vip_df$display_class, "]"), ""))
  vip_df$label <- factor(vip_df$label, levels = rev(vip_df$label))

  pal <- c(metabo_palette("focus_class"), Other = "grey55")
  miss <- setdiff(unique(vip_df$display_class), names(pal))
  if (length(miss) > 0) {
    extra <- rep("grey55", length(miss))
    names(extra) <- miss
    pal <- c(pal, extra)
  }

  p <- ggplot(vip_df, aes(x = VIP, y = label, fill = display_class)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.42, color = "grey35") +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(
      x = "VIP score",
      y = NULL,
      fill = "Class",
      title = title
    ) +
    pub_theme(base_size = base_size, base_family = base_family) +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.18, color = "grey90"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = base_size + 1)
    )

  attr(p, "n_bars") <- nrow(vip_df)
  p
}
