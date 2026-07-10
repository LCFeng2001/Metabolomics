## =========================================================
## 07_Plots.R
## Publication-level plots
## v6:
##   - total heatmap: all metabolites x all samples
##   - volcano uses P.Value and focus classes
##   - boxplot shows significant focus-class metabolites
##   - VIP ranking is generated and combined with volcano figure
## =========================================================

source("config/config.R")
source("functions/utils/common_utils.R")
source("functions/plots/save_plot.R")
source("functions/plots/plot_volcano.R")
source("functions/plots/plot_heatmap.R")
source("functions/plots/plot_boxplot.R")
source("functions/plots/plot_vip.R")

data01 <- readRDS(file.path(dir_rds, "01_clean_aggregated.rds"))
data02 <- readRDS(file.path(dir_rds, "02_preprocessed.rds"))
comparisons <- readRDS(file.path(dir_rds, "04_comparisons.rds"))
limma_results <- readRDS(file.path(dir_rds, "05_limma_results.rds"))
opls_results_file <- file.path(dir_rds, "06_oplsda_results.rds")
opls_results <- if (file.exists(opls_results_file)) readRDS(opls_results_file) else list()

format_comparison_title <- function(comp_row) {
  ## Use the exact comparison name as figure title, e.g. M03_48h_Fg_vs_CK
  if ("comparison" %in% colnames(comp_row)) {
    nm <- as.character(comp_row$comparison)
    if (!is.na(nm) && nm != "") return(nm)
  }

  genotype <- if ("genotype" %in% colnames(comp_row)) as.character(comp_row$genotype) else NA_character_
  time <- if ("time" %in% colnames(comp_row)) as.character(comp_row$time) else NA_character_
  case_treat <- if ("case_treatment" %in% colnames(comp_row)) as.character(comp_row$case_treatment) else NA_character_
  control_treat <- if ("control_treatment" %in% colnames(comp_row)) as.character(comp_row$control_treatment) else NA_character_

  parts <- c(genotype, time, case_treat, "vs", control_treat)
  parts <- parts[!is.na(parts) & parts != ""]
  paste(parts, collapse = "_")
}

get_heatmap_matrix_and_annotation <- function(mode = heatmap_use_matrix) {
  mode <- as.character(mode)

  if (mode == "raw") {
    return(list(
      mat = data01$peak_area,
      annotation = data01$annotation,
      sample_info = data01$sample_info,
      value_label = "Raw\nPeak area"
    ))
  }

  if (mode == "normalized") {
    return(list(
      mat = data02$normalized_area,
      annotation = data02$annotation,
      sample_info = data02$sample_info,
      value_label = "Normalized\nPeak area"
    ))
  }

  if (mode == "log") {
    return(list(
      mat = data02$expr_log,
      annotation = data02$annotation,
      sample_info = data02$sample_info,
      value_label = "log2\nintensity"
    ))
  }

  list(
    mat = data02$expr_scaled,
    annotation = data02$annotation,
    sample_info = data02$sample_info,
    value_label = "Pareto-scaled\nvalue"
  )
}


## ---------------------------------------------------------
## Heatmap output folders
## ---------------------------------------------------------
dir_heatmap_total <- file.path(dir_heatmap, "total")
dir_heatmap_differential <- file.path(dir_heatmap, "differential")
dir.create(dir_heatmap_total, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_heatmap_differential, recursive = TRUE, showWarnings = FALSE)

## remove old root-level heatmap outputs from earlier versions to avoid confusion
old_heatmap_files <- list.files(
  dir_heatmap,
  pattern = "^(Total|Differential)_ComplexHeatmap.*\\.(pdf|png|svg|xlsx)$",
  full.names = TRUE
)
if (length(old_heatmap_files) > 0) {
  suppressWarnings(file.remove(old_heatmap_files))
}

make_differential_heatmap_table <- function(limma_results, stat_col, pvalue_cutoff, logFC_cutoff, annotation = NULL) {
  out <- lapply(names(limma_results), function(nm) {
    res <- as.data.frame(limma_results[[nm]], check.names = FALSE)
    if (!"metabolite_id" %in% colnames(res)) return(data.frame())
    if (!"logFC" %in% colnames(res)) return(data.frame())

    pcol <- if (stat_col %in% colnames(res)) stat_col else if ("P.Value" %in% colnames(res)) "P.Value" else NA_character_
    if (is.na(pcol)) return(data.frame())

    res$logFC <- suppressWarnings(as.numeric(res$logFC))
    res[[pcol]] <- suppressWarnings(as.numeric(res[[pcol]]))

    keep <- which(
      is.finite(res$logFC) &
        is.finite(res[[pcol]]) &
        abs(res$logFC) >= logFC_cutoff &
        res[[pcol]] < pvalue_cutoff
    )

    if (length(keep) == 0) return(data.frame())

    tmp <- res[keep, , drop = FALSE]
    tmp$comparison <- nm
    tmp$stat_col_used <- pcol
    tmp$stat_value_used <- tmp[[pcol]]
    tmp$direction <- ifelse(tmp$logFC > 0, "Up", "Down")
    tmp
  })

  diff_table <- dplyr::bind_rows(out)

  if (nrow(diff_table) > 0 && !is.null(annotation)) {
    annotation <- as.data.frame(annotation, check.names = FALSE)
    keep_anno_cols <- intersect(
      c("metabolite_id", "Compound name", "Class", "source_features"),
      colnames(annotation)
    )
    if ("metabolite_id" %in% keep_anno_cols) {
      diff_table <- dplyr::left_join(
        diff_table,
        annotation[, keep_anno_cols, drop = FALSE],
        by = "metabolite_id",
        suffix = c("", "_annotation")
      )
    }
  }

  diff_table
}

write_differential_heatmap_not_drawn <- function(filename, diff_table, diff_ids, hm_diff, reason) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) return(invisible(NULL))

  summary_df <- data.frame(
    heatmap_type = "differential",
    status = "not_drawn",
    reason = reason,
    n_differential_rows = nrow(diff_table),
    n_unique_differential_metabolites = length(diff_ids),
    stringsAsFactors = FALSE
  )

  unique_ids_df <- data.frame(
    metabolite_id = diff_ids,
    stringsAsFactors = FALSE
  )

  openxlsx::write.xlsx(
    list(
      heatmap_summary = summary_df,
      differential_candidates = diff_table,
      differential_unique_ids = unique_ids_df,
      sample_order_original = hm_diff$sample_info,
      metabolite_annotation = hm_diff$annotation
    ),
    file = filename,
    overwrite = TRUE
  )
}

## ---------------------------------------------------------
## 1. Total ComplexHeatmap: all metabolites x all samples
## ---------------------------------------------------------
if (isTRUE(make_total_heatmap)) {
  tryCatch({
    message("Plotting total ComplexHeatmap: all metabolites x all samples")

    hm <- get_heatmap_matrix_and_annotation(heatmap_use_matrix)
    heatmap_mat <- hm$mat

    heatmap_title <- if (heatmap_use_matrix == "log") "All metabolites heatmap" else "All metabolites heatmap"

    plot_total_heatmap_all_metabolites(
      expr_mat = heatmap_mat,
      sample_info = hm$sample_info,
      annotation = hm$annotation,
      filename_prefix = file.path(dir_heatmap_total, "Total_ComplexHeatmap_all_metabolites_CK_Fg_time_order"),
      width = heatmap_width,
      height = heatmap_height,
      dpi = plot_dpi,
      auto_height = heatmap_auto_height,
      row_height_in = heatmap_row_height_in,
      max_height = heatmap_max_height,
      cluster_rows = heatmap_cluster_rows,
      cluster_columns = heatmap_cluster_columns,
      show_row_names = heatmap_show_metabolite_names,
      show_column_names = heatmap_show_sample_names,
      base_family = base_font_family,
      clip_value = heatmap_clip_value,
      color_quantiles = heatmap_color_quantiles,
      value_label = hm$value_label,
      title = heatmap_title,
      heatmap_type = "total",
      extra_sheets = list(
        total_description = data.frame(
          heatmap_type = "total",
          matrix_type = heatmap_use_matrix,
          description = "All metabolites heatmap matrix and plot-order annotation",
          stringsAsFactors = FALSE
        )
      )
    )
  }, error = function(e) {
    warning("Total ComplexHeatmap failed: ", e$message)
  })
}

## ---------------------------------------------------------
## 1B. Differential metabolite ComplexHeatmap: union of all significant metabolites
## ---------------------------------------------------------
if (isTRUE(make_total_heatmap)) {
  tryCatch({
    message("Plotting differential-metabolite ComplexHeatmap: union of all significant metabolites")

    hm_diff <- get_heatmap_matrix_and_annotation(heatmap_use_matrix)

    diff_table <- make_differential_heatmap_table(
      limma_results = limma_results,
      stat_col = volcano_stat_col,
      pvalue_cutoff = pvalue_cutoff,
      logFC_cutoff = logFC_cutoff,
      annotation = hm_diff$annotation
    )

    diff_ids <- unique(diff_table$metabolite_id)
    diff_ids <- diff_ids[!is.na(diff_ids) & diff_ids != ""]
    diff_ids_in_matrix <- intersect(diff_ids, rownames(hm_diff$mat))

    if (length(diff_ids_in_matrix) >= 2) {
      heatmap_mat_diff <- hm_diff$mat[diff_ids_in_matrix, , drop = FALSE]

      plot_total_heatmap_all_metabolites(
        expr_mat = heatmap_mat_diff,
        sample_info = hm_diff$sample_info,
        annotation = hm_diff$annotation,
        filename_prefix = file.path(dir_heatmap_differential, "Differential_ComplexHeatmap_all_DE_metabolites_CK_Fg_time_order"),
        width = heatmap_width,
        height = heatmap_height,
        dpi = plot_dpi,
        auto_height = TRUE,
        row_height_in = min(heatmap_row_height_in, 0.09),
        max_height = heatmap_max_height,
        cluster_rows = heatmap_cluster_rows,
        cluster_columns = heatmap_cluster_columns,
        show_row_names = FALSE,
        show_column_names = FALSE,
        base_family = base_font_family,
        clip_value = heatmap_clip_value,
        color_quantiles = heatmap_color_quantiles,
        value_label = hm_diff$value_label,
        title = "Differential_ComplexHeatmap_all_DE_metabolites_CK_Fg_time_order",
        heatmap_type = "differential",
        extra_sheets = list(
          differential_candidates = diff_table,
          differential_unique_ids = data.frame(metabolite_id = diff_ids, stringsAsFactors = FALSE),
          differential_ids_in_matrix = data.frame(metabolite_id = diff_ids_in_matrix, stringsAsFactors = FALSE)
        )
      )
    } else {
      reason <- paste0(
        "Too few differential metabolites to draw heatmap. Unique significant IDs = ",
        length(diff_ids),
        "; IDs present in selected heatmap matrix = ",
        length(diff_ids_in_matrix),
        ". Threshold: ",
        volcano_stat_col,
        " < ",
        pvalue_cutoff,
        " and |logFC| >= ",
        logFC_cutoff,
        "."
      )
      message(reason)

      write_differential_heatmap_not_drawn(
        filename = file.path(dir_heatmap_differential, "Differential_heatmap_NOT_DRAWN_matrix_and_annotation.xlsx"),
        diff_table = diff_table,
        diff_ids = diff_ids,
        hm_diff = hm_diff,
        reason = reason
      )
    }
  }, error = function(e) {
    warning("Differential-metabolite ComplexHeatmap failed: ", e$message)

    hm_diff <- tryCatch(get_heatmap_matrix_and_annotation(heatmap_use_matrix), error = function(e2) NULL)
    if (!is.null(hm_diff)) {
      write_differential_heatmap_not_drawn(
        filename = file.path(dir_heatmap_differential, "Differential_heatmap_FAILED_matrix_and_annotation.xlsx"),
        diff_table = data.frame(),
        diff_ids = character(0),
        hm_diff = hm_diff,
        reason = e$message
      )
    }
  })
}

## ---------------------------------------------------------
## 2. Volcano (VIP as point size) + focus-class boxplot
## ---------------------------------------------------------
for (i in seq_len(nrow(comparisons))) {
  comp <- comparisons[i, ]
  nm <- comp$comparison

  if (!nm %in% names(limma_results)) next

  message("Plotting comparison figures: ", nm)
  res <- limma_results[[nm]]
  vip_df <- NULL
  if (length(opls_results) > 0 && nm %in% names(opls_results)) {
    vip_df <- opls_results[[nm]]$vip
  }

  ## Volcano with VIP-based point size
  p_vol <- tryCatch({
    plot_volcano(
      res = res,
      vip_df = vip_df,
      logFC_cutoff = logFC_cutoff,
      pvalue_cutoff = pvalue_cutoff,
      stat_col = volcano_stat_col,
      y_label = volcano_y_label,
      focus_classes = focus_classes,
      focus_bx_patterns = focus_bx_patterns,
      label_focus_only = volcano_label_focus_only,
      label_sig_focus_only = volcano_label_sig_focus_only,
      label_all_focus = volcano_label_all_focus,
      label_top = volcano_label_top_n,
      title = format_comparison_title(comp),
      base_size = base_font_size,
      base_family = base_font_family
    )
  }, error = function(e) {
    warning("Volcano failed for ", nm, ": ", e$message)
    NULL
  })

  if (!is.null(p_vol)) {
    tryCatch({
      save_ggplot_multi(
        p = p_vol,
        filename_prefix = file.path(dir_volcano, paste0(safe_filename(nm), "_volcano")),
        width = volcano_width,
        height = volcano_height,
        dpi = plot_dpi
      )
    }, error = function(e) {
      warning("Saving volcano failed for ", nm, ": ", e$message)
    })
  }

  ## Boxplot for significant focus-class metabolites
  p_box <- tryCatch({
    plot_boxplot_top(
      expr_log = data02$expr_log,
      res = res,
      sample_info = data02$sample_info,
      comparison_row = comp,
      annotation = data02$annotation,
      top_n = boxplot_top_n,
      focus_classes = focus_classes,
      focus_bx_patterns = focus_bx_patterns,
      stat_col = boxplot_sig_stat_col,
      pvalue_cutoff = boxplot_pvalue_cutoff,
      logFC_cutoff = logFC_cutoff,
      focus_only = boxplot_focus_only,
      show_all_sig_focus = boxplot_show_all_sig_focus,
      max_facets = boxplot_max_facets,
      facet_ncol = boxplot_ncol,
      title = format_comparison_title(comp),
      base_size = base_font_size,
      base_family = base_font_family
    )
  }, error = function(e) {
    warning("Boxplot failed for ", nm, ": ", e$message)
    NULL
  })

  if (!is.null(p_box)) {
    n_facets <- attr(p_box, "n_facets")
    ncol_fac <- attr(p_box, "facet_ncol")
    if (is.null(n_facets) || !is.finite(n_facets)) n_facets <- boxplot_top_n
    if (is.null(ncol_fac) || !is.finite(ncol_fac)) ncol_fac <- boxplot_ncol
    nrow_fac <- ceiling(n_facets / ncol_fac)
    this_height <- max(boxplot_min_height, nrow_fac * boxplot_row_height)

    tryCatch({
      save_ggplot_multi(
        p = p_box,
        filename_prefix = file.path(dir_boxplot, paste0(safe_filename(nm), "_boxplot_focus_sig_metabolites")),
        width = boxplot_width,
        height = this_height,
        dpi = plot_dpi
      )
    }, error = function(e) {
      warning("Saving boxplot failed for ", nm, ": ", e$message)
    })
  }
}


## ---------------------------------------------------------
## Heatmap output index
## ---------------------------------------------------------
if (requireNamespace("openxlsx", quietly = TRUE)) {
  total_files <- list.files(dir_heatmap_total, full.names = FALSE)
  differential_files <- list.files(dir_heatmap_differential, full.names = FALSE)

  openxlsx::write.xlsx(
    list(
      total_outputs = data.frame(file = total_files, stringsAsFactors = FALSE),
      differential_outputs = data.frame(file = differential_files, stringsAsFactors = FALSE),
      output_folders = data.frame(
        heatmap_type = c("total", "differential"),
        folder = c(dir_heatmap_total, dir_heatmap_differential),
        stringsAsFactors = FALSE
      )
    ),
    file = file.path(dir_heatmap, "Heatmap_output_index.xlsx"),
    overwrite = TRUE
  )
}


message("07_Plots finished.")
