## =========================================================
## plot_heatmap.R
## ComplexHeatmap publication-level heatmaps
## v11: total matrix heatmap | class-grouped like v7 | square-cell style
## - rows: all metabolites
## - columns: all samples
## - columns sorted by CK/Fg split, then time order, genotype, replicate
## - rows annotated and split by Class
## =========================================================

order_sample_info_for_total_heatmap <- function(
    sample_info,
    treatment_levels = c("CK", "Fg"),
    time_levels = c("0h", "24h", "48h")
) {
  sample_info <- as.data.frame(sample_info, check.names = FALSE)

  if (!"time_num" %in% colnames(sample_info)) {
    sample_info$time_num <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", sample_info$time)))
  }

  sample_info$genotype <- as.character(sample_info$genotype)
  genotype_levels <- unique(sample_info$genotype)
  sample_info$genotype <- factor(sample_info$genotype, levels = genotype_levels)

  sample_info$treatment <- factor(as.character(sample_info$treatment), levels = treatment_levels)
  sample_info$time <- as.character(sample_info$time)
  sample_info$rep <- suppressWarnings(as.numeric(as.character(sample_info$rep)))

  ## New order for clearer comparison within each genotype:
  ## genotype -> CK time order -> Fg time order -> replicate
  ## Example: M03 CK 0dpi, CK 12dpi, CK 24dpi, Fg 12dpi, Fg 24dpi
  sample_info$treatment_order <- ifelse(as.character(sample_info$treatment) == "CK", 0, 1)
  sample_info <- sample_info[order(
    sample_info$genotype,
    sample_info$treatment_order,
    sample_info$time_num,
    sample_info$rep,
    na.last = TRUE
  ), , drop = FALSE]

  sample_info
}

get_metabolite_labels_and_class <- function(
    metabolite_ids,
    annotation = NULL,
    class_col = "Class"
) {
  row_labels <- metabolite_ids
  row_class <- rep("Unknown", length(metabolite_ids))
  source_features <- rep(NA_character_, length(metabolite_ids))

  if (!is.null(annotation)) {
    annotation <- as.data.frame(annotation, check.names = FALSE)

    if ("metabolite_id" %in% colnames(annotation)) {
      anno_tmp <- annotation[match(metabolite_ids, annotation$metabolite_id), , drop = FALSE]

      if ("Compound name" %in% colnames(anno_tmp)) {
        row_labels <- ifelse(
          is.na(anno_tmp[["Compound name"]]) | anno_tmp[["Compound name"]] == "",
          metabolite_ids,
          anno_tmp[["Compound name"]]
        )
      }

      if (class_col %in% colnames(anno_tmp)) {
        row_class <- as.character(anno_tmp[[class_col]])
        row_class[is.na(row_class) | row_class == ""] <- "Unknown"
        ## 如果重复分子合并后 Class 形如 A;B，这里取第一个用于热图分组；完整 Class 仍保留在 annotation 表里。
        row_class <- sub(";.*$", "", row_class)
      }

      if ("source_features" %in% colnames(anno_tmp)) {
        source_features <- as.character(anno_tmp$source_features)
      }
    }
  }

  data.frame(
    metabolite_id = metabolite_ids,
    metabolite_label = make.unique(as.character(row_labels)),
    Class = row_class,
    source_features = source_features,
    stringsAsFactors = FALSE
  )
}

make_named_palette <- function(values, palette = "Set3") {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0) return(character(0))

  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, min(length(values), max_n))), palette)
    if (length(values) > length(base_cols)) {
      cols <- grDevices::colorRampPalette(base_cols)(length(values))
    } else {
      cols <- base_cols[seq_len(length(values))]
    }
  } else {
    cols <- grDevices::rainbow(length(values))
  }

  stats::setNames(cols, values)
}

prepare_total_heatmap_matrix <- function(
    expr_mat,
    sample_info,
    annotation = NULL,
    clip_value = 3,
    treatment_levels = c("CK", "Fg"),
    time_levels = c("0h", "24h", "48h")
) {
  expr_mat <- as_numeric_matrix(expr_mat)
  expr_mat[!is.finite(expr_mat)] <- NA_real_

  sample_info_ordered <- order_sample_info_for_total_heatmap(
    sample_info = sample_info,
    treatment_levels = treatment_levels,
    time_levels = time_levels
  )

  selected_samples <- intersect(sample_info_ordered$sample, colnames(expr_mat))
  sample_info_ordered <- sample_info_ordered[match(selected_samples, sample_info_ordered$sample), , drop = FALSE]

  if (length(selected_samples) < 2) {
    stop("Too few samples for total heatmap.")
  }

  mat <- expr_mat[, selected_samples, drop = FALSE]
  mat <- as_numeric_matrix(mat)

  ## 对残余 NA 做行均值填补；如果整行全 NA，则填 0。
  for (i in seq_len(nrow(mat))) {
    if (any(is.na(mat[i, ]))) {
      row_mean <- mean(mat[i, ], na.rm = TRUE)
      if (!is.finite(row_mean)) row_mean <- 0
      mat[i, is.na(mat[i, ])] <- row_mean
    }
  }

  ## 不删除零方差行：用户要求展示所有代谢物。
  ## 但如果仍有非有限值，统一置 0，避免 ComplexHeatmap 报错。
  mat[!is.finite(mat)] <- 0

  if (!is.null(clip_value) && is.finite(clip_value) && clip_value > 0) {
    mat[mat > clip_value] <- clip_value
    mat[mat < -clip_value] <- -clip_value
  }

  row_anno <- get_metabolite_labels_and_class(
    metabolite_ids = rownames(mat),
    annotation = annotation,
    class_col = "Class"
  )

  rownames(mat) <- row_anno$metabolite_label

  sample_info_ordered$treatment <- factor(as.character(sample_info_ordered$treatment), levels = treatment_levels)
  sample_info_ordered$time <- factor(as.character(sample_info_ordered$time), levels = time_levels)
  sample_info_ordered$column_split <- factor(as.character(sample_info_ordered$treatment), levels = treatment_levels)

  list(
    mat = mat,
    col_info = sample_info_ordered,
    row_anno = row_anno
  )
}

save_complex_heatmap_multi <- function(
    ht,
    filename_prefix,
    width = 16,
    height = 10,
    dpi = 600
) {
  dir.create(dirname(filename_prefix), recursive = TRUE, showWarnings = FALSE)

  grDevices::cairo_pdf(paste0(filename_prefix, ".pdf"), width = width, height = height)
  ht_drawn <- ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = TRUE
  )
  grDevices::dev.off()

  grDevices::png(
    filename = paste0(filename_prefix, ".png"),
    width = width,
    height = height,
    units = "in",
    res = dpi,
    type = "cairo",
    bg = "white"
  )
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = TRUE
  )
  grDevices::dev.off()

  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(paste0(filename_prefix, ".svg"), width = width, height = height, bg = "white")
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legend = TRUE
    )
    grDevices::dev.off()
  }

  invisible(ht_drawn)
}

build_total_complex_heatmap <- function(
    mat,
    col_info,
    row_anno,
    title = "All metabolites heatmap",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    base_family = "Arial",
    clip_value = NA,
    color_quantiles = c(0.05, 0.50, 0.95),
    value_label = "Peak area"
) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
    library(scales)
  })

  genotype_levels <- unique(as.character(col_info$genotype))
  class_levels <- unique(as.character(row_anno$Class))
  class_levels <- class_levels[!is.na(class_levels)]

  genotype_cols <- stats::setNames(
    scales::hue_pal(l = 65, c = 90)(length(genotype_levels)),
    genotype_levels
  )

  class_cols <- make_named_palette(class_levels, palette = "Set3")

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    Genotype = as.character(col_info$genotype),
    Time = as.character(col_info$time),
    Treatment = as.character(col_info$treatment),
    col = list(
      Genotype = genotype_cols,
      Time = metabo_palette("time"),
      Treatment = metabo_palette("treatment")
    ),
    annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold", fontfamily = base_family),
    annotation_name_side = "left",
    annotation_legend_param = list(
      Genotype = list(title_gp = grid::gpar(fontface = "bold", fontfamily = base_family)),
      Time = list(title_gp = grid::gpar(fontface = "bold", fontfamily = base_family)),
      Treatment = list(title_gp = grid::gpar(fontface = "bold", fontfamily = base_family))
    )
  )

  left_anno <- ComplexHeatmap::rowAnnotation(
    Class = as.character(row_anno$Class),
    col = list(Class = class_cols),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Class = list(title_gp = grid::gpar(fontface = "bold", fontfamily = base_family))
    )
  )

  row_split <- factor(as.character(row_anno$Class), levels = class_levels)
  column_split <- factor(as.character(col_info$genotype), levels = genotype_levels)

  row_font_size <- if (nrow(mat) <= 60) 7 else if (nrow(mat) <= 150) 5.5 else 4
  col_font_size <- if (ncol(mat) <= 60) 6 else 4.5

  ComplexHeatmap::Heatmap(
    mat,
    name = value_label,
    col = {
      finite_vals <- as.numeric(mat[is.finite(mat)])
      finite_vals <- finite_vals[!is.na(finite_vals)]

      if (length(finite_vals) == 0) {
        brks <- c(0, 0.5, 1)
      } else if (!is.null(clip_value) && is.finite(clip_value) && clip_value > 0) {
        brks <- c(-clip_value, 0, clip_value)
      } else {
        qs <- suppressWarnings(stats::quantile(
          finite_vals,
          probs = color_quantiles,
          na.rm = TRUE,
          names = FALSE,
          type = 7
        ))
        qs <- as.numeric(qs)
        qs <- qs[is.finite(qs)]

        if (length(qs) < 3 || length(unique(qs)) < 3) {
          rng <- range(finite_vals, na.rm = TRUE)
          if (!all(is.finite(rng)) || rng[1] == rng[2]) {
            brks <- c(rng[1] - 1, rng[1], rng[1] + 1)
          } else {
            brks <- c(rng[1], mean(rng), rng[2])
          }
        } else {
          brks <- qs
        }
      }

      circlize::colorRamp2(
        brks,
        c("#2166AC", "#F7F7F7", "#B2182B")
      )
    },
    top_annotation = top_anno,
    left_annotation = left_anno,
    row_split = row_split,
    column_split = column_split,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_gp = grid::gpar(fontsize = row_font_size, fontfamily = base_family),
    column_names_gp = grid::gpar(fontsize = col_font_size, fontfamily = base_family),
    row_title = NULL,
    row_title_gp = grid::gpar(fontsize = 0.01, col = NA, fontfamily = base_family),
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold", fontfamily = base_family),
    column_title = title,
    column_gap = grid::unit(2.5, "mm"),
    row_gap = grid::unit(1.5, "mm"),
    border = TRUE,
    rect_gp = grid::gpar(col = "white", lwd = 0.08),
    width = grid::unit(ncol(mat) * 2.2, "mm"),
    height = grid::unit(nrow(mat) * 2.2, "mm"),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontface = "bold", fontfamily = base_family),
      labels_gp = grid::gpar(fontfamily = base_family)
    )
  )
}

plot_total_heatmap_all_metabolites <- function(
    expr_mat,
    sample_info,
    annotation = NULL,
    filename_prefix,
    width = 16,
    height = 10,
    dpi = 600,
    auto_height = TRUE,
    row_height_in = 0.11,
    max_height = 60,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    base_family = "Arial",
    clip_value = NA,
    color_quantiles = c(0.05, 0.50, 0.95),
    value_label = "Peak area",
    title = "All metabolites heatmap",
    heatmap_type = "total",
    extra_sheets = NULL
) {
  prep <- prepare_total_heatmap_matrix(
    expr_mat = expr_mat,
    sample_info = sample_info,
    annotation = annotation,
    clip_value = clip_value
  )

  if (isTRUE(auto_height)) {
    height <- max(height, min(max_height, 4 + nrow(prep$mat) * row_height_in))
  }

  message("  [Total heatmap] metabolites: ", nrow(prep$mat))
  message("  [Total heatmap] samples: ", ncol(prep$mat))
  message("  [Total heatmap] column order: genotype -> CK time order -> Fg time order -> replicate")
  message("  [Total heatmap] matrix values: no Z-score; color follows value magnitude")
  message("  [Total heatmap] row annotation: Class legend only")

  ht <- build_total_complex_heatmap(
    mat = prep$mat,
    col_info = prep$col_info,
    row_anno = prep$row_anno,
    title = title,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    base_family = base_family,
    clip_value = clip_value,
    color_quantiles = color_quantiles,
    value_label = value_label
  )

  ht_drawn <- save_complex_heatmap_multi(
    ht = ht,
    filename_prefix = filename_prefix,
    width = width,
    height = height,
    dpi = dpi
  )

  ## 同时导出热图实际使用的矩阵、样本顺序、代谢物注释，以及 ComplexHeatmap 实际绘图顺序。
  if (requireNamespace("openxlsx", quietly = TRUE)) {

    row_order_raw <- tryCatch(ComplexHeatmap::row_order(ht_drawn), error = function(e) NULL)
    col_order_raw <- tryCatch(ComplexHeatmap::column_order(ht_drawn), error = function(e) NULL)

    if (is.null(row_order_raw)) {
      row_order_vec <- seq_len(nrow(prep$mat))
      row_slice <- rep("all", length(row_order_vec))
    } else if (is.list(row_order_raw)) {
      row_order_vec <- unlist(row_order_raw, use.names = FALSE)
      row_slice <- rep(names(row_order_raw), lengths(row_order_raw))
      if (is.null(row_slice) || length(row_slice) == 0) {
        row_slice <- rep(seq_along(row_order_raw), lengths(row_order_raw))
      }
    } else {
      row_order_vec <- as.integer(row_order_raw)
      row_slice <- rep("all", length(row_order_vec))
    }

    if (is.null(col_order_raw)) {
      col_order_vec <- seq_len(ncol(prep$mat))
      col_slice <- rep("all", length(col_order_vec))
    } else if (is.list(col_order_raw)) {
      col_order_vec <- unlist(col_order_raw, use.names = FALSE)
      col_slice <- rep(names(col_order_raw), lengths(col_order_raw))
      if (is.null(col_slice) || length(col_slice) == 0) {
        col_slice <- rep(seq_along(col_order_raw), lengths(col_order_raw))
      }
    } else {
      col_order_vec <- as.integer(col_order_raw)
      col_slice <- rep("all", length(col_order_vec))
    }

    row_order_table <- data.frame(
      heatmap_row_position = seq_along(row_order_vec),
      row_slice = row_slice,
      matrix_row_index = row_order_vec,
      metabolite_label = rownames(prep$mat)[row_order_vec],
      metabolite_id = prep$row_anno$metabolite_id[row_order_vec],
      Class = prep$row_anno$Class[row_order_vec],
      source_features = prep$row_anno$source_features[row_order_vec],
      stringsAsFactors = FALSE
    )

    column_order_table <- data.frame(
      heatmap_column_position = seq_along(col_order_vec),
      column_slice = col_slice,
      matrix_column_index = col_order_vec,
      sample = colnames(prep$mat)[col_order_vec],
      genotype = prep$col_info$genotype[col_order_vec],
      treatment = prep$col_info$treatment[col_order_vec],
      time = prep$col_info$time[col_order_vec],
      rep = prep$col_info$rep[col_order_vec],
      group = prep$col_info$group[col_order_vec],
      stringsAsFactors = FALSE
    )

    heatmap_matrix_plot_order <- data.frame(
      metabolite_label = rownames(prep$mat)[row_order_vec],
      prep$mat[row_order_vec, col_order_vec, drop = FALSE],
      check.names = FALSE
    )

    heatmap_summary <- data.frame(
      heatmap_type = heatmap_type,
      title = title,
      value_label = value_label,
      n_metabolites = nrow(prep$mat),
      n_samples = ncol(prep$mat),
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      stringsAsFactors = FALSE
    )

    row_order_by_class <- row_order_table[order(row_order_table$row_slice, row_order_table$heatmap_row_position), , drop = FALSE]
    column_order_by_group <- column_order_table[order(column_order_table$column_slice, column_order_table$heatmap_column_position), , drop = FALSE]

    ## A matrix with explicit row metadata, in the exact plotted row/column order.
    heatmap_matrix_plot_order_with_info <- data.frame(
      heatmap_row_position = seq_along(row_order_vec),
      row_slice = row_slice,
      metabolite_label = rownames(prep$mat)[row_order_vec],
      metabolite_id = prep$row_anno$metabolite_id[row_order_vec],
      Class = prep$row_anno$Class[row_order_vec],
      source_features = prep$row_anno$source_features[row_order_vec],
      prep$mat[row_order_vec, col_order_vec, drop = FALSE],
      check.names = FALSE
    )

    ## Long-format table: one row = one heatmap cell.
    ## This is the most explicit mapping from color block to metabolite/sample.
    heatmap_cell_long <- do.call(
      rbind,
      lapply(seq_along(row_order_vec), function(i) {
        ridx <- row_order_vec[i]
        data.frame(
          heatmap_row_position = i,
          row_slice = row_slice[i],
          metabolite_label = rownames(prep$mat)[ridx],
          metabolite_id = prep$row_anno$metabolite_id[ridx],
          Class = prep$row_anno$Class[ridx],
          source_features = prep$row_anno$source_features[ridx],
          heatmap_column_position = seq_along(col_order_vec),
          column_slice = col_slice,
          sample = colnames(prep$mat)[col_order_vec],
          genotype = prep$col_info$genotype[col_order_vec],
          treatment = prep$col_info$treatment[col_order_vec],
          time = prep$col_info$time[col_order_vec],
          rep = prep$col_info$rep[col_order_vec],
          group = prep$col_info$group[col_order_vec],
          value = as.numeric(prep$mat[ridx, col_order_vec]),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      })
    )

    make_safe_sheet_name <- function(prefix, value, used_names) {
      value <- as.character(value)
      value[is.na(value) | value == ""] <- "Unknown"
      nm <- paste0(prefix, "_", value)
      nm <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", nm)
      nm <- gsub("\\s+", "_", nm)
      nm <- gsub("_+", "_", nm)
      nm <- substr(nm, 1, 31)

      base_nm <- nm
      k <- 1
      while (nm %in% used_names) {
        suffix <- paste0("_", k)
        nm <- paste0(substr(base_nm, 1, 31 - nchar(suffix)), suffix)
        k <- k + 1
      }
      nm
    }

    grouped_sheets <- list()
    used_sheet_names <- names(grouped_sheets)

    ## Split matrix by column group exactly as shown in the heatmap column slices.
    for (grp in unique(column_order_table$column_slice)) {
      these_cols <- column_order_table$matrix_column_index[column_order_table$column_slice == grp]
      if (length(these_cols) < 1) next

      sheet_nm <- make_safe_sheet_name("colgrp", grp, used_sheet_names)
      used_sheet_names <- c(used_sheet_names, sheet_nm)

      grouped_sheets[[sheet_nm]] <- data.frame(
        heatmap_row_position = seq_along(row_order_vec),
        row_slice = row_slice,
        metabolite_label = rownames(prep$mat)[row_order_vec],
        metabolite_id = prep$row_anno$metabolite_id[row_order_vec],
        Class = prep$row_anno$Class[row_order_vec],
        source_features = prep$row_anno$source_features[row_order_vec],
        prep$mat[row_order_vec, these_cols, drop = FALSE],
        check.names = FALSE
      )
    }

    ## Split matrix by row group exactly as shown in the heatmap row slices.
    for (grp in unique(row_order_table$row_slice)) {
      these_rows <- row_order_table$matrix_row_index[row_order_table$row_slice == grp]
      if (length(these_rows) < 1) next

      sheet_nm <- make_safe_sheet_name("rowgrp", grp, used_sheet_names)
      used_sheet_names <- c(used_sheet_names, sheet_nm)

      grouped_sheets[[sheet_nm]] <- data.frame(
        heatmap_row_position = row_order_table$heatmap_row_position[row_order_table$row_slice == grp],
        row_slice = grp,
        metabolite_label = rownames(prep$mat)[these_rows],
        metabolite_id = prep$row_anno$metabolite_id[these_rows],
        Class = prep$row_anno$Class[these_rows],
        source_features = prep$row_anno$source_features[these_rows],
        prep$mat[these_rows, col_order_vec, drop = FALSE],
        check.names = FALSE
      )
    }

    group_sheet_index <- data.frame(
      sheet_name = names(grouped_sheets),
      split_type = ifelse(grepl("^colgrp_", names(grouped_sheets)), "column_group", "row_group"),
      group_name = sub("^(colgrp|rowgrp)_", "", names(grouped_sheets)),
      stringsAsFactors = FALSE
    )

    xlsx_sheets <- c(
      list(
        heatmap_summary = heatmap_summary,
        group_sheet_index = group_sheet_index,
        heatmap_matrix_original_order = data.frame(
          metabolite_label = rownames(prep$mat),
          prep$mat,
          check.names = FALSE
        ),
        heatmap_matrix_plot_order = heatmap_matrix_plot_order,
        matrix_plot_order_with_info = heatmap_matrix_plot_order_with_info,
        heatmap_cell_long = heatmap_cell_long,
        heatmap_row_order = row_order_table,
        heatmap_column_order = column_order_table,
        row_order_by_Class = row_order_by_class,
        column_order_by_group = column_order_by_group,
        sample_order_original = prep$col_info,
        metabolite_annotation = prep$row_anno
      ),
      grouped_sheets
    )

    if (!is.null(extra_sheets) && length(extra_sheets) > 0) {
      xlsx_sheets <- c(
        list(heatmap_summary = heatmap_summary),
        extra_sheets,
        xlsx_sheets[names(xlsx_sheets) != "heatmap_summary"]
      )
    }

    openxlsx::write.xlsx(
      xlsx_sheets,
      file = paste0(filename_prefix, "_matrix_and_annotation.xlsx"),
      overwrite = TRUE
    )
  }

  invisible(prep)
}

## ---------------------------------------------------------
## Backward-compatible wrappers from v3
## These are kept so old scripts do not break, but v4 uses
## plot_total_heatmap_all_metabolites() as the main heatmap.
## ---------------------------------------------------------
select_heatmap_ids_from_result <- function(
    res,
    top_n = 50,
    logFC_cutoff = log2(1.5),
    adjP_cutoff = 0.05
) {
  suppressPackageStartupMessages({ library(dplyr) })
  res <- as.data.frame(res, check.names = FALSE)
  res$logFC <- suppressWarnings(as.numeric(res$logFC))
  res$adj.P.Val <- suppressWarnings(as.numeric(res$adj.P.Val))

  sig_ids <- res %>%
    filter(
      is.finite(adj.P.Val),
      is.finite(logFC),
      adj.P.Val < adjP_cutoff,
      abs(logFC) >= logFC_cutoff
    ) %>%
    arrange(adj.P.Val) %>%
    pull(metabolite_id)

  if (length(sig_ids) < 2) {
    res2 <- res %>% filter(is.finite(adj.P.Val)) %>% arrange(adj.P.Val)
    n_top <- min(top_n, nrow(res2))
    if (n_top < 2) return(character(0))
    sig_ids <- res2$metabolite_id[seq_len(n_top)]
  } else {
    sig_ids <- sig_ids[seq_len(min(top_n, length(sig_ids)))]
  }

  unique(sig_ids)
}

plot_heatmap_comparison <- function(
    expr_scaled,
    res,
    sample_info,
    comparison_row,
    annotation = NULL,
    top_n = 50,
    logFC_cutoff = log2(1.5),
    adjP_cutoff = 0.05,
    filename_prefix,
    show_all_groups = TRUE,
    width = 16,
    height = 10,
    dpi = 600,
    show_row_names_limit = Inf,
    base_family = "Arial"
) {
  selected_ids <- select_heatmap_ids_from_result(
    res = res,
    top_n = top_n,
    logFC_cutoff = logFC_cutoff,
    adjP_cutoff = adjP_cutoff
  )

  if (length(selected_ids) < 2) {
    warning("Too few metabolites for comparison heatmap: ", comparison_row$comparison)
    return(NULL)
  }

  expr_scaled <- as_numeric_matrix(expr_scaled)
  expr_sub <- expr_scaled[selected_ids, , drop = FALSE]

  plot_total_heatmap_all_metabolites(
    expr_mat = expr_sub,
    sample_info = sample_info,
    annotation = annotation,
    filename_prefix = filename_prefix,
    width = width,
    height = height,
    dpi = dpi,
    auto_height = TRUE,
    row_height_in = 0.11,
    max_height = 60,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = nrow(expr_sub) <= show_row_names_limit,
    show_column_names = TRUE,
    base_family = base_family,
    clip_value = 3,
    title = paste0(comparison_row$comparison, " | selected metabolites | all samples")
  )
}

plot_heatmap_global_all_groups <- function(
    expr_scaled,
    limma_results,
    sample_info,
    annotation = NULL,
    top_n = 100,
    logFC_cutoff = log2(1.5),
    adjP_cutoff = 0.05,
    filename_prefix,
    width = 16,
    height = 10,
    dpi = 600,
    show_row_names_limit = Inf,
    base_family = "Arial"
) {
  plot_total_heatmap_all_metabolites(
    expr_mat = expr_scaled,
    sample_info = sample_info,
    annotation = annotation,
    filename_prefix = filename_prefix,
    width = width,
    height = height,
    dpi = dpi,
    auto_height = TRUE,
    row_height_in = 0.11,
    max_height = 60,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    base_family = base_family,
    clip_value = 3,
    title = "All metabolites heatmap"
  )
}
