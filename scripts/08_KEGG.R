## =========================================================
## 08_KEGG.R
## KEGG compound enrichment
## v10:
##   - completely safe skip when kegg_annotation.xlsx is absent/malformed
##   - no top-level return() statements
##   - no pipeline interruption when KEGG annotation is not ready
## =========================================================

source("config/config.R")
source("functions/utils/common_utils.R")
source("functions/pathway/kegg_enrichment.R")
source("functions/plots/save_plot.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

data02 <- readRDS(file.path(dir_rds, "02_preprocessed.rds"))
limma_results <- readRDS(file.path(dir_rds, "05_limma_results.rds"))

annotation <- as.data.frame(data02$annotation, check.names = FALSE)

## ---------------------------------------------------------
## Helper: write KEGG annotation template
## ---------------------------------------------------------
write_kegg_template <- function(annotation, outfile) {
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

  template <- as.data.frame(annotation, check.names = FALSE)

  if (!"Compound name" %in% colnames(template)) {
    if ("metabolite_id" %in% colnames(template)) {
      template[["Compound name"]] <- template$metabolite_id
    } else {
      template[["Compound name"]] <- NA_character_
    }
  }

  if (!"Class" %in% colnames(template)) {
    template[["Class"]] <- NA_character_
  }

  keep_cols <- intersect(
    c("metabolite_id", "Compound name", "Class", "source_features"),
    colnames(template)
  )

  template <- template[, keep_cols, drop = FALSE]
  template[[kegg_id_col]] <- NA_character_
  template$note <- "Fill KEGG Compound ID such as C00065, C00148, C00106. Leave blank if unknown."

  openxlsx::write.xlsx(
    list(KEGG_annotation_template = template),
    file = outfile,
    overwrite = TRUE
  )

  invisible(outfile)
}

write_kegg_skip_output <- function(reason, extra = NULL) {
  msg <- data.frame(
    status = "skipped",
    reason = reason,
    required_columns = paste(c("Compound name", kegg_id_col), collapse = ", "),
    kegg_annotation_file = kegg_annotation_file,
    stringsAsFactors = FALSE
  )

  if (!is.null(extra)) {
    extra <- as.data.frame(extra, check.names = FALSE)
    out <- list(message = msg, extra = extra)
  } else {
    out <- list(message = msg)
  }

  write_xlsx_list(
    out,
    file.path(dir_tables, "08_KEGG_output.xlsx")
  )

  saveRDS(list(), file.path(dir_rds, "08_kegg_results.rds"))
  invisible(NULL)
}

template_file <- file.path(dir_tables, "08_KEGG_annotation_template.xlsx")

## ---------------------------------------------------------
## 1. Decide whether KEGG can run
## ---------------------------------------------------------
run_kegg_now <- TRUE
skip_reason <- NULL
provided_preview <- NULL

if (!file.exists(kegg_annotation_file)) {
  run_kegg_now <- FALSE
  skip_reason <- "kegg_annotation.xlsx not found"

  message("KEGG annotation file not found:")
  message("  ", kegg_annotation_file)
  message("KEGG enrichment skipped.")
  message("A template has been written to:")
  message("  ", template_file)

  write_kegg_template(annotation, template_file)
}

kegg_anno <- NULL

if (isTRUE(run_kegg_now)) {
  kegg_anno <- tryCatch(
    readxl::read_excel(kegg_annotation_file) |> as.data.frame(check.names = FALSE),
    error = function(e) {
      run_kegg_now <<- FALSE
      skip_reason <<- paste0("kegg_annotation.xlsx cannot be read: ", e$message)
      NULL
    }
  )
}

if (isTRUE(run_kegg_now)) {
  names(kegg_anno) <- trimws(names(kegg_anno))
  required_cols <- c("Compound name", kegg_id_col)

  if (!all(required_cols %in% colnames(kegg_anno))) {
    run_kegg_now <- FALSE
    skip_reason <- "kegg_annotation.xlsx lacks required columns"
    provided_preview <- head(kegg_anno, 50)

    message("KEGG annotation file exists but lacks required columns.")
    message("Detected columns:")
    message("  ", paste(colnames(kegg_anno), collapse = ", "))
    message("Required columns:")
    message("  ", paste(required_cols, collapse = ", "))
    message("KEGG enrichment skipped.")
    message("A template has been written to:")
    message("  ", template_file)

    write_kegg_template(annotation, template_file)
  }
}

if (!isTRUE(run_kegg_now)) {
  if (is.null(skip_reason)) skip_reason <- "KEGG annotation is not ready"
  write_kegg_skip_output(skip_reason, provided_preview)
  message("08_KEGG skipped.")
} else {

  ## -------------------------------------------------------
  ## 2. Merge KEGG annotation
  ## -------------------------------------------------------
  kegg_anno[["Compound name"]] <- trimws(as.character(kegg_anno[["Compound name"]]))
  kegg_anno[[kegg_id_col]] <- trimws(as.character(kegg_anno[[kegg_id_col]]))
  kegg_anno[[kegg_id_col]][kegg_anno[[kegg_id_col]] == ""] <- NA_character_
  kegg_anno[[kegg_id_col]] <- gsub("^cpd:", "", kegg_anno[[kegg_id_col]])

  kegg_anno <- kegg_anno[
    !is.na(kegg_anno[["Compound name"]]) &
      kegg_anno[["Compound name"]] != "",
    ,
    drop = FALSE
  ]

  kegg_anno <- kegg_anno[!duplicated(kegg_anno[["Compound name"]]), , drop = FALSE]

  annotation2 <- annotation

  if (!"Compound name" %in% colnames(annotation2)) {
    run_kegg_now <- FALSE
    skip_reason <- "internal annotation table lacks Compound name"
  } else {
    annotation2 <- dplyr::left_join(
      annotation2,
      kegg_anno[, c("Compound name", kegg_id_col), drop = FALSE],
      by = "Compound name"
    )

    n_with_kegg <- sum(!is.na(annotation2[[kegg_id_col]]) & annotation2[[kegg_id_col]] != "")
    message("Metabolites with KEGG ID: ", n_with_kegg, " / ", nrow(annotation2))

    if (n_with_kegg < 3) {
      run_kegg_now <- FALSE
      skip_reason <- paste0("too few metabolites with KEGG ID: ", n_with_kegg)
    }
  }

  if (!isTRUE(run_kegg_now)) {
    write_xlsx_list(
      list(
        message = data.frame(
          status = "skipped",
          reason = skip_reason,
          stringsAsFactors = FALSE
        ),
        merged_annotation = if (exists("annotation2")) annotation2 else annotation
      ),
      file.path(dir_tables, "08_KEGG_output.xlsx")
    )
    saveRDS(list(), file.path(dir_rds, "08_kegg_results.rds"))
    message("08_KEGG skipped.")
  } else {

    ## -----------------------------------------------------
    ## 3. Run enrichment for each comparison
    ## -----------------------------------------------------
    kegg_results <- list()

    for (nm in names(limma_results)) {
      message("Running KEGG enrichment: ", nm)

      res <- as.data.frame(limma_results[[nm]], check.names = FALSE)

      res2 <- dplyr::left_join(
        res,
        annotation2[, c("metabolite_id", kegg_id_col), drop = FALSE],
        by = "metabolite_id"
      )

      universe_kegg <- unique(na.omit(annotation2[[kegg_id_col]]))

      sig_kegg <- res2 |>
        dplyr::filter(
          is.finite(P.Value),
          is.finite(logFC),
          P.Value < pvalue_cutoff,
          abs(logFC) >= logFC_cutoff
        ) |>
        dplyr::pull(!!kegg_id_col) |>
        unique() |>
        na.omit()

      if (length(sig_kegg) < 3) {
        warning("Too few significant KEGG compounds for ", nm, ". Skipped.")
        next
      }

      enr <- tryCatch(
        run_kegg_compound_enrichment(
          sig_kegg = sig_kegg,
          universe_kegg = universe_kegg
        ),
        error = function(e) {
          warning("KEGG enrichment failed for ", nm, ": ", e$message)
          NULL
        }
      )

      if (is.null(enr)) next

      enr_df <- as.data.frame(enr)
      if (nrow(enr_df) == 0) next

      enr_df$comparison <- nm
      kegg_results[[nm]] <- enr_df

      openxlsx::write.xlsx(
        enr_df,
        file = file.path(dir_kegg, paste0(safe_filename(nm), "_KEGG_enrichment.xlsx")),
        overwrite = TRUE
      )

      p <- tryCatch({
        clusterProfiler::dotplot(enr, showCategory = 20) +
          ggplot2::ggtitle("KEGG enrichment") +
          pub_theme(base_size = base_font_size, base_family = base_font_family)
      }, error = function(e) NULL)

      if (!is.null(p)) {
        tryCatch({
          save_ggplot_multi(
            p = p,
            filename_prefix = file.path(dir_kegg, paste0(safe_filename(nm), "_KEGG_dotplot")),
            width = 7,
            height = 6,
            dpi = plot_dpi
          )
        }, error = function(e) {
          warning("Saving KEGG dotplot failed for ", nm, ": ", e$message)
        })
      }
    }

    kegg_all <- dplyr::bind_rows(kegg_results)

    write_xlsx_list(
      list(
        merged_annotation = annotation2,
        KEGG_all = kegg_all
      ),
      file.path(dir_tables, "08_KEGG_output.xlsx")
    )

    saveRDS(kegg_results, file.path(dir_rds, "08_kegg_results.rds"))

    message("08_KEGG finished.")
  }
}
