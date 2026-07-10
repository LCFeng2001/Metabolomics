## =========================================================
## read_clean_aggregate.R
## Read Skyline raw Peak Area table, clean columns,
## and aggregate duplicate metabolites by Compound name.
## =========================================================

read_raw_table <- function(file, sheet = 1) {
  suppressPackageStartupMessages({
    library(readxl)
    library(readr)
    library(dplyr)
    library(stringr)
  })

  if (!file.exists(file)) {
    stop("Input file not found: ", file)
  }

  ext <- tolower(tools::file_ext(file))

  if (ext %in% c("xlsx", "xls")) {
    df <- readxl::read_excel(
      file,
      sheet = sheet,
      .name_repair = "unique"
    )
  } else if (ext %in% c("txt", "tsv")) {
    df <- readr::read_tsv(
      file,
      show_col_types = FALSE,
      name_repair = "unique"
    )
  } else if (ext %in% c("csv")) {
    df <- readr::read_csv(
      file,
      show_col_types = FALSE,
      name_repair = "unique"
    )
  } else {
    stop("Unsupported file extension: ", ext)
  }

  df <- as.data.frame(df, check.names = FALSE)
  names(df) <- stringr::str_squish(names(df))

  df
}

parse_sample_info <- function(sample_cols, sample_regex) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(tibble)
  })

  m <- stringr::str_match(sample_cols, sample_regex)

  if (any(is.na(m[, 1]))) {
    bad <- sample_cols[is.na(m[, 1])]
    stop("Unparsed sample columns:\n", paste(bad, collapse = "\n"))
  }

  tibble(
    sample = sample_cols,
    genotype = m[, 2],
    treatment = m[, 3],
    time = m[, 4],
    rep = as.integer(m[, 5])
  ) %>%
    mutate(
      time_num = as.numeric(stringr::str_remove(time, "h")),
      group = paste(genotype, time, treatment, sep = "_")
    ) %>%
    arrange(genotype, time_num, treatment, rep)
}

clean_and_aggregate_metabolites <- function(
    df,
    compound_col = "Compound name",
    feature_col = "分子",
    class_col = "Class",
    sample_regex,
    blank_regex,
    drop_exact_cols = c("序号", "母离子质荷比", "子离子质荷比", "blank均值"),
    aggregate_by = "Compound name"
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(tibble)
  })

  original_colnames <- names(df)
  base_colnames <- sub("\\.\\.\\.[0-9]+$", "", original_colnames)

  ## Identify columns
  sample_cols <- original_colnames[stringr::str_detect(original_colnames, sample_regex)]
  blank_cols  <- original_colnames[stringr::str_detect(base_colnames, blank_regex)]

  mean_cols <- original_colnames[base_colnames == "均值" | base_colnames == "blank均值"]
  useless_cols <- original_colnames[base_colnames %in% drop_exact_cols]
  useless_cols <- unique(c(useless_cols, mean_cols))

  keep_meta <- intersect(c(feature_col, compound_col, class_col), original_colnames)
  keep_cols <- unique(c(keep_meta, blank_cols, sample_cols))

  if (!compound_col %in% original_colnames) {
    stop("Cannot find compound column: ", compound_col)
  }
  if (!feature_col %in% original_colnames) {
    warning("Cannot find feature column: ", feature_col, ". It will be generated.")
    df[[feature_col]] <- paste0("feature_", seq_len(nrow(df)))
    keep_meta <- unique(c(feature_col, keep_meta))
    keep_cols <- unique(c(feature_col, keep_cols))
  }
  if (!class_col %in% original_colnames) {
    df[[class_col]] <- NA_character_
    keep_meta <- unique(c(keep_meta, class_col))
    keep_cols <- unique(c(keep_meta, blank_cols, sample_cols))
  }

  if (length(sample_cols) == 0) {
    stop("No biological sample columns detected. Please check sample_regex.")
  }

  clean_df <- df[, keep_cols, drop = FALSE]

  numeric_cols <- unique(c(blank_cols, sample_cols))

  for (cc in numeric_cols) {
    clean_df[[cc]] <- suppressWarnings(as.numeric(as.character(clean_df[[cc]])))
  }

  ## Remove rows without compound names
  clean_df[[compound_col]] <- as.character(clean_df[[compound_col]])
  clean_df[[feature_col]]  <- as.character(clean_df[[feature_col]])
  clean_df[[class_col]]    <- as.character(clean_df[[class_col]])

  clean_df[[compound_col]][is.na(clean_df[[compound_col]]) | clean_df[[compound_col]] == ""] <-
    clean_df[[feature_col]][is.na(clean_df[[compound_col]]) | clean_df[[compound_col]] == ""]

  clean_df <- clean_df[!is.na(clean_df[[compound_col]]) & clean_df[[compound_col]] != "", , drop = FALSE]

  ## Aggregate by Compound name
  aggregate_col <- aggregate_by
  if (!aggregate_col %in% colnames(clean_df)) {
    stop("aggregate_by column not found: ", aggregate_col)
  }

  message("Detected sample columns: ", length(sample_cols))
  message("Detected blank columns: ", length(blank_cols))
  message("Dropped useless / mean columns: ", length(useless_cols))
  message("Rows before duplicate aggregation: ", nrow(clean_df))

  agg_df <- clean_df %>%
    group_by(.data[[aggregate_col]]) %>%
    summarise(
      source_features = paste_unique_non_empty(.data[[feature_col]], sep = ";"),
      `Compound name` = first_non_empty(.data[[compound_col]]),
      Class = paste_unique_non_empty(.data[[class_col]], sep = ";"),
      across(all_of(numeric_cols), sum_na),
      .groups = "drop"
    )

  message("Rows after duplicate aggregation: ", nrow(agg_df))
  message("Merged duplicate rows: ", nrow(clean_df) - nrow(agg_df))

  ## Create metabolite_id
  agg_df$metabolite_id <- make.unique(as.character(agg_df[["Compound name"]]))
  agg_df <- agg_df[, c("metabolite_id", "source_features", "Compound name", "Class", numeric_cols), drop = FALSE]

  ## Important: agg_df is a tibble after dplyr::summarise().
  ## Do NOT set rownames directly on a tibble; convert to plain data.frame first.
  annotation <- as.data.frame(
    agg_df[, c("metabolite_id", "source_features", "Compound name", "Class"), drop = FALSE],
    check.names = FALSE
  )

  peak_area_df <- as.data.frame(agg_df[, sample_cols, drop = FALSE], check.names = FALSE)
  rownames(peak_area_df) <- agg_df$metabolite_id
  peak_area <- as_numeric_matrix(peak_area_df)
  rownames(peak_area) <- agg_df$metabolite_id

  if (length(blank_cols) > 0) {
    blank_area_df <- as.data.frame(agg_df[, blank_cols, drop = FALSE], check.names = FALSE)
    rownames(blank_area_df) <- agg_df$metabolite_id
    blank_area <- as_numeric_matrix(blank_area_df)
    rownames(blank_area) <- agg_df$metabolite_id
  } else {
    blank_area <- NULL
  }

  sample_info <- parse_sample_info(sample_cols, sample_regex)
  peak_area <- peak_area[, sample_info$sample, drop = FALSE]

  list(
    clean_df = clean_df,
    aggregated_df = agg_df,
    annotation = annotation,
    peak_area = peak_area,
    blank_area = blank_area,
    sample_info = sample_info,
    sample_cols = sample_cols,
    blank_cols = blank_cols,
    removed_cols = useless_cols
  )
}
