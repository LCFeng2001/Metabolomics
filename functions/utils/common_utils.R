## =========================================================
## common_utils.R
## v5 utilities: robust numeric conversion, zero-as-missing,
## focus-class detection, publication theme
## =========================================================

safe_filename <- function(x) {
  x <- as.character(x)
  x <- gsub("[/\\\\:*?\"<>| ]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

sum_na <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  if (all(is.na(x))) return(NA_real_)
  sum(x, na.rm = TRUE)
}

first_non_empty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

paste_unique_non_empty <- function(x, sep = ";") {
  x <- as.character(x)
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = sep)
}

as_numeric_matrix <- function(x) {
  if (is.data.frame(x) || inherits(x, "tbl_df")) {
    rn <- rownames(x)
    x <- as.data.frame(x, check.names = FALSE)
    rownames(x) <- rn
  }
  rn <- rownames(x)
  cn <- colnames(x)
  x <- as.matrix(x)
  suppressWarnings(storage.mode(x) <- "numeric")
  rownames(x) <- rn
  colnames(x) <- cn
  x
}

apply_zero_as_missing <- function(mat, zero_as_missing = TRUE) {
  mat <- as_numeric_matrix(mat)
  if (isTRUE(zero_as_missing)) {
    mat[mat <= 0] <- NA_real_
  }
  mat
}

make_missing_report <- function(mat, sample_info = NULL, zero_as_missing = TRUE) {
  raw <- as_numeric_matrix(mat)
  zero_count <- rowSums(raw == 0, na.rm = TRUE)
  negative_count <- rowSums(raw < 0, na.rm = TRUE)
  na_count <- rowSums(is.na(raw))

  mat2 <- apply_zero_as_missing(raw, zero_as_missing = zero_as_missing)
  missing_count <- rowSums(is.na(mat2))
  valid_count <- rowSums(!is.na(mat2))

  out <- data.frame(
    metabolite_id = rownames(raw),
    zero_count = zero_count,
    negative_count = negative_count,
    original_na_count = na_count,
    missing_count_after_zero_rule = missing_count,
    valid_count_after_zero_rule = valid_count,
    missing_rate_after_zero_rule = missing_count / ncol(raw),
    stringsAsFactors = FALSE
  )

  if (!is.null(sample_info)) {
    sample_info <- as.data.frame(sample_info, check.names = FALSE)
    sample_info <- sample_info[match(colnames(raw), sample_info$sample), , drop = FALSE]
    groups <- unique(sample_info$group)
    for (g in groups) {
      cols <- sample_info$sample[sample_info$group == g]
      out[[paste0("valid_", g)]] <- rowSums(!is.na(mat2[, cols, drop = FALSE]))
      out[[paste0("missing_", g)]] <- rowSums(is.na(mat2[, cols, drop = FALSE]))
    }
  }

  out
}

write_xlsx_list <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  openxlsx::write.xlsx(x, file = file, overwrite = TRUE)
}

## ---------------------------------------------------------
## Focus-class utilities
## ---------------------------------------------------------
normalize_class_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  ## Replace non-breaking spaces and repeated whitespace.
  x <- gsub("\\u00A0", " ", x, fixed = FALSE)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

split_class_terms <- function(x) {
  x <- normalize_class_text(x)
  out <- strsplit(x, ";|,|/", perl = TRUE)
  lapply(out, function(v) trimws(v[v != ""]))
}

standardize_focus_class <- function(class_vec, focus_classes, bx_patterns = c("^BX$", "Benzoxazinoid", "Benzoxazinoids")) {
  class_vec <- normalize_class_text(class_vec)
  focus_classes <- normalize_class_text(focus_classes)

  out <- rep(NA_character_, length(class_vec))

  for (i in seq_along(class_vec)) {
    terms <- split_class_terms(class_vec[i])[[1]]
    if (length(terms) == 0) next

    ## Exact class matching first.
    hit <- focus_classes[match(tolower(focus_classes), tolower(terms), nomatch = 0) > 0]

    ## Partial robust matching for terms merged with semicolon or naming variants.
    if (length(hit) == 0) {
      for (fc in focus_classes) {
        if (any(tolower(terms) == tolower(fc))) {
          hit <- fc
          break
        }
      }
    }

    ## BX may appear as BX, Benzoxazinoid, Benzoxazinoids, etc.
    bx_hit <- any(vapply(bx_patterns, function(p) any(grepl(p, terms, ignore.case = TRUE)), logical(1)))
    if (bx_hit) {
      hit <- c(hit, "BX")
    }

    hit <- unique(hit)
    if (length(hit) > 0) out[i] <- paste(hit, collapse = ";")
  }

  out
}

add_focus_class_columns <- function(df, focus_classes, class_col = "Class", bx_patterns = c("^BX$", "Benzoxazinoid", "Benzoxazinoids")) {
  df <- as.data.frame(df, check.names = FALSE)
  if (!class_col %in% colnames(df)) {
    df[[class_col]] <- NA_character_
  }
  df$Class_clean <- normalize_class_text(df[[class_col]])
  df$focus_class <- standardize_focus_class(
    df[[class_col]],
    focus_classes = focus_classes,
    bx_patterns = bx_patterns
  )
  df$is_focus_class <- !is.na(df$focus_class) & df$focus_class != ""
  df
}

## Publication-level ggplot theme
pub_theme <- function(base_size = 12, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = base_size, color = "grey30"),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(color = "black", size = base_size - 1),
      axis.line = ggplot2::element_line(linewidth = 0.45, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.4, color = "black"),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key.size = grid::unit(0.55, "cm"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "black", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold", color = "black"),
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.45)
    )
}

## Consistent palettes. These are color-blind friendly Okabe-Ito inspired colors.
metabo_palette <- function(type = c("treatment", "significance", "time", "focus_class")) {
  type <- match.arg(type)
  if (type == "treatment") {
    return(c(CK = "#0072B2", Fg = "#D55E00"))
  }
  if (type == "significance") {
    return(c(Down = "#0072B2", NS = "grey78", Up = "#D55E00"))
  }
  if (type == "time") {
    return(c(`0h` = "#999999", `24h` = "#009E73", `48h` = "#CC79A7"))
  }
  if (type == "focus_class") {
    return(c(
      Phenolamides = "#D55E00",
      `Benzoic acid derivatives` = "#0072B2",
      `Hydroxycinnamoyl derivatives` = "#009E73",
      Flavone = "#CC79A7",
      Flavonol = "#E69F00",
      Flavanone = "#56B4E9",
      `Flavone C-glycosides` = "#F0E442",
      BX = "#000000",
      Other_focus = "#7F7F7F"
    ))
  }
}
