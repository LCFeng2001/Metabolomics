## =========================================================
## preprocess_utils.R
## v3: 0 means missing; blank filter, normalization, missing filter,
## imputation, log transformation and scaling
## =========================================================

blank_filter <- function(
    peak_area,
    blank_area = NULL,
    blank_fc = 3,
    zero_as_missing = TRUE
) {
  mat <- apply_zero_as_missing(peak_area, zero_as_missing = zero_as_missing)

  sample_mean <- rowMeans(mat, na.rm = TRUE)
  sample_mean[is.nan(sample_mean)] <- NA_real_

  sample_max <- apply(mat, 1, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    max(x, na.rm = TRUE)
  })

  if (is.null(blank_area)) {
    blank_mean <- rep(NA_real_, nrow(mat))
    names(blank_mean) <- rownames(mat)
    keep_blank <- rep(TRUE, nrow(mat))
  } else {
    bmat <- apply_zero_as_missing(blank_area, zero_as_missing = zero_as_missing)
    blank_mean <- rowMeans(bmat, na.rm = TRUE)
    blank_mean[is.nan(blank_mean)] <- NA_real_

    keep_blank <- is.na(blank_mean) |
      blank_mean == 0 |
      sample_mean >= blank_fc * blank_mean
  }

  keep_signal <- !is.na(sample_mean) & sample_mean > 0
  keep <- keep_blank & keep_signal

  report <- data.frame(
    metabolite_id = rownames(mat),
    sample_mean = sample_mean,
    sample_max = sample_max,
    blank_mean = blank_mean,
    sample_blank_ratio = sample_mean / blank_mean,
    keep_blank = keep_blank,
    keep_signal = keep_signal,
    keep = keep,
    stringsAsFactors = FALSE
  )

  list(
    peak_area = mat[keep, , drop = FALSE],
    report = report
  )
}

normalize_peak_area <- function(
    mat,
    method = c("none", "TIC", "median"),
    zero_as_missing = TRUE
) {
  method <- match.arg(method)
  mat <- apply_zero_as_missing(mat, zero_as_missing = zero_as_missing)

  if (method == "none") {
    factor <- rep(1, ncol(mat))
    names(factor) <- colnames(mat)
    norm <- mat
  }

  if (method == "TIC") {
    sample_sum <- colSums(mat, na.rm = TRUE)
    target <- median(sample_sum[sample_sum > 0], na.rm = TRUE)
    factor <- sample_sum / target
    factor[!is.finite(factor) | factor == 0] <- 1
    norm <- sweep(mat, 2, factor, "/")
  }

  if (method == "median") {
    sample_median <- apply(mat, 2, function(x) {
      x <- x[is.finite(x) & x > 0]
      if (length(x) == 0) return(NA_real_)
      median(x, na.rm = TRUE)
    })
    target <- median(sample_median, na.rm = TRUE)
    factor <- sample_median / target
    factor[!is.finite(factor) | factor == 0] <- 1
    norm <- sweep(mat, 2, factor, "/")
  }

  norm <- apply_zero_as_missing(norm, zero_as_missing = zero_as_missing)

  list(normalized = norm, factor = factor, method = method)
}

missing_filter <- function(
    mat,
    sample_info,
    min_valid_per_group = 2,
    zero_as_missing = TRUE
) {
  mat <- apply_zero_as_missing(mat, zero_as_missing = zero_as_missing)

  sample_info <- as.data.frame(sample_info, check.names = FALSE)
  sample_info <- sample_info[match(colnames(mat), sample_info$sample), , drop = FALSE]

  det <- !is.na(mat) & mat > 0
  groups <- unique(sample_info$group)

  valid_count <- sapply(groups, function(g) {
    cols <- sample_info$sample[sample_info$group == g]
    rowSums(det[, cols, drop = FALSE])
  })

  valid_count <- as.data.frame(valid_count, check.names = FALSE)
  valid_count$metabolite_id <- rownames(mat)

  keep <- apply(valid_count[, groups, drop = FALSE] >= min_valid_per_group, 1, any)

  report <- valid_count
  report$total_valid <- rowSums(det)
  report$total_missing <- rowSums(!det)
  report$missing_rate <- report$total_missing / ncol(mat)
  report$keep <- keep

  list(
    filtered = mat[keep, , drop = FALSE],
    report = report
  )
}

impute_missing <- function(
    mat,
    method = c("half_min", "none"),
    zero_as_missing = TRUE
) {
  method <- match.arg(method)

  mat <- apply_zero_as_missing(mat, zero_as_missing = zero_as_missing)

  if (method == "none") {
    return(list(imputed = mat, method = method))
  }

  global_min <- min(mat[is.finite(mat) & mat > 0], na.rm = TRUE)

  if (!is.finite(global_min)) {
    stop("No positive values available for imputation.")
  }

  imp <- t(apply(mat, 1, function(x) {
    min_pos <- min(x[is.finite(x) & x > 0], na.rm = TRUE)
    if (!is.finite(min_pos)) min_pos <- global_min
    x[is.na(x)] <- min_pos / 2
    x
  }))

  rownames(imp) <- rownames(mat)
  colnames(imp) <- colnames(mat)

  list(imputed = imp, method = method)
}

log_transform <- function(mat, base = 2, pseudo_count = 1) {
  mat <- as_numeric_matrix(mat)
  mat[!is.finite(mat)] <- NA_real_
  log(mat + pseudo_count, base = base)
}

scale_metabolomics <- function(mat, method = c("pareto", "auto", "none")) {
  method <- match.arg(method)
  mat <- as_numeric_matrix(mat)

  scaled <- t(apply(mat, 1, function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)

    if (!is.finite(s) || s == 0) {
      return(rep(0, length(x)))
    }

    if (method == "pareto") {
      return((x - m) / sqrt(s))
    }

    if (method == "auto") {
      return((x - m) / s)
    }

    x
  }))

  rownames(scaled) <- rownames(mat)
  colnames(scaled) <- colnames(mat)
  scaled[!is.finite(scaled)] <- 0

  list(scaled = scaled, method = method)
}
