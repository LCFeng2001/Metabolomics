## =========================================================
## run_oplsda.R
## Robust OPLS-DA + VIP calculation
## v14:
##   - try ropls::opls first
##   - if ropls fails, generate a robust fallback VIP-like ranking
##   - never returns empty when enough samples and variables exist
## =========================================================

prepare_oplsda_input <- function(
    expr_scaled,
    sample_info,
    case_group,
    control_group
) {
  expr_scaled <- as_numeric_matrix(expr_scaled)
  sample_info <- as.data.frame(sample_info, check.names = FALSE)

  selected_samples <- sample_info$sample[
    sample_info$group %in% c(control_group, case_group)
  ]
  selected_samples <- intersect(selected_samples, colnames(expr_scaled))

  if (length(selected_samples) < 4) {
    stop("Too few selected samples: ", length(selected_samples))
  }

  sub_info <- sample_info[match(selected_samples, sample_info$sample), , drop = FALSE]
  sub_info <- sub_info[!is.na(sub_info$sample), , drop = FALSE]

  y <- factor(
    ifelse(sub_info$group == case_group, "Case", "Control"),
    levels = c("Control", "Case")
  )

  if (length(unique(y)) < 2) {
    stop("Only one group found after sample selection.")
  }

  n_case <- sum(y == "Case", na.rm = TRUE)
  n_control <- sum(y == "Control", na.rm = TRUE)

  if (n_case < 2 || n_control < 2) {
    stop("Each group needs at least 2 samples for OPLS-DA/VIP. n_case=", n_case, ", n_control=", n_control)
  }

  x <- t(expr_scaled[, sub_info$sample, drop = FALSE])
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"

  ## Replace any remaining non-finite values by feature median, then 0.
  for (j in seq_len(ncol(x))) {
    v <- x[, j]
    if (any(!is.finite(v))) {
      med <- median(v[is.finite(v)], na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      v[!is.finite(v)] <- med
      x[, j] <- v
    }
  }

  keep <- apply(x, 2, function(v) {
    s <- sd(v, na.rm = TRUE)
    is.finite(s) && s > 0
  })

  x <- x[, keep, drop = FALSE]

  if (ncol(x) < 3) {
    stop("Too few non-constant variables after filtering: ", ncol(x))
  }

  list(
    x = x,
    y = y,
    sample_info = sub_info,
    n_case = n_case,
    n_control = n_control
  )
}

calculate_fallback_vip <- function(x, y) {
  ## A robust VIP-like ranking based on standardized case-control separation.
  ## This is not a formal OPLS-DA VIP, but keeps downstream prioritization usable
  ## when ropls fails because of small sample size or unstable model fitting.
  y <- factor(y, levels = c("Control", "Case"))

  case_idx <- y == "Case"
  ctrl_idx <- y == "Control"

  mean_case <- colMeans(x[case_idx, , drop = FALSE], na.rm = TRUE)
  mean_ctrl <- colMeans(x[ctrl_idx, , drop = FALSE], na.rm = TRUE)

  sd_case <- apply(x[case_idx, , drop = FALSE], 2, sd, na.rm = TRUE)
  sd_ctrl <- apply(x[ctrl_idx, , drop = FALSE], 2, sd, na.rm = TRUE)

  n_case <- sum(case_idx)
  n_ctrl <- sum(ctrl_idx)

  pooled_sd <- sqrt(((n_case - 1) * sd_case^2 + (n_ctrl - 1) * sd_ctrl^2) / max(1, n_case + n_ctrl - 2))
  pooled_sd[!is.finite(pooled_sd) | pooled_sd == 0] <- NA_real_

  score <- abs(mean_case - mean_ctrl) / pooled_sd
  score[!is.finite(score)] <- abs(mean_case[!is.finite(score)] - mean_ctrl[!is.finite(score)])

  score[!is.finite(score)] <- 0

  ## Scale to a VIP-like range where average squared score is approximately 1.
  if (sum(score^2, na.rm = TRUE) > 0) {
    vip_like <- score / sqrt(mean(score^2, na.rm = TRUE))
  } else {
    vip_like <- rep(0, length(score))
  }

  data.frame(
    metabolite_id = names(vip_like),
    VIP = as.numeric(vip_like),
    VIP_method = "fallback_standardized_difference",
    stringsAsFactors = FALSE
  )[order(vip_like, decreasing = TRUE), , drop = FALSE]
}

run_oplsda_single <- function(
    expr_scaled,
    sample_info,
    case_group,
    control_group,
    permI = 200,
    fallback = TRUE
) {
  prep <- prepare_oplsda_input(
    expr_scaled = expr_scaled,
    sample_info = sample_info,
    case_group = case_group,
    control_group = control_group
  )

  x <- prep$x
  y <- prep$y

  model <- NULL
  vip_df <- NULL
  summary_df <- NULL
  status <- "not_run"
  message_txt <- NA_character_

  if (requireNamespace("ropls", quietly = TRUE)) {
    ## Large permI can fail or be slow for small n. Keep it sane.
    permI_use <- permI
    if (!is.finite(permI_use) || is.na(permI_use)) permI_use <- 0
    permI_use <- as.integer(max(0, permI_use))

    model <- tryCatch(
      ropls::opls(
        x,
        y,
        predI = 1,
        orthoI = NA,
        permI = permI_use,
        fig.pdfC = "none",
        plotL = FALSE,
        printL = FALSE
      ),
      error = function(e) {
        message_txt <<- e$message
        NULL
      },
      warning = function(w) {
        ## keep going on warnings
        invokeRestart("muffleWarning")
      }
    )

    if (!is.null(model)) {
      vip <- tryCatch(
        ropls::getVipVn(model),
        error = function(e) NULL
      )

      if (is.null(vip)) {
        vip <- tryCatch(model@vipVn, error = function(e) NULL)
      }

      if (!is.null(vip) && length(vip) > 0) {
        vip_df <- data.frame(
          metabolite_id = names(vip),
          VIP = as.numeric(vip),
          VIP_method = "ropls",
          stringsAsFactors = FALSE
        )
        vip_df <- vip_df[is.finite(vip_df$VIP), , drop = FALSE]
        vip_df <- vip_df[order(vip_df$VIP, decreasing = TRUE), , drop = FALSE]
      }

      summary_df <- tryCatch(
        as.data.frame(ropls::getSummaryDF(model)),
        error = function(e) NULL
      )

      status <- if (!is.null(vip_df) && nrow(vip_df) > 0) "ropls_success" else "ropls_no_vip"
      message_txt <- ifelse(is.na(message_txt), "ropls model fitted", message_txt)
    } else {
      status <- "ropls_failed"
      if (is.na(message_txt)) message_txt <- "ropls model failed"
    }
  } else {
    status <- "ropls_not_installed"
    message_txt <- "Package ropls is not installed"
  }

  if ((is.null(vip_df) || nrow(vip_df) == 0) && isTRUE(fallback)) {
    vip_df <- calculate_fallback_vip(x, y)
    status <- paste0(status, "_fallback_vip")
    if (is.na(message_txt)) message_txt <- "fallback VIP-like score generated"
  }

  if (is.null(summary_df) || nrow(summary_df) == 0) {
    summary_df <- data.frame(
      status = status,
      message = message_txt,
      n_samples = nrow(x),
      n_variables = ncol(x),
      n_case = prep$n_case,
      n_control = prep$n_control,
      stringsAsFactors = FALSE
    )
  } else {
    summary_df$status <- status
    summary_df$message <- message_txt
    summary_df$n_samples <- nrow(x)
    summary_df$n_variables <- ncol(x)
    summary_df$n_case <- prep$n_case
    summary_df$n_control <- prep$n_control
  }

  list(
    model = model,
    vip = vip_df,
    summary = summary_df,
    status = status,
    message = message_txt
  )
}
