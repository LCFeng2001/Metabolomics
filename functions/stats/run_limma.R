## =========================================================
## run_limma.R
## =========================================================

run_limma_all <- function(
    expr_log,
    sample_info,
    comparisons,
    annotation = NULL,
    adj_method = "BH",
    logFC_cutoff = log2(1.5),
    adjP_cutoff = 0.05,
    pvalue_cutoff = 0.05,
    focus_classes = NULL,
    focus_bx_patterns = c("^BX$", "Benzoxazinoid", "Benzoxazinoids")
) {
  suppressPackageStartupMessages({
    library(limma)
    library(dplyr)
  })

  expr_log <- as_numeric_matrix(expr_log)
  sample_info <- as.data.frame(sample_info, check.names = FALSE)
  comparisons <- as.data.frame(comparisons, check.names = FALSE)

  sample_info <- sample_info[match(colnames(expr_log), sample_info$sample), , drop = FALSE]

  group <- factor(sample_info$group)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- make.names(levels(group))

  fit <- limma::lmFit(expr_log, design)

  out <- list()

  for (i in seq_len(nrow(comparisons))) {
    comp <- comparisons[i, ]

    case_design <- make.names(comp$case_group)
    ctrl_design <- make.names(comp$control_group)

    contrast <- paste0(case_design, " - ", ctrl_design)

    cont <- limma::makeContrasts(
      contrasts = contrast,
      levels = design
    )

    fit2 <- limma::contrasts.fit(fit, cont)
    fit2 <- limma::eBayes(fit2)

    res <- limma::topTable(
      fit2,
      number = Inf,
      adjust.method = adj_method,
      sort.by = "P"
    )

    res$metabolite_id <- rownames(res)
    res$comparison <- comp$comparison
    res$genotype <- comp$genotype
    res$time <- comp$time
    res$case_group <- comp$case_group
    res$control_group <- comp$control_group

    case_samples <- sample_info$sample[sample_info$group == comp$case_group]
    ctrl_samples <- sample_info$sample[sample_info$group == comp$control_group]

    res$mean_case_log2 <- rowMeans(expr_log[res$metabolite_id, case_samples, drop = FALSE], na.rm = TRUE)
    res$mean_control_log2 <- rowMeans(expr_log[res$metabolite_id, ctrl_samples, drop = FALSE], na.rm = TRUE)
    res$FC <- 2^res$logFC

    ## FDR-based regulation is retained for stringent reference.
    res$FDR_regulation <- ifelse(
      res$adj.P.Val < adjP_cutoff & res$logFC >= logFC_cutoff,
      "Up",
      ifelse(
        res$adj.P.Val < adjP_cutoff & res$logFC <= -logFC_cutoff,
        "Down",
        "NS"
      )
    )

    ## P.Value-based preliminary regulation is used for volcano screening.
    res$PValue_regulation <- ifelse(
      res$P.Value < pvalue_cutoff & res$logFC >= logFC_cutoff,
      "Up",
      ifelse(
        res$P.Value < pvalue_cutoff & res$logFC <= -logFC_cutoff,
        "Down",
        "NS"
      )
    )

    ## Backward-compatible column; now uses preliminary P.Value regulation by default.
    res$regulation <- res$PValue_regulation

    if (!is.null(annotation)) {
      annotation <- as.data.frame(annotation, check.names = FALSE)
      res <- dplyr::left_join(res, annotation, by = "metabolite_id")
    }

    if (!is.null(focus_classes)) {
      res <- add_focus_class_columns(
        res,
        focus_classes = focus_classes,
        class_col = "Class",
        bx_patterns = focus_bx_patterns
      )
    }

    out[[comp$comparison]] <- res
  }

  out
}
