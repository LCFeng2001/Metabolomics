## =========================================================
## 05_Limma.R
## v5: limma + P.Value preliminary screening + focus-class tables
## =========================================================

source("config/config.R")
source("functions/utils/common_utils.R")
source("functions/stats/run_limma.R")

data02 <- readRDS(file.path(dir_rds, "02_preprocessed.rds"))
comparisons <- readRDS(file.path(dir_rds, "04_comparisons.rds"))

limma_results <- run_limma_all(
  expr_log = data02$expr_log,
  sample_info = data02$sample_info,
  comparisons = comparisons,
  annotation = data02$annotation,
  adj_method = limma_adj_method,
  logFC_cutoff = logFC_cutoff,
  adjP_cutoff = adjP_cutoff,
  pvalue_cutoff = pvalue_cutoff,
  focus_classes = focus_classes,
  focus_bx_patterns = focus_bx_patterns
)

saveRDS(limma_results, file.path(dir_rds, "05_limma_results.rds"))

limma_all <- dplyr::bind_rows(limma_results)

## FDR-stringent significant metabolites for reference.
limma_sig_FDR <- limma_all[
  limma_all$adj.P.Val < adjP_cutoff & abs(limma_all$logFC) >= logFC_cutoff,
  ,
  drop = FALSE
]

## P.Value preliminary significant metabolites for volcano screening.
limma_sig_PValue <- limma_all[
  limma_all$P.Value < pvalue_cutoff & abs(limma_all$logFC) >= logFC_cutoff,
  ,
  drop = FALSE
]

## Focus classes requested by user.
limma_focus_all <- limma_all[
  !is.na(limma_all$is_focus_class) & limma_all$is_focus_class,
  ,
  drop = FALSE
]

limma_focus_PValue <- limma_sig_PValue[
  !is.na(limma_sig_PValue$is_focus_class) & limma_sig_PValue$is_focus_class,
  ,
  drop = FALSE
]

## Summary table: focus class x comparison.
if (nrow(limma_focus_all) > 0) {
  focus_summary <- limma_focus_all |>
    dplyr::mutate(
      prelim_significant = P.Value < pvalue_cutoff & abs(logFC) >= logFC_cutoff,
      stringent_FDR = adj.P.Val < adjP_cutoff & abs(logFC) >= logFC_cutoff
    ) |>
    dplyr::group_by(comparison, genotype, time, focus_class) |>
    dplyr::summarise(
      n_focus = dplyr::n(),
      n_PValue_sig = sum(prelim_significant, na.rm = TRUE),
      n_FDR_sig = sum(stringent_FDR, na.rm = TRUE),
      n_up_PValue = sum(PValue_regulation == "Up", na.rm = TRUE),
      n_down_PValue = sum(PValue_regulation == "Down", na.rm = TRUE),
      min_P.Value = min(P.Value, na.rm = TRUE),
      min_adj.P.Val = min(adj.P.Val, na.rm = TRUE),
      .groups = "drop"
    ) |>
    as.data.frame()
} else {
  focus_summary <- data.frame()
}

for (nm in names(limma_results)) {
  one <- limma_results[[nm]]
  one_focus <- one[!is.na(one$is_focus_class) & one$is_focus_class, , drop = FALSE]
  one_focus_P <- one_focus[
    one_focus$P.Value < pvalue_cutoff & abs(one_focus$logFC) >= logFC_cutoff,
    ,
    drop = FALSE
  ]

  write_xlsx_list(
    list(
      limma_all = one,
      PValue_preliminary = one[one$P.Value < pvalue_cutoff & abs(one$logFC) >= logFC_cutoff, , drop = FALSE],
      FDR_stringent = one[one$adj.P.Val < adjP_cutoff & abs(one$logFC) >= logFC_cutoff, , drop = FALSE],
      focus_all = one_focus,
      focus_PValue_preliminary = one_focus_P
    ),
    file.path(dir_limma, paste0(safe_filename(nm), "_limma.xlsx"))
  )
}

write_xlsx_list(
  list(
    comparisons = comparisons,
    limma_all = limma_all,
    PValue_preliminary_significant = limma_sig_PValue,
    FDR_stringent_significant = limma_sig_FDR,
    focus_all = limma_focus_all,
    focus_PValue_preliminary = limma_focus_PValue,
    focus_summary = focus_summary
  ),
  file.path(dir_tables, "05_Limma_output.xlsx")
)

message("05_Limma finished.")
message("Total comparisons: ", length(limma_results))
message("P.Value preliminary significant metabolites: ", nrow(limma_sig_PValue))
message("FDR stringent significant metabolites: ", nrow(limma_sig_FDR))
message("Focus-class metabolites: ", nrow(limma_focus_all))
message("Focus-class P.Value preliminary significant metabolites: ", nrow(limma_focus_PValue))
