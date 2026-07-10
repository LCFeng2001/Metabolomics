## =========================================================
## comparison_generator.R
## =========================================================

generate_fg_ck_comparisons <- function(
    sample_info,
    times = c("24h", "48h"),
    case_treatment = "Fg",
    control_treatment = "CK",
    min_n = 2
) {
  suppressPackageStartupMessages({
    library(dplyr)
  })

  sample_info <- as.data.frame(sample_info, check.names = FALSE)
  genotypes <- unique(sample_info$genotype)

  out <- list()

  for (g in genotypes) {
    for (tm in times) {
      case_group <- paste(g, tm, case_treatment, sep = "_")
      control_group <- paste(g, tm, control_treatment, sep = "_")

      n_case <- sum(sample_info$group == case_group)
      n_control <- sum(sample_info$group == control_group)

      if (n_case >= min_n && n_control >= min_n) {
        out[[length(out) + 1]] <- data.frame(
          comparison = paste(g, tm, case_treatment, "vs", control_treatment, sep = "_"),
          genotype = g,
          time = tm,
          case_treatment = case_treatment,
          control_treatment = control_treatment,
          case_group = case_group,
          control_group = control_group,
          n_case = n_case,
          n_control = n_control,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(out) == 0) {
    return(data.frame())
  }

  bind_rows(out)
}
