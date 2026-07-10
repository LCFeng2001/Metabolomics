## =========================================================
## 04_Comparisons.R
## =========================================================

source("config/config.R")
source("functions/utils/common_utils.R")
source("functions/utils/comparison_generator.R")

data02 <- readRDS(file.path(dir_rds, "02_preprocessed.rds"))

comparisons <- generate_fg_ck_comparisons(
  sample_info = data02$sample_info,
  times = comparison_times,
  case_treatment = case_treatment,
  control_treatment = control_treatment,
  min_n = min_valid_per_group
)

saveRDS(comparisons, file.path(dir_rds, "04_comparisons.rds"))

write_xlsx_list(
  list(comparisons = comparisons),
  file.path(dir_tables, "04_Comparisons_output.xlsx")
)

message("04_Comparisons finished.")
message("Total comparisons: ", nrow(comparisons))
