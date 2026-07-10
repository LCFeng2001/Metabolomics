## =========================================================
## 01_ReadCleanAggregate.R
## Read raw file, drop useless columns, aggregate duplicates
## =========================================================

source("config/config.R")
source("functions/utils/common_utils.R")
source("functions/io/read_clean_aggregate.R")

raw_df <- read_raw_table(
  file = raw_file,
  sheet = raw_sheet
)

data01 <- clean_and_aggregate_metabolites(
  df = raw_df,
  compound_col = compound_col,
  feature_col = feature_col,
  class_col = class_col,
  sample_regex = sample_regex,
  blank_regex = blank_regex,
  drop_exact_cols = drop_exact_cols,
  aggregate_by = aggregate_by
)

saveRDS(data01, file.path(dir_rds, "01_clean_aggregated.rds"))

write_xlsx_list(
  list(
    removed_columns = data.frame(removed_cols = data01$removed_cols),
    sample_info = data01$sample_info,
    annotation_aggregated = data01$annotation,
    peak_area_aggregated = data.frame(
      metabolite_id = data01$annotation$metabolite_id,
      data01$peak_area,
      check.names = FALSE
    ),
    blank_area_aggregated = if (!is.null(data01$blank_area)) {
      data.frame(
        metabolite_id = data01$annotation$metabolite_id,
        data01$blank_area,
        check.names = FALSE
      )
    } else {
      data.frame()
    }
  ),
  file.path(dir_tables, "01_ReadCleanAggregate_output.xlsx")
)

message("01_ReadCleanAggregate finished.")
