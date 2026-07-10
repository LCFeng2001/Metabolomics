## =========================================================
## 02_Preprocess.R
## Blank filter, normalization, missing filtering, imputation,
## log2 transform and Pareto scaling
## v3: 0 is treated as missing value
## =========================================================

source("config/config.R")
source("functions/utils/common_utils.R")
source("functions/preprocess/preprocess_utils.R")

data01 <- readRDS(file.path(dir_rds, "01_clean_aggregated.rds"))

raw_missing_report <- make_missing_report(
  mat = data01$peak_area,
  sample_info = data01$sample_info,
  zero_as_missing = zero_as_missing
)

bf <- blank_filter(
  peak_area = data01$peak_area,
  blank_area = data01$blank_area,
  blank_fc = blank_fc,
  zero_as_missing = zero_as_missing
)

norm <- normalize_peak_area(
  mat = bf$peak_area,
  method = normalization_method,
  zero_as_missing = zero_as_missing
)

mf <- missing_filter(
  mat = norm$normalized,
  sample_info = data01$sample_info,
  min_valid_per_group = min_valid_per_group,
  zero_as_missing = zero_as_missing
)

post_filter_missing_report <- make_missing_report(
  mat = mf$filtered,
  sample_info = data01$sample_info,
  zero_as_missing = zero_as_missing
)

imp <- impute_missing(
  mat = mf$filtered,
  method = impute_method,
  zero_as_missing = zero_as_missing
)

expr_log <- log_transform(
  mat = imp$imputed,
  base = log_base,
  pseudo_count = pseudo_count
)

sc <- scale_metabolomics(
  mat = expr_log,
  method = scale_method
)

annotation2 <- data01$annotation[
  match(rownames(sc$scaled), data01$annotation$metabolite_id),
  ,
  drop = FALSE
]

data02 <- list(
  annotation = annotation2,
  sample_info = data01$sample_info,
  zero_as_missing = zero_as_missing,
  raw_missing_report = raw_missing_report,
  blank_filter_report = bf$report,
  normalization_factor = norm$factor,
  missing_filter_report = mf$report,
  post_filter_missing_report = post_filter_missing_report,
  filtered_area = mf$filtered,
  normalized_area = norm$normalized[rownames(mf$filtered), , drop = FALSE],
  imputed_area = imp$imputed,
  expr_log = expr_log,
  expr_scaled = sc$scaled,
  parameters = list(
    zero_as_missing = zero_as_missing,
    blank_fc = blank_fc,
    normalization_method = normalization_method,
    min_valid_per_group = min_valid_per_group,
    impute_method = impute_method,
    scale_method = scale_method
  )
)

saveRDS(data02, file.path(dir_rds, "02_preprocessed.rds"))

write_xlsx_list(
  list(
    annotation = annotation2,
    raw_missing_report_zero_as_NA = raw_missing_report,
    blank_filter_report = bf$report,
    normalization_factor = data.frame(
      sample = names(norm$factor),
      factor = as.numeric(norm$factor)
    ),
    missing_filter_report = mf$report,
    post_filter_missing_report = post_filter_missing_report,
    normalized_area = data.frame(
      metabolite_id = rownames(data02$normalized_area),
      data02$normalized_area,
      check.names = FALSE
    ),
    imputed_area = data.frame(
      metabolite_id = rownames(imp$imputed),
      imp$imputed,
      check.names = FALSE
    ),
    log2_expression = data.frame(
      metabolite_id = rownames(expr_log),
      expr_log,
      check.names = FALSE
    ),
    pareto_scaled = data.frame(
      metabolite_id = rownames(sc$scaled),
      sc$scaled,
      check.names = FALSE
    )
  ),
  file.path(dir_tables, "02_Preprocess_output.xlsx")
)

message("02_Preprocess finished.")
message("Zero as missing: ", zero_as_missing)
message("Metabolites after aggregation: ", nrow(data01$peak_area))
message("Metabolites after blank filter: ", nrow(bf$peak_area))
message("Metabolites after missing filter: ", nrow(mf$filtered))
message("Raw zero count: ", sum(raw_missing_report$zero_count, na.rm = TRUE))
