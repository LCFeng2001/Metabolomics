## =========================================================
## config/config.R
## Rebuilt configuration for Skyline raw Peak Area metabolomics
## v6: P.Value volcano + focus-class boxplot + VIP panel combined with volcano
## =========================================================

project_dir <- getwd()

## Input file
## For real data, change this to:
raw_file <- file.path(project_dir, "data/raw/skyline_peak_area.xlsx")
## raw_file <- file.path(project_dir, "data/raw/example_peak_area.tsv")
raw_sheet <- 1

## Column rules
compound_col <- "Compound name"
feature_col  <- "分子"
class_col    <- "Class"

## These useless columns will be removed at the first cleaning step
drop_exact_cols <- c("序号", "母离子质荷比", "子离子质荷比", "blank均值")

## Biological sample format:
## D3140-CK 24h-1
## M03-Fg 48h-3
sample_regex <- "^(.+)-(CK|Fg)\\s*(0h|24h|48h)-(\\d+)$"

## Blank replicate format:
## blank004, blank005
blank_regex <- "^blank\\d+$"

## Duplicate aggregation method
## "Compound name" is recommended because the same metabolite may have multiple transitions / IDs.
aggregate_by <- "Compound name"

## =========================================================
## Missing-value logic
## =========================================================
## Your Skyline empty values are encoded as 0.
## TRUE means 0 or negative Peak Area will be treated as missing NA from preprocessing onward.
zero_as_missing <- TRUE

## After normalization and filtering, missing values are imputed by half of the row-wise minimum positive value.
impute_method <- "half_min"  ## "half_min" or "none"

## Blank filtering
blank_fc <- 3

## Missing filtering
## Retain a metabolite if any biological group has at least this many non-missing values.
min_valid_per_group <- 2

## Normalization method:
## "none", "TIC", or "median"
## For raw Skyline Peak Area without internal standards, median normalization is usually a safe default.
normalization_method <- "median"

## Transformation and scaling
log_base <- 2
pseudo_count <- 1
scale_method <- "pareto"  ## "pareto", "auto", "none"

## Differential analysis thresholds
logFC_cutoff <- 1
adjP_cutoff <- 0.05
limma_adj_method <- "BH"
## Volcano preliminary screening threshold
## For broad targeted metabolomics screening, volcano plots use nominal P.Value by default.
## FDR remains in limma result tables for reference, but is not used to color/label volcano points.
pvalue_cutoff <- 0.05
volcano_stat_col <- "P.Value"       ## "P.Value" or "adj.P.Val"
volcano_y_label <- "P.Value"        ## display label on y-axis

## Focus metabolite classes for volcano highlighting and labeling.
## Non-breaking spaces in copied class names are normalized automatically.
focus_classes <- c(
  "Phenolamides",
  "Benzoic acid derivatives",
  "Hydroxycinnamoyl derivatives",
  "Flavone",
  "Flavonol",
  "Flavanone",
  "Flavone C-glycosides",
  "BX"
)

## BX class naming can vary among annotation files.
focus_bx_patterns <- c("^BX$", "Benzoxazinoid", "Benzoxazinoids")

## Volcano label logic
volcano_label_focus_only <- TRUE
volcano_label_sig_focus_only <- TRUE
volcano_label_all_focus <- TRUE
volcano_label_max <- Inf


## Automatically generated comparisons
comparison_times <- c("24h", "48h")
case_treatment <- "Fg"
control_treatment <- "CK"

## OPLS-DA
run_oplsda <- TRUE
opls_permI <- 200
opls_top_vip <- 30
opls_fallback_vip <- TRUE


## Focus-class boxplot settings
boxplot_sig_stat_col <- "P.Value"
boxplot_pvalue_cutoff <- pvalue_cutoff
boxplot_focus_only <- TRUE
boxplot_show_all_sig_focus <- TRUE
boxplot_ncol <- 4
boxplot_max_facets <- Inf
boxplot_min_height <- 7.5
boxplot_row_height <- 2.0

## VIP plot settings
vip_focus_only <- TRUE
vip_sig_only <- TRUE
vip_add_to_volcano <- TRUE
vip_panel_top_n <- 15
vip_panel_width_ratio <- 0.42
vip_standalone_top_n <- 20

## =========================================================
## Publication-grade plotting
## =========================================================
plot_dpi <- 600
base_font_family <- "Arial"
base_font_size <- 12

## Figure sizes
pca_width <- 10
pca_height <- 6
volcano_width <- 6.2
volcano_height <- 5.4
boxplot_width <- 11
boxplot_height <- 7.5
vip_width <- 7
vip_height <- 8
heatmap_width <- 16
heatmap_height <- 10

## Heatmap
## v5 主热图逻辑：
## - 使用总数据矩阵，不再只画每个比较的 top 差异代谢物
## - 行 = 所有代谢物
## - 列 = 所有样本
## - 列顺序 = CK / Fg 分组内按照 0h -> 24h -> 48h，再按品种和重复排序
## - 行按 Class 注释并分块
make_total_heatmap <- TRUE
heatmap_use_matrix <- "scaled"     ## "raw", "normalized", "log", or "scaled"; scaled = Pareto-scaled matrix for heatmap
heatmap_show_all_metabolites <- TRUE
heatmap_column_split <- "treatment" ## 当前按 CK / Fg 分组
heatmap_row_split <- "Class"
heatmap_cluster_rows <- TRUE
heatmap_cluster_columns <- FALSE
heatmap_show_metabolite_names <- FALSE
heatmap_show_sample_names <- FALSE
heatmap_clip_value <- 3            ## for scaled heatmap, clip around 0 with a symmetric color scale
heatmap_color_quantiles <- c(0.05, 0.50, 0.95)  ## not used when heatmap_clip_value is finite

## 热图尺寸：如果代谢物很多，建议自动按行数增加 PDF/PNG 高度
heatmap_width <- 16
heatmap_height <- 10
heatmap_auto_height <- TRUE
heatmap_row_height_in <- 0.11
heatmap_max_height <- 60

## 兼容旧流程：是否继续输出每个比较的 heatmap
## 现在建议 FALSE，只输出一张全量 ComplexHeatmap 主图
make_comparison_heatmaps <- FALSE
heatmap_show_all_groups <- TRUE
heatmap_top_n <- 50
make_global_heatmap <- FALSE
heatmap_global_top_n <- 100
heatmap_show_row_names_limit <- Inf

## Other plot limits
boxplot_top_n <- 12
volcano_label_top_n <- 5

## KEGG annotation file, optional.
## Required columns if used:
## Compound name    KEGG
kegg_annotation_file <- file.path(project_dir, "data/annotation/kegg_annotation.xlsx")
kegg_id_col <- "KEGG"

## Output directories
dir_results <- file.path(project_dir, "results")
dir_rds     <- file.path(dir_results, "rds")
dir_tables  <- file.path(dir_results, "tables")
dir_pca     <- file.path(dir_results, "pca")
dir_oplsda  <- file.path(dir_results, "oplsda")
dir_limma   <- file.path(dir_results, "limma")
dir_volcano <- file.path(dir_results, "volcano")
dir_heatmap <- file.path(dir_results, "heatmap")
dir_boxplot <- file.path(dir_results, "boxplot")
dir_kegg    <- file.path(dir_results, "kegg")

dirs <- c(
  dir_results, dir_rds, dir_tables,
  dir_pca, dir_oplsda, dir_limma,
  dir_volcano, dir_heatmap, dir_boxplot, dir_kegg
)

invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
