## =========================================================
## 03_PCA.R
## Publication-level PCA
## =========================================================

source("config/config.R")
source("functions/utils/common_utils.R")
source("functions/plots/save_plot.R")
source("functions/plots/plot_pca.R")

data02 <- readRDS(file.path(dir_rds, "02_preprocessed.rds"))

pca_res <- run_pca_analysis(
  expr_scaled = data02$expr_scaled,
  sample_info = data02$sample_info
)

p <- plot_pca(
  pca_res,
  base_size = base_font_size,
  base_family = base_font_family
)

save_ggplot_multi(
  p = p,
  filename_prefix = file.path(dir_pca, "PCA_all_samples_publication"),
  width = pca_width,
  height = pca_height,
  dpi = plot_dpi
)

saveRDS(pca_res, file.path(dir_rds, "03_pca_result.rds"))

write_xlsx_list(
  list(
    PCA_scores = pca_res$scores,
    PCA_variance = data.frame(
      PC = paste0("PC", seq_along(pca_res$var_exp)),
      variance = pca_res$var_exp,
      variance_percent = pca_res$var_exp * 100
    )
  ),
  file.path(dir_tables, "03_PCA_output.xlsx")
)

message("03_PCA finished.")
