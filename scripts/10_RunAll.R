## =========================================================
## 10_RunAll.R
## Complete rebuilt targeted metabolomics pipeline
## =========================================================

scripts <- c(
  "scripts/01_ReadCleanAggregate.R",
  "scripts/02_Preprocess.R",
  "scripts/03_PCA.R",
  "scripts/04_Comparisons.R",
  "scripts/05_Limma.R",
  "scripts/06_OPLSDA.R",
  "scripts/07_Plots.R",
  "scripts/08_KEGG.R"
)

for (s in scripts) {
  message("\n======================================")
  message("Running: ", s)
  message("======================================\n")

  tryCatch({
    source(s)
  }, error = function(e) {
    stop("Pipeline stopped at ", s, "\nError: ", e$message)
  })
}

message("\nAll analysis finished.")
