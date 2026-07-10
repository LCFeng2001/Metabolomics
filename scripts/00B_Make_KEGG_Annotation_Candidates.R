## =========================================================
## 00B_Make_KEGG_Annotation_Candidates.R
## Search KEGG candidates by Compound name
## Output must be manually checked.
## =========================================================

source("config/config.R")
source("functions/utils/common_utils.R")

suppressPackageStartupMessages({
  library(KEGGREST)
  library(openxlsx)
  library(dplyr)
})

data01_file <- file.path(dir_rds, "01_clean_aggregated.rds")

if (!file.exists(data01_file)) {
  stop("Please run scripts/01_ReadCleanAggregate.R first.")
}

data01 <- readRDS(data01_file)
anno <- as.data.frame(data01$annotation, check.names = FALSE)

if (!"Compound name" %in% colnames(anno)) {
  stop("annotation table lacks Compound name.")
}

compound_names <- unique(anno[["Compound name"]])
compound_names <- compound_names[!is.na(compound_names) & compound_names != ""]

find_kegg_one <- function(x) {
  message("Searching KEGG: ", x)

  ans <- tryCatch(
    KEGGREST::keggFind("compound", x),
    error = function(e) NULL
  )

  if (is.null(ans) || length(ans) == 0) {
    return(data.frame(
      `Compound name` = x,
      KEGG = NA_character_,
      KEGG_name = NA_character_,
      match_rank = NA_integer_,
      use_as_final = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  data.frame(
    `Compound name` = x,
    KEGG = gsub("^cpd:", "", names(ans)),
    KEGG_name = as.character(ans),
    match_rank = seq_along(ans),
    use_as_final = ifelse(seq_along(ans) == 1, "CHECK", ""),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

kegg_candidates <- bind_rows(lapply(compound_names, find_kegg_one))

template <- anno[, intersect(c("metabolite_id", "Compound name", "Class", "source_features"), colnames(anno)), drop = FALSE]
template[[kegg_id_col]] <- NA_character_
template$note <- "Fill KEGG Compound ID such as C00065. Use candidate sheet for reference."

outfile <- file.path(dir_tables, "00B_KEGG_annotation_candidates.xlsx")

openxlsx::write.xlsx(
  list(
    KEGG_annotation_template = template,
    KEGG_candidates = kegg_candidates
  ),
  file = outfile,
  overwrite = TRUE
)

message("KEGG candidate annotation table generated:")
message(outfile)
message("After manual checking, save final file as:")
message(kegg_annotation_file)
