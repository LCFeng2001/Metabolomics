## =========================================================
## kegg_enrichment.R
## Optional KEGG compound enrichment
## =========================================================

run_kegg_compound_enrichment <- function(sig_kegg, universe_kegg) {
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    warning("KEGGREST is not installed. KEGG skipped.")
    return(NULL)
  }

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    warning("clusterProfiler is not installed. KEGG skipped.")
    return(NULL)
  }

  sig_kegg <- unique(na.omit(as.character(sig_kegg)))
  universe_kegg <- unique(na.omit(as.character(universe_kegg)))

  sig_kegg <- gsub("^cpd:", "", sig_kegg)
  universe_kegg <- gsub("^cpd:", "", universe_kegg)

  if (length(sig_kegg) < 3) {
    warning("Too few KEGG compounds for enrichment.")
    return(NULL)
  }

  link <- KEGGREST::keggLink("pathway", "cpd")

  term2gene <- data.frame(
    pathway = gsub("^path:", "", unname(link)),
    KEGG = gsub("^cpd:", "", names(link)),
    stringsAsFactors = FALSE
  )

  pathway_names <- KEGGREST::keggList("pathway")

  term2name <- data.frame(
    pathway = gsub("^path:", "", names(pathway_names)),
    name = as.character(pathway_names),
    stringsAsFactors = FALSE
  )

  clusterProfiler::enricher(
    gene = sig_kegg,
    universe = universe_kegg,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
}
