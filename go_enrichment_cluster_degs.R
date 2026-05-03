################################################################################
# GO enrichment analysis (mouse; clusterProfiler) for Up/Down DEGs 
#
# - No hard-coded absolute/local paths
# - Input/output provided via command-line arguments
# - Mouse annotation: org.Mm.eg.db (SYMBOL -> ENTREZID)
#
# Usage:
#   Rscript go_enrichment_mouse_public_safe.R \
#     --in   path/to/deg_table.csv \
#     --out  results/go_enrichment \
#     --prefix WT-axo_vs_KO-axo_ChAT-posi
#
# Input table requirements (default column names):
#   - Gene       (mouse gene symbol)
#   - avg_log2FC (log2 fold-change)
#   - p_val_adj  (adjusted p-value)
#
# Optional overrides:
#   --gene_col "Gene" --fc_col "avg_log2FC" --padj_col "p_val_adj"
################################################################################

suppressPackageStartupMessages({
  # CRAN
  pkgs_cran <- c("dplyr", "ggplot2", "patchwork", "stringr", "readr")
  to_install <- pkgs_cran[!vapply(pkgs_cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

  # Bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  pkgs_bioc <- c("clusterProfiler", "AnnotationDbi", "org.Mm.eg.db")
  miss_bioc <- pkgs_bioc[!vapply(pkgs_bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_bioc) > 0) BiocManager::install(miss_bioc, ask = FALSE, update = FALSE)

  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(readr)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

set.seed(1)

#-----------------------------#
# 0) Argument parsing         #
#-----------------------------#
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 1 && length(args) >= hit + 1) return(args[hit + 1])
  default
}

in_file <- get_arg("--in")
out_dir <- get_arg("--out", default = file.path("results", "go_enrichment"))
prefix  <- get_arg("--prefix", default = "DEG_GO")

gene_col <- get_arg("--gene_col", default = "Gene")
fc_col   <- get_arg("--fc_col",   default = "avg_log2FC")
padj_col <- get_arg("--padj_col", default = "p_val_adj")

padj_cut  <- as.numeric(get_arg("--padj_cut", default = "0.05"))
show_n    <- as.integer(get_arg("--show_n",   default = "10"))
show_cnet <- as.integer(get_arg("--show_cnet", default = "5"))

if (is.null(in_file)) stop("Missing required argument: --in <path_to_deg_csv>")
if (!file.exists(in_file)) stop("Input DEG file not found (path not printed).")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Running GO enrichment pipeline (mouse; paths are not printed).")
message("Outputs will be written under the specified output directory.")

#-----------------------------#
# 1) Load DEG table           #
#-----------------------------#
deg <- readr::read_csv(in_file, show_col_types = FALSE)

req <- c(gene_col, fc_col, padj_col)
if (!all(req %in% names(deg))) {
  stop(
    paste0(
      "Required columns are missing.\n",
      "Expected: ", paste(req, collapse = ", "), "\n",
      "Found columns include: ", paste(head(names(deg), 25), collapse = ", "),
      ifelse(ncol(deg) > 25, ", ...", "")
    )
  )
}

# Standardize columns internally (do not modify original names)
data <- deg %>%
  transmute(
    Gene       = toupper(as.character(.data[[gene_col]])),
    avg_log2FC = as.numeric(.data[[fc_col]]),
    p_val_adj  = as.numeric(.data[[padj_col]])
  ) %>%
  filter(!is.na(Gene), Gene != "", is.finite(avg_log2FC), is.finite(p_val_adj))

#-----------------------------#
# 2) SYMBOL -> ENTREZ         #
#-----------------------------#
convert_to_entrez <- function(gene_symbols) {
  gene_symbols <- unique(gene_symbols)
  ids <- AnnotationDbi::mapIds(
    x         = org.Mm.eg.db,
    keys      = gene_symbols,
    column    = "ENTREZID",
    keytype   = "SYMBOL",
    multiVals = "first"
  )
  ids <- as.character(ids)
  ids <- ids[!is.na(ids)]
  unique(ids)
}

# Universe/background genes
universe_ids <- convert_to_entrez(data$Gene)
if (length(universe_ids) == 0) stop("No Entrez IDs could be mapped from the input gene symbols.")

# Split Up/Down (padj < threshold)
up_df <- data %>% filter(avg_log2FC > 0, p_val_adj < padj_cut)
dn_df <- data %>% filter(avg_log2FC < 0, p_val_adj < padj_cut)

up_ids <- convert_to_entrez(up_df$Gene)
dn_ids <- convert_to_entrez(dn_df$Gene)

# Fold-change vector for cnetplot (named by gene symbol)
geneList <- setNames(data$avg_log2FC, data$Gene)

#-----------------------------#
# 3) GO enrichment (BP/MF/CC) #
#-----------------------------#
perform_go <- function(entrez_ids, universe_ids) {
  if (length(entrez_ids) == 0) return(list(BP = NULL, MF = NULL, CC = NULL))

  list(
    BP = enrichGO(
      gene          = entrez_ids,
      universe      = universe_ids,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = padj_cut,
      qvalueCutoff  = padj_cut,
      readable      = TRUE
    ),
    MF = enrichGO(
      gene          = entrez_ids,
      universe      = universe_ids,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff  = padj_cut,
      qvalueCutoff  = padj_cut,
      readable      = TRUE
    ),
    CC = enrichGO(
      gene          = entrez_ids,
      universe      = universe_ids,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff  = padj_cut,
      qvalueCutoff  = padj_cut,
      readable      = TRUE
    )
  )
}

up_go <- perform_go(up_ids, universe_ids)
dn_go <- perform_go(dn_ids, universe_ids)

#-----------------------------#
# 4) Plot + save              #
#-----------------------------#
plot_save_one <- function(go_obj, ontology, direction_label) {
  if (is.null(go_obj) || nrow(as.data.frame(go_obj)) == 0) {
    message("[SKIP] ", direction_label, " / ", ontology, ": enrichment result is empty.")
    return(invisible(NULL))
  }

  g1 <- barplot(go_obj, drop = TRUE, showCategory = show_n)
  g2 <- clusterProfiler::dotplot(go_obj, showCategory = show_n)

  g3 <- try(
    clusterProfiler::cnetplot(
      go_obj,
      showCategory = show_cnet,
      color.params = list(foldChange = geneList)
    ),
    silent = TRUE
  )

  if (inherits(g3, "try-error")) {
    message("[WARN] cnetplot failed for ", direction_label, " / ", ontology, " (saving barplot+dotplot).")
    g <- g1 + g2 + plot_layout(ncol = 2)
    w <- 18; h <- 6
  } else {
    g <- g1 + g2 + g3 + plot_layout(ncol = 3)
    w <- 25; h <- 10
  }

  out_file <- file.path(out_dir, sprintf("%s_%s_%s.tiff", prefix, ontology, direction_label))

  ggsave(
    filename    = out_file,
    plot        = g,
    width       = w,
    height      = h,
    dpi         = 300,
    device      = "tiff",
    compression = "lzw"
  )
  message("[SAVED] ", basename(out_file))
}

plot_save_set <- function(go_list, direction_label) {
  plot_save_one(go_list$BP, "BP", direction_label)
  plot_save_one(go_list$MF, "MF", direction_label)
  plot_save_one(go_list$CC, "CC", direction_label)
}

plot_save_set(up_go, "Upregulated")
plot_save_set(dn_go, "Downregulated")

# Reproducibility (public-safe)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done.")
################################################################################
