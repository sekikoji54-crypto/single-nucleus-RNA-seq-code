################################################################################
# TailScore computation and visualization
#
# TailScore definition:
#   TailScore = mean(z(up genes)) - mean(z(down genes))
# where z is computed per gene across all cells (mean=0, sd=1) using log-normalized
# expression (Seurat "data" slot).
#
# Outputs (in --out):
#   - UMAP overlay of TailScore (TIFF)
#   - Violin + box plots of TailScore by cluster (TIFF)
#   - Violin + box plots of TailScore by group (TIFF; if group column exists)
#   - sessionInfo.txt
#
# Usage:
#   Rscript tailscore_public_safe.R \
#     --rds  path/to/seurat_object.rds \
#     --out  results/tailscore \
#     --up   path/to/tail_up_genes.txt \
#     --down path/to/tail_down_genes.txt \
#     --group_col group
#
# Notes:
#   - This script does NOT print local paths (safe for public GitHub).
#   - Ensure that the Seurat object has UMAP embeddings (reduction "umap") if you
#     want the UMAP plot. If missing, you can add UMAP upstream.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
  if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

  library(Seurat)
  library(Matrix)
  library(ggplot2)
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

rds_path   <- get_arg("--rds")
out_dir    <- get_arg("--out", default = file.path("results", "tailscore"))
up_file    <- get_arg("--up",   default = file.path(out_dir, "tail_up_genes.txt"))
down_file  <- get_arg("--down", default = file.path(out_dir, "tail_down_genes.txt"))
group_col  <- get_arg("--group_col", default = "group")

if (is.null(rds_path)) stop("Missing required argument: --rds <path_to_seurat_rds>")
if (!file.exists(rds_path)) stop("Input RDS file not found (path not printed).")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Running TailScore pipeline (paths are not printed).")
message("Outputs will be written under the specified output directory.")

#-----------------------------#
# 1) Load Seurat object       #
#-----------------------------#
obj <- readRDS(rds_path)
if (!inherits(obj, "Seurat")) stop("Input object is not a Seurat object.")

# Ensure Idents are set (for cluster-wise violin)
if (is.null(Idents(obj))) {
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    Idents(obj) <- "seurat_clusters"
  } else {
    stop("Idents(obj) is not set and 'seurat_clusters' was not found in meta.data.")
  }
}

#-----------------------------#
# 2) Read tail gene lists     #
#-----------------------------#
read_gene_list <- function(path) {
  if (!file.exists(path)) return(character(0))
  x <- readLines(path, warn = FALSE)
  x <- unique(trimws(x))
  x <- x[nzchar(x)]
  x
}

tail_up   <- read_gene_list(up_file)
tail_down <- read_gene_list(down_file)

genes_in_obj <- rownames(obj)
tail_up      <- intersect(tail_up, genes_in_obj)
tail_down    <- intersect(tail_down, genes_in_obj)

if (length(tail_up) + length(tail_down) < 3) {
  stop("Too few genes in tail_up/tail_down after intersecting with the object. Check your input lists.")
}

#-----------------------------#
# 3) Compute TailScore        #
#-----------------------------#
# Use log-normalized expression by default (Seurat 'data' slot)
E <- GetAssayData(obj, slot = "data")  # genes x cells (dgCMatrix)

# Efficient z-scoring per gene for a small gene subset (keeps sparsity)
z_subset <- function(mat_dgC, genes) {
  sub <- mat_dgC[genes, , drop = FALSE]
  m  <- Matrix::rowMeans(sub)
  m2 <- Matrix::rowMeans(sub^2)
  v  <- pmax(m2 - m^2, 0)
  s  <- sqrt(v)
  s[s == 0] <- 1

  # Center and scale (works with dgCMatrix; may densify slightly depending on ops)
  sub_centered <- sub - m
  sub_scaled   <- sub_centered / s
  sub_scaled
}

Z_up   <- if (length(tail_up)   > 0) z_subset(E, tail_up)   else NULL
Z_down <- if (length(tail_down) > 0) z_subset(E, tail_down) else NULL

pos_score <- if (!is.null(Z_up))   Matrix::colMeans(Z_up)   else rep(0, ncol(E))
neg_score <- if (!is.null(Z_down)) Matrix::colMeans(Z_down) else rep(0, ncol(E))

tail_score <- as.numeric(pos_score - neg_score)
names(tail_score) <- colnames(obj)

obj$TailScore <- tail_score

#-----------------------------#
# 4) UMAP overlay (optional)  #
#-----------------------------#
if ("umap" %in% names(obj@reductions)) {
  p_umap <- FeaturePlot(obj, features = "TailScore", reduction = "umap") +
    theme_bw(base_family = "Arial") +
    ggtitle("TailScore (mean z[up] - mean z[down])")

  ggsave(
    filename = file.path(out_dir, "UMAP_TailScore.tiff"),
    plot = p_umap,
    width = 6.8, height = 5.6, dpi = 400,
    device = "tiff", compression = "lzw"
  )
} else {
  message("UMAP reduction not found; skipping UMAP_TailScore.tiff.")
}

#-----------------------------#
# 5) Violin by cluster        #
#-----------------------------#
df_cl <- data.frame(
  TailScore = obj$TailScore,
  cluster   = as.character(Idents(obj))
)

# Natural ordering if clusters are numeric-like
lev_cl <- unique(df_cl$cluster)
suppressWarnings({
  if (all(grepl("^[0-9]+$", lev_cl))) lev_cl <- as.character(sort(as.integer(lev_cl)))
})
df_cl$cluster <- factor(df_cl$cluster, levels = lev_cl)

p_vln_cl <- ggplot(df_cl, aes(x = cluster, y = TailScore, fill = cluster)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
  labs(x = "Cluster", y = "TailScore", title = "TailScore by cluster") +
  theme_bw(base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )

ggsave(
  filename = file.path(out_dir, "TailScore_violin_byCluster.tiff"),
  plot = p_vln_cl,
  width = 7.4, height = 4.8, dpi = 400,
  device = "tiff", compression = "lzw"
)

#-----------------------------#
# 6) Violin by group (optional)
#-----------------------------#
if (group_col %in% colnames(obj@meta.data)) {

  df_gp <- data.frame(
    TailScore = obj$TailScore,
    group     = obj@meta.data[[group_col]]
  )

  # If your study uses canonical group labels, define them here (optional).
  canonical_levels <- c("WT-sham", "WT-axotomy", "KO-sham", "KO-axotomy")
  if (all(canonical_levels %in% unique(as.character(df_gp$group)))) {
    df_gp$group <- factor(as.character(df_gp$group), levels = canonical_levels)
  } else {
    df_gp$group <- factor(as.character(df_gp$group))
  }

  p_vln_gp <- ggplot(df_gp, aes(x = group, y = TailScore, fill = group)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
    labs(x = NULL, y = "TailScore", title = "TailScore by group") +
    theme_bw(base_family = "Arial") +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    filename = file.path(out_dir, "TailScore_violin_byGroup.tiff"),
    plot = p_vln_gp,
    width = 7.4, height = 4.8, dpi = 400,
    device = "tiff", compression = "lzw"
  )

} else {
  message("Specified group column not found; skipping TailScore_violin_byGroup.tiff.")
}

#-----------------------------#
# 7) Reproducibility          #
#-----------------------------#
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Done.")
