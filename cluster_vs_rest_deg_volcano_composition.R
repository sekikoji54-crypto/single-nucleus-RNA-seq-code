################################################################################
# Cluster-vs-rest DEG analysis (Seurat FindMarkers) + volcano plots
# and cluster composition across 4 experimental groups (100% stacked bar plot)
#
# - No hard-coded absolute/local paths
# - Input/output are provided via command-line arguments
# - Designed for PUBLIC GitHub (anonymized)
#
# Usage:
#   Rscript cluster_deg_volcano_public_safe.R \
#     --rds path/to/seurat_object.rds \
#     --out results/cluster_deg_volcano
#
# Notes:
# - Cluster identity: uses Idents(obj) if set; otherwise falls back to 'seurat_clusters'
# - Group column: auto-detected from common candidates; see group_col_candidates
# - Volcano thresholds:
#     padj < 0.05 and |avg_log2FC| >= 0.25
# - Exports:
#     * Cluster<id>_vs_rest_DEG.csv
#     * Cluster<id>_volcano.tiff
#     * Cluster_composition_100pct.tiff
#     * Cluster_composition_table.csv
#     * sessionInfo.txt
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
  if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
  if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
  if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggrepel)
  library(tidyr)
  library(scales)
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

rds_path <- get_arg("--rds")
out_dir  <- get_arg("--out", default = file.path("results", "cluster_deg_volcano"))

if (is.null(rds_path)) stop("Missing required argument: --rds <path_to_seurat_rds>")
if (!file.exists(rds_path)) stop("Input RDS file not found (path not printed).")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Running cluster DEG + volcano pipeline (paths are not printed).")
message("Outputs will be written under the specified output directory.")

#-----------------------------#
# 1) User-configurable params #
#-----------------------------#
group_col_candidates <- c("group", "Group", "condition", "Condition", "Group4", "genotype_surgery")

# Group order and colors (edit if your project uses different labels)
group_order <- c("WT-sham", "WT-axotomy", "KO-sham", "KO-axotomy")
group_colors <- c(
  "WT-sham"    = "#88CCEE",
  "WT-axotomy" = "#4477AA",
  "KO-sham"    = "#CCBB44",
  "KO-axotomy" = "#AA7733"
)

# DEG thresholds
padj_cut <- 0.05
lfc_cut  <- 0.25  # corresponds to Seurat's avg_log2FC

# Plot theme
theme_set(theme_bw(base_family = "Arial"))

#-----------------------------#
# 2) Load object + identities #
#-----------------------------#
obj <- readRDS(rds_path)
if (!inherits(obj, "Seurat")) stop("Input object is not a Seurat object.")

# Ensure cluster identities
if (is.null(Idents(obj))) {
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    Idents(obj) <- "seurat_clusters"
  } else {
    stop("Cluster identity not found: neither Idents(obj) nor 'seurat_clusters' exists.")
  }
}

# Detect group column
grp_col <- NULL
for (cc in group_col_candidates) {
  if (cc %in% colnames(obj@meta.data)) { grp_col <- cc; break }
}
if (is.null(grp_col)) {
  stop(
    paste0(
      "No group column found. Expected one of: ",
      paste(group_col_candidates, collapse = ", ")
    )
  )
}

# If labels match, enforce intended order
if (all(group_order %in% obj@meta.data[[grp_col]])) {
  obj@meta.data[[grp_col]] <- factor(obj@meta.data[[grp_col]], levels = group_order)
}

#-----------------------------#
# 3) Helper functions         #
#-----------------------------#
classify_status <- function(lfc, padj, lfc_cut = 0.25, padj_cut = 0.05) {
  if (is.na(padj)) return("NS")
  if (padj < padj_cut && lfc >=  lfc_cut) return("UP")
  if (padj < padj_cut && lfc <= -lfc_cut) return("DOWN")
  "NS"
}

top_labels <- function(df, top_n = 30, padj_cut = 0.05) {
  df %>%
    filter(p_val_adj < padj_cut) %>%
    arrange(p_val_adj) %>%
    slice_head(n = top_n) %>%
    pull(gene)
}

#-----------------------------#
# 4) Cluster-vs-rest DEG      #
#-----------------------------#
all_clusters <- sort(unique(Idents(obj)))
message("Number of clusters: ", length(all_clusters))

for (cl in all_clusters) {
  message("Processing cluster ", cl, " ...")

  cells_in  <- WhichCells(obj, idents = cl)
  cells_out <- setdiff(colnames(obj), cells_in)

  # FindMarkers (Wilcoxon, two-sided)
  deg <- FindMarkers(
    object = obj,
    ident.1 = cells_in,
    ident.2 = cells_out,
    test.use = "wilcox",
    logfc.threshold = 0,   # keep all genes; thresholding is done afterwards
    min.pct = 0.05,
    only.pos = FALSE
  )

  deg <- as.data.frame(deg) %>%
    rownames_to_column("gene")

  # Compatibility: Seurat may output avg_logFC instead of avg_log2FC
  if ("avg_logFC" %in% names(deg) && !("avg_log2FC" %in% names(deg))) {
    deg <- deg %>% rename(avg_log2FC = avg_logFC)
  }

  deg <- deg %>%
    mutate(
      status = mapply(
        classify_status, avg_log2FC, p_val_adj,
        MoreArgs = list(lfc_cut = lfc_cut, padj_cut = padj_cut)
      ),
      neglog10_padj = -log10(p_val_adj)
    )

  # Save DEG table
  csv_path <- file.path(out_dir, paste0("Cluster", cl, "_vs_rest_DEG.csv"))
  write_csv(deg, csv_path)

  # Volcano plot (label top 30 significant genes by adjusted p-value)
  lab_genes <- top_labels(deg, top_n = 30, padj_cut = padj_cut)

  p_vol <- ggplot(deg, aes(x = avg_log2FC, y = neglog10_padj)) +
    geom_point(aes(color = status), size = 1.1, alpha = 0.8) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = subset(deg, gene %in% lab_genes),
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.size = 0.2
    ) +
    scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NS" = "grey70")) +
    labs(
      title = paste0("Cluster ", cl, " vs Rest: Volcano plot"),
      x = "log2 fold-change (avg_log2FC)",
      y = "-log10(adjusted p-value)",
      color = NULL
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )

  # Save volcano (TIFF only)
  tiff_path <- file.path(out_dir, paste0("Cluster", cl, "_volcano.tiff"))
  ggsave(
    filename = tiff_path,
    plot = p_vol,
    width = 6, height = 5,
    dpi = 400,
    device = "tiff",
    compression = "lzw"
  )
}

#-----------------------------#
# 5) Cluster composition plot #
#-----------------------------#
meta <- obj@meta.data %>%
  mutate(cluster = as.factor(Idents(obj)))

# Enforce group order if labels match; otherwise keep as-is with a warning.
if (all(group_order %in% meta[[grp_col]])) {
  meta[[grp_col]] <- factor(meta[[grp_col]], levels = group_order)
} else {
  warning("Group labels do not perfectly match the expected group_order; colors/order may differ.")
}

comp <- meta %>%
  count(cluster, !!sym(grp_col), name = "n") %>%
  group_by(cluster) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

p_prop <- ggplot(comp, aes(x = cluster, y = freq, fill = !!sym(grp_col))) +
  geom_col(width = 0.75, color = "white") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(
    title = "Cluster-wise composition across experimental groups",
    x = "Cluster",
    y = "Proportion",
    fill = "Group"
  ) +
  theme(
    axis.text.x = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

prop_tiff <- file.path(out_dir, "Cluster_composition_100pct.tiff")
ggsave(
  filename = prop_tiff,
  plot = p_prop,
  width = 7, height = 5,
  dpi = 400,
  device = "tiff",
  compression = "lzw"
)

# Save numeric composition table
comp_wide <- comp %>%
  mutate(Group = as.character(!!sym(grp_col))) %>%
  select(cluster, Group, freq) %>%
  pivot_wider(names_from = Group, values_from = freq)

write_csv(comp_wide, file.path(out_dir, "Cluster_composition_table.csv"))

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Done.")
