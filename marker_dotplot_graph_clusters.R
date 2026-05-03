################################################################################
# Marker-based dot plot for graph-based clusters
#
# - No absolute/local paths are hard-coded or printed
# - Input/output locations are provided via command-line arguments
# - Cluster labels remain unlabeled (numeric IDs only)
#
# Usage:
#   Rscript dotplot_markers_public_safe.R \
#     --h5   path/to/raw_feature_bc_matrix.h5 \
#     --csv  path/to/graph_based_clusters.csv \
#     --out  results/dotplot
#
# Required CSV columns:
#   - Barcode
#   - Graph-based   (e.g., "Cluster 23" or any string containing an integer cluster ID)
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

  library(Seurat)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
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

h5_path  <- get_arg("--h5")
csv_path <- get_arg("--csv")
out_dir  <- get_arg("--out", default = file.path("results", "dotplot"))

if (is.null(h5_path))  stop("Missing required argument: --h5 <path_to_10x_h5>")
if (is.null(csv_path)) stop("Missing required argument: --csv <path_to_cluster_csv>")

if (!file.exists(h5_path))  stop("10x H5 file not found (path not printed).")
if (!file.exists(csv_path)) stop("Cluster CSV file not found (path not printed).")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Running marker-based dot plot pipeline (paths are not printed).")
message("Outputs will be written under the specified output directory.")

#-----------------------------#
# 1) Load snRNA-seq counts    #
#-----------------------------#
obj <- CreateSeuratObject(counts = Read10X_h5(h5_path))

# Basic preprocessing for visualization
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)

#-----------------------------#
# 2) Load cluster assignments #
#-----------------------------#
cldf_raw <- readr::read_csv(csv_path, show_col_types = FALSE)
required_cols <- c("Barcode", "Graph-based")
if (!all(required_cols %in% names(cldf_raw))) {
  stop("CSV must contain columns: 'Barcode' and 'Graph-based'.")
}

df <- cldf_raw %>%
  transmute(
    barcode = as.character(Barcode),
    cl_num  = as.integer(stringr::str_extract(`Graph-based`, "\\d+"))
  ) %>%
  filter(!is.na(cl_num))

# Map barcode -> cluster ID and attach to Seurat object
cl_map <- setNames(df$cl_num, df$barcode)
obj$GraphClust <- cl_map[colnames(obj)]

# Keep only cells with a matched cluster assignment
obj <- subset(obj, cells = colnames(obj)[!is.na(obj$GraphClust)])

obj$GraphClust <- factor(obj$GraphClust)
Idents(obj)    <- obj$GraphClust

# NOTE (public-safe): Do not print tables that might indirectly reveal dataset specifics.
# If you need QC locally, you can uncomment the next line:
# print(table(obj$GraphClust))

#-----------------------------#
# 3) Marker sets              #
#-----------------------------#
markers_neuron <- c(
  "Gad1", "Gad2",                         # GABAergic
  "Slc17a6", "Slc17a7", "Slc17a8",        # Glutamatergic
  "Slc6a5", "Slc7a8",                     # Glycinergic
  "Chat", "Slc5a7", "Ache",               # Cholinergic
  "Slc6a4", "Fev",                        # Serotonergic
  "Th", "Slc18a2",                        # Dopaminergic
  "Slc6a2",                               # Noradrenergic
  "Slc16a2"                               # Auxiliary marker (optional)
)
# Optional helper for low nuclear signal:
# markers_neuron <- union(markers_neuron, "Ddc")

markers_non_neuron <- c(
  "Aldh1a1",            # Astrocyte
  "Ctss", "Tmem119",    # Microglia
  "Mbp",                # Oligodendrocyte
  "Vcan",               # OPC
  "Ccdc170",            # Ependymal
  "Cldn5",              # Endothelial
  "Pdgfrb",             # Pericyte
  "Col8a1",             # Fibroblast
  "Ttr"                 # Choroid plexus
)

markers <- c(markers_neuron, markers_non_neuron)

features_use <- intersect(markers, rownames(obj))
if (length(features_use) == 0) stop("None of the specified marker genes were found in the dataset.")

#-----------------------------#
# 4) DotPlot data -> ggplot   #
#-----------------------------#
dp <- DotPlot(obj, features = features_use, group.by = "GraphClust")$data

# Remove non-finite values to avoid warnings and keep output clean
dp <- dp[is.finite(dp$avg.exp.scaled) & is.finite(dp$pct.exp), , drop = FALSE]

# Preserve intended ordering
dp$features.plot <- factor(dp$features.plot, levels = rev(features_use))
dp$id            <- factor(dp$id, levels = levels(obj$GraphClust))

p <- ggplot(dp, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  coord_flip() +
  scale_size_area(max_size = 6) +
  scale_y_discrete(drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_color_gradient(low = "navy", high = "goldenrod") +
  theme_bw(base_size = 10, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold")
  ) +
  labs(
    title   = "Marker-based dot plot (graph-based clusters, unlabeled)",
    x       = "Graph-based cluster ID",
    y       = "Marker genes",
    caption = "Dot size: % expressing; Color: scaled mean expression (z-score)."
  ) +
  geom_hline(
    yintercept = length(markers_neuron) + 0.5,
    linetype = "dashed",
    linewidth = 0.4
  )

#-----------------------------#
# 5) Save (journal-friendly)  #
#-----------------------------#
w_in <- 180 / 25.4  # 180 mm
h_in <- 120 / 25.4  # 120 mm (adjust as needed)

ggsave(
  filename = file.path(out_dir, "DotPlot_unlabeled.pdf"),
  plot = p, width = w_in, height = h_in, dpi = 300, useDingbats = FALSE
)

ggsave(
  filename = file.path(out_dir, "DotPlot_unlabeled.tiff"),
  plot = p, width = w_in, height = h_in, dpi = 300, compression = "lzw"
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Done.")
