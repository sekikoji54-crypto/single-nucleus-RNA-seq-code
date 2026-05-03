################################################################################
# Microglia analysis pipeline (PUBLIC-SAFE, anonymized)
#
# What this script does:
#   1) Load 10x Genomics HDF5 (raw_feature_bc_matrix.h5)
#   2) Subset cells using a barcode→Group CSV
#   3) Preprocess: SCTransform → fallback to LogNormalize if SCT fails
#   4) PCA / UMAP
#   5) Resolution sweep → choose a clustering resolution near target_k (default ~10)
#   6) Compute module scores (HM / DAM / IRM) and heatmap
#   7) FGSEA-based coarse annotation (NES>0 only; HM/DAM/IRM/Other)
#   8) DotPlots for HM/DAM/IRM gene lists (paged)
#   9) Optional violin plots (paged)
#  10) Composition stacked bar by Group + chi-square and pairwise Fisher tests (BH)
#
# Public-safety (anonymization) features:
#   - No absolute/local paths are hard-coded
#   - Script does not print input file paths
#   - Summary files do not include file paths
#
# Usage:
#   Rscript microglia_pipeline_public_safe.R \
#     --h5        path/to/raw_feature_bc_matrix.h5 \
#     --barcodes  path/to/barcode_group.csv \
#     --out       results/microglia_k10 \
#     --target_k  10 \
#     --seed      42
#
# Expected barcode CSV:
#   - First column: barcode (e.g., "AAAC...-1")
#   - Second column: group label (e.g., WT-sham, WT-axotomy, KO-sham, KO-axotomy)
#   - If column names exist, the script tries to detect common names automatically.
################################################################################

suppressPackageStartupMessages({
  # Minimal "install if missing" helper (safe for shared scripts)
  install_if_missing <- function(pkgs) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")
  }

  cran_pkgs <- c(
    "Seurat", "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble",
    "ggplot2", "patchwork", "scales", "pheatmap", "Matrix",
    "igraph", "cluster", "ragg"
  )
  # fgsea is commonly installed from Bioconductor, but can be present already.
  # We'll attempt BiocManager install only if needed.
  install_if_missing(setdiff(cran_pkgs, "fgsea"))
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("fgsea", ask = FALSE, update = FALSE)
  }

  library(Seurat)
  library(dplyr); library(tidyr); library(purrr)
  library(readr); library(stringr); library(tibble)
  library(ggplot2); library(patchwork); library(scales)
  library(fgsea); library(pheatmap); library(Matrix)
  library(igraph); library(cluster)
  library(ragg)
})

#-----------------------------#
# 0) Argument parsing         #
#-----------------------------#
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 1 && length(args) >= hit + 1) return(args[hit + 1])
  default
}

h5_path      <- get_arg("--h5")
barc_csv     <- get_arg("--barcodes")
out_dir      <- get_arg("--out", default = file.path("results", "microglia_k10"))
target_k     <- as.integer(get_arg("--target_k", default = "10"))
seed         <- as.integer(get_arg("--seed", default = "42"))

if (is.null(h5_path))  stop("Missing required argument: --h5 <path_to_10x_h5>")
if (is.null(barc_csv)) stop("Missing required argument: --barcodes <path_to_barcode_group_csv>")
if (!file.exists(h5_path))  stop("10x H5 file not found (path not printed).")
if (!file.exists(barc_csv)) stop("Barcode CSV file not found (path not printed).")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

message("Running microglia pipeline (public-safe; file paths are not printed).")
message("Outputs will be written under the specified output directory.")

#-----------------------------#
# 1) Plot defaults            #
#-----------------------------#
theme_set(theme_bw(base_size = 11) + theme(text = element_text(family = "sans")))

ggsave_safe <- function(filename, plot, width, height, dpi = 300) {
  ext <- tolower(tools::file_ext(filename))
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)

  if (ext %in% c("tif", "tiff")) {
    ggsave(filename, plot = plot, width = width, height = height,
           dpi = dpi, device = ragg::agg_tiff, compression = "lzw")
  } else if (ext %in% c("png")) {
    ggsave(filename, plot = plot, width = width, height = height,
           dpi = dpi, device = ragg::agg_png)
  } else if (ext %in% c("pdf")) {
    ggsave(filename, plot = plot, width = width, height = height, device = cairo_pdf)
  } else {
    ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
  }
}

#-----------------------------#
# 2) Helpers                  #
#-----------------------------#
read_barcodes_groups <- function(path, expected_levels = c("WT-sham","WT-axotomy","KO-sham","KO-axotomy")) {
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))

  bar_col <- intersect(c("Barcode","barcode","Barcodes","barcodes","cell","Cell","CELL","cell_id"), colnames(df))
  grp_col <- intersect(c("Group","group","Condition","condition","cond","Cluster","cluster","Label","label"), colnames(df))

  if (length(bar_col) == 0) bar_col <- colnames(df)[1]
  if (length(grp_col) == 0) grp_col <- if (ncol(df) >= 2) colnames(df)[2] else colnames(df)[1]

  out <- df %>%
    transmute(
      Barcode = as.character(.data[[bar_col[1]]]),
      Group   = str_trim(as.character(.data[[grp_col[1]]]))
    ) %>%
    filter(!is.na(Barcode), nzchar(Barcode), !is.na(Group), nzchar(Group)) %>%
    distinct()

  # If expected levels are present, enforce their order; otherwise keep as-is
  if (all(expected_levels %in% unique(out$Group))) {
    out$Group <- factor(out$Group, levels = expected_levels)
  } else {
    out$Group <- factor(out$Group)
  }
  out
}

# Gene symbol "make.names" mapping helper (handles hyphens, etc.)
map_to_rownames <- function(cands, rn) {
  rn_norm <- setNames(rn, make.names(rn))
  ok <- intersect(make.names(cands), names(rn_norm))
  unname(rn_norm[ok])
}

plot_violin_paged <- function(obj, genes, tag, out_dir, by = "seurat_clusters", genes_per_page = 12) {
  genes <- intersect(genes, rownames(obj))
  if (!length(genes)) return(invisible(NULL))

  pages <- split(genes, ceiling(seq_along(genes) / genes_per_page))
  for (i in seq_along(pages)) {
    pg <- pages[[i]]
    plist <- Seurat::VlnPlot(obj, features = pg, group.by = by, pt.size = 0.05, combine = FALSE)
    plist <- lapply(plist, function(p) p + theme_bw() + theme(legend.position = "none"))
    pwrap <- wrap_plots(plist, ncol = 3)
    ggsave_safe(file.path(out_dir, sprintf("Violin_%s_p%02d.tiff", tag, i)), pwrap, width = 12, height = 8, dpi = 300)
    ggsave_safe(file.path(out_dir, sprintf("Violin_%s_p%02d.pdf",  tag, i)), pwrap, width = 12, height = 8)
  }
}

plot_heatmap_module_scores <- function(avg_mat, out_dir) {
  pheatmap::pheatmap(
    mat = avg_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
    fontsize = 10,
    filename = file.path(out_dir, "ModuleScores_byCluster_heatmap.pdf"),
    width = 6, height = 5
  )
}

pick_snn_name <- function(obj) {
  gnames <- names(obj@graphs)
  pick <- gnames[grepl("snn$", gnames)]
  pick <- pick[order(grepl("^SCT", pick), decreasing = TRUE)]
  if (length(pick) == 0) stop("SNN graph not found in obj@graphs.")
  pick[1]
}

evaluate_resolutions <- function(obj,
                                 res_grid = seq(0.2, 2.0, by = 0.05),
                                 dims = 30,
                                 min_pct_marker = 0.10,
                                 logfc_thr = 0.25,
                                 tiny_frac = 0.005,
                                 out_dir = ".",
                                 prefer_k = 10) {
  obj <- FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
  snn_name <- pick_snn_name(obj)
  emb <- Embeddings(obj, "pca")[, 1:dims, drop = FALSE]

  res_tbl <- lapply(res_grid, function(res) {
    tmp <- FindClusters(obj, resolution = res, verbose = FALSE)
    cl  <- tmp$seurat_clusters
    k   <- nlevels(cl)

    tab <- table(cl)
    n_cells <- length(cl)
    tiny_rate <- sum(tab < (tiny_frac * n_cells)) / max(1, k)

    sil <- NA_real_
    if (k > 1) {
      d <- dist(emb, method = "euclidean")
      sil <- tryCatch(mean(cluster::silhouette(as.integer(cl), d)[, "sil_width"]),
                      error = function(e) NA_real_)
    }

    mod <- NA_real_
    ig <- tryCatch(
      igraph::graph_from_adjacency_matrix(tmp@graphs[[snn_name]], mode = "undirected", weighted = TRUE),
      error = function(e) NULL
    )
    if (!is.null(ig) && k > 1) {
      mem <- as.integer(cl)
      mod <- tryCatch(igraph::modularity(ig, membership = mem, weights = E(ig)$weight),
                      error = function(e) NA_real_)
    }

    mrk_n_med <- NA_real_
    mrk_try <- tryCatch({
      mrk <- FindAllMarkers(tmp, only.pos = TRUE, test.use = "wilcox",
                            min.pct = min_pct_marker, logfc.threshold = logfc_thr)
      mrk %>%
        filter(p_val_adj < 0.05) %>%
        count(cluster, name = "n_sig") %>%
        summarise(med = ifelse(n() > 0, median(n_sig), 0)) %>%
        pull(med)
    }, error = function(e) NA_real_)
    if (is.numeric(mrk_try)) mrk_n_med <- mrk_try

    data.frame(
      resolution = res,
      n_clusters = k,
      silhouette = sil,
      modularity = mod,
      marker_median = mrk_n_med,
      tiny_rate = tiny_rate
    )
  }) %>% bind_rows()

  scale01 <- function(x) {
    if (all(is.na(x))) return(x)
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(x)))
    (x - rng[1]) / diff(rng)
  }

  res_tbl <- res_tbl %>%
    mutate(
      sil_n = scale01(silhouette),
      mod_n = scale01(modularity),
      mrk_n = scale01(marker_median),
      tiny_pen = 1 - pmin(1, tiny_rate * 5)
    ) %>%
    mutate(score = rowMeans(cbind(sil_n, mod_n, mrk_n, tiny_pen), na.rm = TRUE))

  readr::write_csv(res_tbl, file.path(out_dir, sprintf("resolution_sweep_metrics_targetK%d.csv", prefer_k)))

  p <- ggplot(res_tbl, aes(x = resolution)) +
    geom_line(aes(y = score), linewidth = 1.1) +
    geom_point(aes(y = score)) +
    geom_line(aes(y = sil_n), linetype = 2, alpha = .7) +
    geom_line(aes(y = mod_n), linetype = 3, alpha = .7) +
    geom_line(aes(y = mrk_n), linetype = 4, alpha = .7) +
    geom_line(aes(y = 1 - tiny_rate), linetype = 1, alpha = .4) +
    labs(
      y = "score (scaled metrics)",
      title = sprintf("Resolution sweep diagnostics (target ≈ %d clusters)", prefer_k),
      subtitle = "solid=overall; dashed=Silhouette; dotted=Modularity; dotdash=Markers; thin=1−tinyRate"
    ) +
    theme_bw()

  ggsave_safe(file.path(out_dir, sprintf("resolution_sweep_diagnostics_targetK%d.pdf", prefer_k)),
             p, width = 7.5, height = 5.5)

  cand <- res_tbl %>% filter(score >= (max(score, na.rm = TRUE) - 0.03))
  cand2 <- cand %>%
    mutate(pref = abs(n_clusters - prefer_k)) %>%
    arrange(pref, desc(silhouette), desc(modularity), desc(marker_median), tiny_rate, resolution)

  picked <- cand2 %>% slice(1)
  list(metrics = res_tbl, picked = picked)
}

dot_paged <- function(obj, genes, title_prefix, out_dir, per_page = 25) {
  genes <- intersect(genes, rownames(obj))
  if (!length(genes)) return(invisible(NULL))

  pages <- split(genes, ceiling(seq_along(genes) / per_page))
  for (i in seq_along(pages)) {
    feats <- pages[[i]]
    p <- Seurat::DotPlot(obj, features = feats, group.by = "seurat_clusters", scale = TRUE) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(sprintf("%s DotPlot (gene-wise z), page %02d", title_prefix, i))

    safe_tag <- gsub("[^A-Za-z0-9_\\-]+", "_", title_prefix)
    ggsave_safe(file.path(out_dir, sprintf("DotPlot_%s_p%02d.tiff", safe_tag, i)),
                p, width = 11, height = 6.5, dpi = 300)
    ggsave_safe(file.path(out_dir, sprintf("DotPlot_%s_p%02d.pdf", safe_tag, i)),
                p, width = 11, height = 6.5)
  }
}

#-----------------------------#
# 3) Load data + subset cells #
#-----------------------------#
mat <- Read10X_h5(h5_path)
obj <- CreateSeuratObject(counts = mat, project = "SeuratProject")

barc_tbl <- read_barcodes_groups(barc_csv)
valid_barcodes <- intersect(barc_tbl$Barcode, colnames(obj))
if (length(valid_barcodes) == 0) stop("No overlapping barcodes between H5 matrix and barcode CSV.")

obj <- subset(obj, cells = valid_barcodes)
obj$Group <- barc_tbl$Group[match(colnames(obj), barc_tbl$Barcode)]
obj$Group <- factor(obj$Group, levels = levels(barc_tbl$Group))

#-----------------------------#
# 4) Preprocess + PCA + UMAP  #
#-----------------------------#
obj_proc <- try({
  obj %>%
    SCTransform(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
}, silent = TRUE)

if (inherits(obj_proc, "try-error")) {
  obj_proc <- obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
}

#-----------------------------#
# 5) Resolution sweep + k~K   #
#-----------------------------#
eval <- evaluate_resolutions(
  obj_proc,
  res_grid = seq(0.2, 2.0, by = 0.05),
  dims = 30,
  out_dir = out_dir,
  prefer_k = target_k
)

best_res <- eval$picked$resolution

# Ensure neighbors / SNN graph exists, then cluster on a selected SNN graph
if (!any(grepl("snn$", names(obj_proc@graphs)))) {
  obj_proc <- FindNeighbors(obj_proc, dims = 1:30, verbose = FALSE)
}
snn_candidates <- names(obj_proc@graphs)[grepl("snn$", names(obj_proc@graphs))]
if (length(snn_candidates) == 0) stop("SNN graph not found after FindNeighbors().")
snn_name_use <- snn_candidates[1]

set.seed(seed)
obj_proc <- FindClusters(obj_proc, graph.name = snn_name_use, resolution = best_res, verbose = FALSE)

# Stabilize cluster level ordering
obj_proc$seurat_clusters <- factor(obj_proc$seurat_clusters, levels = sort(unique(obj_proc$seurat_clusters)))

# Save cluster sizes (safe: counts only)
readr::write_csv(
  as.data.frame(table(obj_proc$seurat_clusters)) |> setNames(c("cluster","n_cells")),
  file.path(out_dir, sprintf("cluster_sizes_targetK%d.csv", target_k))
)

#-----------------------------#
# 6) UMAP plots               #
#-----------------------------#
p1 <- DimPlot(obj_proc, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  theme_bw() +
  ggtitle(sprintf("UMAP by cluster (resolution=%.2f; target≈%d)", best_res, target_k))

p2 <- DimPlot(obj_proc, reduction = "umap", group.by = "Group") +
  theme_bw() +
  ggtitle("UMAP by group")

ggsave_safe(file.path(out_dir, sprintf("UMAP_byCluster_targetK%d.tiff", target_k)), p1, width = 7, height = 6, dpi = 300)
ggsave_safe(file.path(out_dir, sprintf("UMAP_byCluster_targetK%d.pdf",  target_k)), p1, width = 7, height = 6)
ggsave_safe(file.path(out_dir, sprintf("UMAP_byGroup_targetK%d.tiff",   target_k)), p2, width = 7, height = 6, dpi = 300)
ggsave_safe(file.path(out_dir, sprintf("UMAP_byGroup_targetK%d.pdf",    target_k)), p2, width = 7, height = 6)

#-----------------------------#
# 7) Marker sets (HM/DAM/IRM) #
#-----------------------------#
HM_genes_raw <- c(
  "Selplg","Cd164","Malat1","Ccr5","P2ry13","P2ry12","Maf","Lpcat2","Srgap2","Cmtm6","Tmem173","Ifngr1","Tmem119","Ecscr",
  "Zfhx3","Marcks","Calm2","Txnip","Gpr56","Alox5ap","Slco2b1","Fam105a","Crybb1","Elmo1","Il6ra","Rnase4","Sparc","Susd3",
  "Cx3cr1","Stab1","Slc2a5","Lpar6","Rgs2","Ptgs1","Fscn1","Sall1","Ssh2","Csf1r","Lrba","Pag1","Pmepa1","Clec4a2","Ltc4s",
  "Rps6ka1","Nrip1","Frmd4a","Bin1","Adap2","Siglech","Tgfbr1","Gpr34","Fgd2","Il10ra","Rab3il1","Ckb","Mef2a","Golm1",
  "Entpd1","Lgmn","Upk1b","Rhob","Olfml3","Nr3c1","Pnp"
)

DAM_genes_raw <- c(
  "Lpl","Apoe","Ctsb","Cst7","Clec7a","Cd9","H2-D1","Ctsd","B2m","Axl","Cd63","Ctsz","Fth1","Lyz2","Serpine2","Ank","Tyrobp",
  "Ctsl","Lgals3bp","Cd52","Trem2","Ccl3","Ccl6","Npc2","H2-K1","Cd68","Csf1","Itgax","Cadm1","Cd63-ps","Cd83","Olfr110","Hexa",
  "Ctsa","Gusb","Rps14","Rplp1","Capg","Baiap2l2","Igf1","Pld3","Rpl23","Mif","Rpl32","Rps28","Rpl18a","Ramp1","Ccl4","Nceh1",
  "Apbb2","Spp1","Gnas"
)

IRM_genes_raw <- c(
  "Ifit1","Ifi213","Ifit3b","Ifi206","Mx1","Oas1g","Isg15","Oas2","Ifi211","Zbp1","Tgtp2","Slfn5","Ifi208","Gm7676","RP24-84E18.1",
  "Oasl1","Rsad2","Iigp1","Cmpk2","Phf11a","Irf7","Ifi207","Ifit3","Rtp4","Ifi214","Gm6904","Isg20","Oas3","Ifit2","Gbp6","Ifitm3",
  "Oas1a","Lgals3bp","Usp18","Ly6e","Stat1","Xaf1","Oasl2","Gm7609","Ms4a4c","Slfn9","Ifi204","Helz2","Ccl2","Ifi47","Slfn2","Sp100",
  "Ifi209","Ifi27l2a","Bst2","Cxcl10","Oas1b","Gm13822","Trim30a","B2m","Ifit1bl1","Ly6a","Rnf213","Adar","Stat2","Gm7592","Phf11b",
  "4930512H18Rik","Gm45418","Ifi44","BC147527","F830016B08Rik","Parp14","Gm8953","H2-T23","Gm44935","Gm4951","Slfn1","Trim30c","Phf11d",
  "Ccl12","Fcgr1","Gm20559","Gbp4","Siglec1","Samd9l","Npsr1","H2-K1","C4b","H2-D1","Irgm1","Slfn8","Gbp2","Lgals9","Ifi35","Gm26797",
  "Trim30d","Parp9","Csprs","Tor3a","Gm15433","Tgtp1","9330175E14Rik","A530064D06Rik","Cxcl9"
)

HM_in  <- map_to_rownames(HM_genes_raw,  rownames(obj_proc))
DAM_in <- map_to_rownames(DAM_genes_raw, rownames(obj_proc))
IRM_in <- map_to_rownames(IRM_genes_raw, rownames(obj_proc))

#-----------------------------#
# 8) Module scores + heatmap  #
#-----------------------------#
if (length(HM_in)  > 0) { obj_proc <- AddModuleScore(obj_proc, list(HM_in),  name = "HM_Score");  obj_proc$HM  <- obj_proc$HM_Score1  }
if (length(DAM_in) > 0) { obj_proc <- AddModuleScore(obj_proc, list(DAM_in), name = "DAM_Score"); obj_proc$DAM <- obj_proc$DAM_Score1 }
if (length(IRM_in) > 0) { obj_proc <- AddModuleScore(obj_proc, list(IRM_in), name = "IRM_Score"); obj_proc$IRM <- obj_proc$IRM_Score1 }

mod_cols <- intersect(colnames(obj_proc@meta.data), c("HM","DAM","IRM"))
avg_scores <- obj_proc@meta.data %>%
  mutate(cluster = obj_proc$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(all_of(mod_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

rownames(avg_scores) <- avg_scores$cluster

if (length(mod_cols) >= 1) {
  # Heatmap expects rows = modules, cols = clusters
  plot_heatmap_module_scores(t(as.matrix(avg_scores[, mod_cols, drop = FALSE])), out_dir)
}

#-----------------------------#
# 9) FindAllMarkers + FGSEA   #
#-----------------------------#
markers <- FindAllMarkers(
  obj_proc, only.pos = TRUE, test.use = "wilcox",
  min.pct = 0.1, logfc.threshold = 0.1
)
write_csv(markers, file.path(out_dir, sprintf("Markers_FindAllMarkers_targetK%d.csv", target_k)))

score_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"

pathways <- list(HM = HM_in, DAM = DAM_in, IRM = IRM_in)
pathways_filtered <- pathways[sapply(pathways, length) >= 10]

clu_levels <- levels(obj_proc$seurat_clusters)
fgsea_list <- list()

for (cl in clu_levels) {
  v <- markers %>%
    filter(cluster == cl) %>%
    select(gene, !!score_col) %>%
    rename(score = !!score_col) %>%
    filter(is.finite(score)) %>%
    group_by(gene) %>%
    summarise(score = max(score), .groups = "drop")

  ranks <- v$score
  names(ranks) <- v$gene
  ranks <- sort(ranks, decreasing = TRUE)

  if (length(ranks) < 20 || length(pathways_filtered) == 0) next

  fg <- suppressWarnings(
    fgsea(pathways = pathways_filtered, stats = ranks,
          nperm = 5000, minSize = 10, maxSize = 5000)
  )
  if (nrow(fg) > 0) {
    fg$cluster <- cl
    fgsea_list[[cl]] <- fg
  }
}

fgsea_res <- if (length(fgsea_list)) bind_rows(fgsea_list) else tibble()
fgsea_res_pos <- fgsea_res %>%
  filter(NES > 0) %>%
  arrange(cluster, padj, desc(NES))

write_csv(fgsea_res_pos, file.path(out_dir, sprintf("FGSEA_byCluster_positiveOnly_targetK%d.csv", target_k)))

# Annotation rule:
#   - Prefer smallest padj (significant first), then highest NES; else "Other"
annot_tbl <- {
  out <- tibble(cluster = clu_levels, annotation = "Other", NES = NA_real_, padj = NA_real_)
  if (nrow(fgsea_res_pos) > 0) {
    best <- fgsea_res_pos %>%
      group_by(cluster) %>%
      arrange(if_else(padj < 0.05, 0L, 1L), desc(NES)) %>%
      slice(1) %>%
      ungroup() %>%
      transmute(cluster, annotation = pathway, NES, padj)

    out$annotation[match(best$cluster, out$cluster)] <- best$annotation
    out$NES[match(best$cluster, out$cluster)] <- best$NES
    out$padj[match(best$cluster, out$cluster)] <- best$padj
  }
  out
}

obj_proc$MicrogliaClass <- factor(
  annot_tbl$annotation[match(obj_proc$seurat_clusters, annot_tbl$cluster)],
  levels = c("HM","DAM","IRM","Other")
)

# UMAP by class
p3 <- DimPlot(obj_proc, reduction = "umap", group.by = "MicrogliaClass") +
  theme_bw() +
  ggtitle("UMAP by microglia class (FGSEA; NES>0)")

ggsave_safe(file.path(out_dir, sprintf("UMAP_byMicrogliaClass_targetK%d.tiff", target_k)), p3, width = 7, height = 6, dpi = 300)
ggsave_safe(file.path(out_dir, sprintf("UMAP_byMicrogliaClass_targetK%d.pdf",  target_k)), p3, width = 7, height = 6)

#-----------------------------#
# 10) DotPlots (paged)        #
#-----------------------------#
DefaultAssay(obj_proc) <- if ("RNA" %in% names(obj_proc@assays)) "RNA" else DefaultAssay(obj_proc)

dot_paged(obj_proc, HM_in,  "HM",  out_dir, per_page = 25)
dot_paged(obj_proc, DAM_in, "DAM", out_dir, per_page = 25)
dot_paged(obj_proc, IRM_in, "IRM", out_dir, per_page = 25)

# Optional: violin plots (paged)
plot_violin_paged(obj_proc, HM_in,  sprintf("HM_targetK%d",  target_k), out_dir)
plot_violin_paged(obj_proc, DAM_in, sprintf("DAM_targetK%d", target_k), out_dir)
plot_violin_paged(obj_proc, IRM_in, sprintf("IRM_targetK%d", target_k), out_dir)

#-----------------------------#
# 11) Composition + statistics#
#-----------------------------#
comp_tbl <- obj_proc@meta.data %>%
  as_tibble(rownames = "cell") %>%
  filter(!is.na(Group), !is.na(MicrogliaClass)) %>%
  count(Group, MicrogliaClass, name = "n")

class_levels <- c("HM","DAM","IRM","Other")
cat_cols <- c("HM"="#88CCEE","DAM"="#CC6677","IRM"="#332288","Other"="#999999")

present_classes <- comp_tbl %>%
  group_by(MicrogliaClass) %>% summarise(total = sum(n), .groups = "drop") %>%
  filter(total > 0) %>% pull(MicrogliaClass) %>% as.character()

ordered_levels <- class_levels[class_levels %in% present_classes]
colors_sub <- cat_cols[ordered_levels]

comp_f <- comp_tbl %>%
  filter(MicrogliaClass %in% ordered_levels) %>%
  group_by(Group) %>% mutate(prop = n / sum(n)) %>% ungroup()

comp_f$MicrogliaClass <- factor(comp_f$MicrogliaClass, levels = ordered_levels)
comp_f$Group <- factor(comp_f$Group, levels = levels(obj_proc$Group))

p_stack <- ggplot(comp_f, aes(x = Group, y = prop, fill = MicrogliaClass)) +
  geom_col(color = "grey20", width = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = colors_sub, breaks = ordered_levels, limits = ordered_levels) +
  labs(x = "Group", y = "Proportion", fill = "Class") +
  theme_bw() +
  theme(legend.position = "right")

ggsave_safe(file.path(out_dir, sprintf("StackedBar_Class_byGroup_targetK%d.tiff", target_k)), p_stack, width = 7.5, height = 5.5, dpi = 300)
ggsave_safe(file.path(out_dir, sprintf("StackedBar_Class_byGroup_targetK%d.pdf",  target_k)), p_stack, width = 7.5, height = 5.5)

# Chi-square test (overall)
tab <- comp_f %>%
  select(Group, MicrogliaClass, n) %>%
  tidyr::pivot_wider(names_from = MicrogliaClass, values_from = n, values_fill = 0) %>%
  arrange(Group)

mat_ct <- as.matrix(tab[, -1, drop = FALSE])
rownames(mat_ct) <- tab$Group

chisq_res <- chisq.test(mat_ct)

capture.output({
  cat("=== Chi-square test of independence ===\n")
  print(addmargins(mat_ct))
  cat("\n")
  print(chisq_res)
}, file = file.path(out_dir, sprintf("Composition_chisq_summary_targetK%d.txt", target_k)))

# Pairwise Fisher tests per class (BH correction)
pairs <- combn(rownames(mat_ct), 2, simplify = FALSE)
res_list <- list()

for (pr in pairs) {
  g1 <- pr[1]; g2 <- pr[2]
  for (cl in colnames(mat_ct)) {
    m2 <- matrix(
      c(
        mat_ct[g1, cl], sum(mat_ct[g1, ]) - mat_ct[g1, cl],
        mat_ct[g2, cl], sum(mat_ct[g2, ]) - mat_ct[g2, cl]
      ),
      nrow = 2, byrow = TRUE
    )
    ft <- fisher.test(m2)
    res_list[[paste(g1, g2, cl, sep = "|")]] <- data.frame(
      group1 = g1, group2 = g2, class = cl,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      n_g1 = sum(mat_ct[g1, ]), n_g2 = sum(mat_ct[g2, ]),
      k_g1 = mat_ct[g1, cl],    k_g2 = mat_ct[g2, cl]
    )
  }
}

res_df <- bind_rows(res_list) %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
write_csv(res_df, file.path(out_dir, sprintf("Composition_Fisher_pairwise_targetK%d.csv", target_k)))

#-----------------------------#
# 12) Save outputs            #
#-----------------------------#
saveRDS(obj_proc, file = file.path(out_dir, sprintf("Seurat_processed_targetK%d.rds", target_k)))

# Public-safe summary (no file paths)
summary_lines <- c(
  sprintf("Cells (subset): %d", ncol(obj_proc)),
  sprintf("Genes: %d", nrow(obj_proc)),
  sprintf("Clusters: %d (resolution=%.2f; target≈%d)", nlevels(obj_proc$seurat_clusters), as.numeric(best_res), target_k),
  sprintf("Module genes present: HM=%d, DAM=%d, IRM=%d", length(HM_in), length(DAM_in), length(IRM_in)),
  sprintf("Seed: %d", seed)
)
writeLines(summary_lines, con = file.path(out_dir, sprintf("_summary_targetK%d.txt", target_k)))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Done.")
