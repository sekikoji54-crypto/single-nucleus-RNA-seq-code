################################################################################
# Slingshot trajectory inference (UMAP) + pseudotime analysis
#
# - No hard-coded absolute paths
# - No printing of local file paths
# - Designed to be uploaded to a PUBLIC GitHub repo
#
# Usage:
#   Rscript trajectory_slingshot_public_safe.R --rds path/to/input.rds --out results/trajectory
#
# Inputs:
#   --rds : Seurat object (.rds) containing metadata for group and cluster, and UMAP (optional)
#   --out : output directory (default: results/trajectory)
#
# Outputs (in --out):
#   - UMAP lineage plots (TIFF, 600 dpi)
#   - Pseudotime tables (CSV)
#   - Nonparametric tests (CSV)
#   - Barplots for winner-lineage composition and lineage coverage (TIFF)
#   - sessionInfo.txt (no local paths printed)
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

  pkgs <- c(
    "Seurat", "SingleCellExperiment", "slingshot",
    "dplyr", "tibble", "tidyr", "stringr", "readr",
    "ggplot2", "viridis", "grid"
  )
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) BiocManager::install(missing, ask = FALSE, update = FALSE)

  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(viridis)
  library(grid)
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
out_dir  <- get_arg("--out", default = file.path("results", "trajectory"))

if (is.null(rds_path)) stop("Missing required argument: --rds <path_to_seurat_rds>")
if (!file.exists(rds_path)) stop("Input RDS file not found (path not printed for safety).")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Running trajectory pipeline (paths are not printed for public safety).")
message("Outputs will be written under the specified output directory.")

#-----------------------------#
# 1) Helpers (public-safe)    #
#-----------------------------#
to_char <- function(v) {
  if (is.factor(v)) as.character(v) else as.character(v)
}

norm_group <- function(x) {
  # Normalize group labels to lower-case, use '-' separators, and remove trailing '-chat'
  x <- tolower(gsub("_", "-", to_char(x)))
  sub("-chat$", "", x)
}

is_k4_like <- function(u) {
  u <- sort(unique(to_char(u)))
  length(u) == 4 && (all(grepl("^[0-3]$", u)) || all(grepl("^C[0-3]$", toupper(u))))
}

cluster_to_C <- function(v) {
  v <- toupper(trimws(to_char(v)))
  v <- ifelse(grepl("^[0-3]$", v), paste0("C", v), v)
  out <- ifelse(grepl("^C[0-3]$", v), v, NA_character_)
  out
}

find_group_col <- function(meta, target_levels) {
  for (nm in names(meta)) {
    u <- sort(unique(norm_group(meta[[nm]])))
    if (setequal(u, sort(target_levels))) return(nm)
  }
  stop("No metadata column matched the expected group labels (see target_group_levels).")
}

find_cluster_col <- function(meta, exclude = NULL) {
  cand <- setdiff(names(meta), exclude)
  pri  <- cand[grepl("cluster", cand, ignore.case = TRUE)]
  for (nm in c(pri, setdiff(cand, pri))) {
    if (is_k4_like(meta[[nm]])) return(nm)
  }
  stop("No metadata column matched the expected 4-cluster pattern (C0–C3 or 0–3).")
}

save_tiff <- function(filename, width_in = 7, height_in = 6, res_dpi = 600) {
  tiff(
    filename,
    width = width_in, height = height_in, units = "in",
    res = res_dpi, compression = "lzw"
  )
}

#-----------------------------#
# 2) Load input Seurat object #
#-----------------------------#
obj <- readRDS(rds_path)
if (!inherits(obj, "Seurat")) stop("Input object is not a Seurat object.")

# Expected group labels for this project (edit if your project uses other labels)
target_group_levels <- c("wt-sham", "wt-axotomy", "ko-sham", "ko-axotomy")

#-----------------------------#
# 3) Restore group/cluster    #
#-----------------------------#
grp_col <- find_group_col(obj@meta.data, target_group_levels)
cl_col  <- find_cluster_col(obj@meta.data, exclude = grp_col)

obj$.__group__ <- factor(norm_group(obj@meta.data[[grp_col]]), levels = target_group_levels)
obj$.__cluster__ <- cluster_to_C(obj@meta.data[[cl_col]])

stopifnot(!anyNA(obj$.__group__), !anyNA(obj$.__cluster__))

#-----------------------------#
# 4) Ensure UMAP exists       #
#-----------------------------#
DefaultAssay(obj) <- "RNA"
if (!"umap" %in% names(obj@reductions)) {
  obj <- NormalizeData(obj, verbose = FALSE) |>
    FindVariableFeatures(verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = 30, verbose = FALSE) |>
    RunUMAP(dims = 1:30, verbose = FALSE)
}

#-----------------------------#
# 5) Seurat -> SCE -> Slingshot
#-----------------------------#
sce <- as.SingleCellExperiment(obj)
if (!"UMAP" %in% names(reducedDims(sce))) {
  reducedDims(sce)$UMAP <- Embeddings(obj, "umap")
}
colData(sce)$clusterC <- obj$.__cluster__
colData(sce)$group    <- obj$.__group__

# Start cluster fixed to C0, end cluster not specified
sce_sl <- slingshot(
  sce,
  clusterLabels = "clusterC",
  reducedDim    = "UMAP",
  start.clus    = "C0",
  reweight      = TRUE,
  reassign      = TRUE
)

#-----------------------------#
# 6) Plotting + summarization #
#-----------------------------#
plot_and_summarize <- function(sceX, out_dir) {

  pt  <- slingPseudotime(sceX)
  cur <- slingCurves(sceX)

  um <- as.data.frame(reducedDim(sceX, "UMAP"))
  colnames(um) <- c("UMAP1", "UMAP2")

  df <- um |>
    mutate(
      cell    = rownames(um),
      cluster = as.character(colData(sceX)$clusterC),
      group   = as.character(colData(sceX)$group)
    )

  centers <- df |>
    group_by(cluster) |>
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop") |>
    arrange(cluster)

  curve_df <- bind_rows(lapply(seq_along(cur), function(i) {
    data.frame(
      lineage = paste0("L", i),
      x = cur[[i]]$s[, 1],
      y = cur[[i]]$s[, 2],
      order = seq_len(nrow(cur[[i]]$s))
    )
  }))

  arrow_df <- bind_rows(lapply(seq_along(cur), function(i) {
    cr <- cur[[i]]$s
    n  <- nrow(cr)
    data.frame(
      lineage = paste0("L", i),
      x0 = cr[n - 1, 1], y0 = cr[n - 1, 2],
      x1 = cr[n, 1],     y1 = cr[n, 2]
    )
  }))

  # Okabe-Ito palette (color-blind friendly), recycled if needed
  cols_line <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00")
  cols_line <- rep(cols_line, length.out = max(1, length(cur)))

  # Map lineage index -> route (e.g., C0→C2→C3)
  lin <- slingLineages(sceX)
  lin_names <- sapply(lin, function(v) paste(v, collapse = "→"))
  key_df <- data.frame(
    lineage = paste0("Lineage", seq_along(lin_names)),
    path    = unname(lin_names)
  )
  write_csv(key_df, file.path(out_dir, "Lineage_route_key.csv"))

  # ---- Fig 1: lineages over UMAP (colored by cluster) ----
  save_tiff(file.path(out_dir, "UMAP_lineages_byCluster_with_arrows.tiff"))
  print(
    ggplot(df, aes(UMAP1, UMAP2)) +
      geom_point(aes(color = cluster), size = 0.7, alpha = 0.6) +
      lapply(split(curve_df, curve_df$lineage), function(dd) {
        idx <- as.integer(sub("L", "", unique(dd$lineage)))
        geom_path(data = dd, aes(x, y), linewidth = 1.2, color = cols_line[idx])
      }) +
      lapply(split(arrow_df, arrow_df$lineage), function(dd) {
        idx <- as.integer(sub("L", "", unique(dd$lineage)))
        geom_segment(
          data = dd, aes(x = x0, y = y0, xend = x1, yend = y1),
          linewidth = 1.2,
          arrow = arrow(type = "closed", length = unit(0.18, "cm")),
          color = cols_line[idx]
        )
      }) +
      geom_label(
        data = centers, aes(UMAP1, UMAP2, label = cluster),
        size = 4, label.size = 0.2, fill = "white"
      ) +
      theme_bw(base_family = "Arial") +
      labs(title = "Slingshot lineages (colored by cluster)", color = "Cluster")
  )
  dev.off()

  # ---- Fig 2: lineages over UMAP (colored by group) ----
  save_tiff(file.path(out_dir, "UMAP_lineages_byGroup_with_arrows.tiff"))
  print(
    ggplot(df, aes(UMAP1, UMAP2)) +
      geom_point(aes(color = group), size = 0.7, alpha = 0.7) +
      lapply(split(curve_df, curve_df$lineage), function(dd) {
        idx <- as.integer(sub("L", "", unique(dd$lineage)))
        geom_path(data = dd, aes(x, y), linewidth = 1.2, color = cols_line[idx])
      }) +
      lapply(split(arrow_df, arrow_df$lineage), function(dd) {
        idx <- as.integer(sub("L", "", unique(dd$lineage)))
        geom_segment(
          data = dd, aes(x = x0, y = y0, xend = x1, yend = y1),
          linewidth = 1.2,
          arrow = arrow(type = "closed", length = unit(0.18, "cm")),
          color = cols_line[idx]
        )
      }) +
      geom_label(
        data = centers, aes(UMAP1, UMAP2, label = cluster),
        size = 4, label.size = 0.2, fill = "white"
      ) +
      theme_bw(base_family = "Arial") +
      labs(title = "Slingshot lineages (colored by group)", color = "Group")
  )
  dev.off()

  # ---- Fig 3+: pseudotime focus per lineage ----
  for (i in seq_len(ncol(pt))) {
    li <- paste0("L", i)
    df_i <- df
    df_i$pt_i <- pt[, i]
    df_i$in_i <- is.finite(df_i$pt_i)

    route_label <- if (i <= length(lin_names)) lin_names[[i]] else li

    save_tiff(file.path(out_dir, sprintf("UMAP_pseudotime_L%d_focus.tiff", i)))
    print(
      ggplot() +
        geom_point(
          data = subset(df_i, !in_i),
          aes(UMAP1, UMAP2),
          color = "grey85", size = 0.6, alpha = 0.6
        ) +
        geom_point(
          data = subset(df_i, in_i),
          aes(UMAP1, UMAP2, color = pt_i),
          size = 0.8, alpha = 0.9
        ) +
        { dd <- subset(curve_df, lineage == li)
          geom_path(data = dd, aes(x, y), linewidth = 1.3, color = cols_line[i]) } +
        { ad <- subset(arrow_df, lineage == li)
          geom_segment(
            data = ad, aes(x = x0, y = y0, xend = x1, yend = y1),
            linewidth = 1.3,
            arrow = arrow(type = "closed", length = unit(0.18, "cm")),
            color = cols_line[i]
          ) } +
        geom_label(
          data = centers, aes(UMAP1, UMAP2, label = cluster),
          size = 4, label.size = 0.2, fill = "white"
        ) +
        scale_color_viridis_c(option = "plasma", na.value = "grey80") +
        theme_bw(base_family = "Arial") +
        labs(
          title = sprintf("Pseudotime (focused): %s (Lineage %d)", route_label, i),
          color = "Pseudotime"
        )
    )
    dev.off()
  }

  # ---- Tables (long format pseudotime) ----
  pt_long <- as.data.frame(pt) |>
    rownames_to_column("cell") |>
    pivot_longer(-cell, names_to = "lineage", values_to = "pseudotime") |>
    left_join(df[, c("cell", "group", "cluster")], by = "cell") |>
    filter(is.finite(pseudotime))

  write_csv(pt_long, file.path(out_dir, "Pseudotime_allLineages_long.csv"))

  # Summary by group
  by_group <- pt_long |>
    group_by(lineage, group) |>
    summarise(
      n      = n(),
      mean   = mean(pseudotime),
      median = median(pseudotime),
      sd     = sd(pseudotime),
      .groups = "drop"
    )
  write_csv(by_group, file.path(out_dir, "Pseudotime_byGroup_summary.csv"))

  # Kruskal-Wallis per lineage
  kw <- pt_long |>
    group_by(lineage) |>
    summarise(KW_p = kruskal.test(pseudotime ~ group)$p.value, .groups = "drop")
  write_csv(kw, file.path(out_dir, "KW_test_pvalues.csv"))

  # Pairwise Wilcoxon (BH-adjusted)
  for (lin in unique(pt_long$lineage)) {
    dat <- subset(pt_long, lineage == lin)
    pw  <- pairwise.wilcox.test(dat$pseudotime, dat$group, p.adjust.method = "BH")
    m <- as.data.frame(pw$p.value) |>
      rownames_to_column("grp1") |>
      pivot_longer(-grp1, names_to = "grp2", values_to = "p_adj") |>
      filter(!is.na(p_adj)) |>
      mutate(lineage = lin)
    write_csv(m, file.path(out_dir, paste0("PairwiseWilcox_", lin, ".csv")))
  }

  # Violin plots per lineage
  for (lin in unique(pt_long$lineage)) {
    dat <- subset(pt_long, lineage == lin)
    save_tiff(file.path(out_dir, paste0("Violin_pseudotime_byGroup_", lin, ".tiff")))
    print(
      ggplot(dat, aes(x = group, y = pseudotime, fill = group)) +
        geom_violin(trim = FALSE, alpha = 0.7) +
        geom_boxplot(width = 0.15, outlier.size = 0.7, alpha = 0.9) +
        theme_bw(base_family = "Arial") +
        labs(title = paste0("Pseudotime by group (", lin, ")"), x = NULL, y = "Pseudotime") +
        theme(legend.position = "none")
    )
    dev.off()
  }

  # Winner-lineage assignment:
  # - Normalize each lineage pseudotime to 0–1 within each lineage column
  # - Assign each cell to the lineage with the maximum normalized pseudotime
  wide <- pt_long |>
    select(cell, group, lineage, pseudotime) |>
    pivot_wider(names_from = lineage, values_from = pseudotime)

  lin_cols <- setdiff(names(wide), c("cell", "group"))

  scale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (is.finite(rng[1]) && diff(rng) > 0) (x - rng[1]) / diff(rng) else rep(0, length(x))
  }

  wide_norm <- wide
  wide_norm[lin_cols] <- lapply(wide[lin_cols], scale01)

  mat  <- as.matrix(wide_norm[lin_cols])
  mat2 <- mat
  mat2[!is.finite(mat2)] <- -Inf

  winner_idx <- max.col(mat2, ties.method = "first")
  row_all_na <- apply(is.finite(as.matrix(wide[lin_cols])), 1, function(v) !any(v))

  winner_lineage <- rep(NA_character_, nrow(mat))
  winner_lineage[!row_all_na] <- lin_cols[winner_idx[!row_all_na]]

  assign_df <- wide |>
    select(cell, group) |>
    mutate(winner_lineage = winner_lineage)

  order_groups <- c("wt-sham", "wt-axotomy", "ko-sham", "ko-axotomy")
  assign_df$group <- factor(as.character(assign_df$group), levels = order_groups)

  counts <- assign_df |>
    filter(!is.na(winner_lineage)) |>
    group_by(group, winner_lineage) |>
    summarise(n_cells = n(), .groups = "drop") |>
    group_by(group) |>
    mutate(total = sum(n_cells), percent = 100 * n_cells / total) |>
    ungroup()

  write_csv(counts, file.path(out_dir, "WinnerLineage_counts_perc.csv"))

  # Coverage (not mutually exclusive): % cells with finite pseudotime per lineage
  cover <- pt_long |>
    mutate(has = is.finite(pseudotime)) |>
    group_by(group, lineage) |>
    summarise(coverage_percent = 100 * mean(has), .groups = "drop")

  write_csv(cover, file.path(out_dir, "Coverage_percent_byGroup.csv"))

  # Barplots
  plt1 <- ggplot(counts, aes(x = group, y = percent, fill = winner_lineage)) +
    geom_col(width = 0.75, color = "white") +
    theme_bw(base_family = "Arial") +
    labs(title = "Winner-lineage composition by group", x = NULL, y = "% of cells", fill = "Winner")

  save_tiff(file.path(out_dir, "Bar_WinnerLineage_byGroup.tiff"))
  print(plt1)
  dev.off()

  plt2 <- ggplot(cover, aes(x = group, y = coverage_percent, fill = lineage)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    theme_bw(base_family = "Arial") +
    labs(title = "Lineage coverage by group", x = NULL, y = "% cells with finite pseudotime", fill = "Lineage")

  save_tiff(file.path(out_dir, "Bar_Coverage_byGroup.tiff"))
  print(plt2)
  dev.off()
}

plot_and_summarize(sce_sl, out_dir)

# Save session info (no local file paths printed here)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Done.")
