################################################################################
# NicheNet pipeline (Receiver: ChAT; Sender: Microglia) 
#
# Comparisons (default):
#   1) KO-axotomy vs KO-sham
#   2) WT-axotomy vs WT-sham
#
# Main outputs (per comparison):
#   - DEG tables (receiver)
#   - Ligand activity (AUPR corrected)
#   - Ligand→target links (top)
#   - Ligand–receptor candidates (filtered by expression)
#   - Barplot / Heatmap / Chord / Sankey / DotPlot (receptors)
#
# Cross-comparison outputs:
#   - WT vs KO ligand activity scatter
#   - Top15 overlap CSV + barplot + UpSet
#
# Usage:
#   Rscript nichenet_chat_from_mgl_public_safe.R \
#     --chat_rds path/to/seurat_ChAT_4groups.rds \
#     --mgl_rds  path/to/seurat_Microglia_4groups.rds \
#     --group_col group \
#     --out results/NicheNet_ChAT_from_MGL
#
# Notes (public safety):
#   - No absolute paths are hard-coded or printed.
#   - Output filenames do not include local paths.
################################################################################

set.seed(1)

#-----------------------------#
# 0) Dependencies             #
#-----------------------------#
need <- c(
  "Seurat","dplyr","tibble","readr","stringr","tidyr","purrr","Matrix",
  "nichenetr","AnnotationDbi","org.Mm.eg.db","babelgene",
  "ggplot2","ggrepel","ggalluvial","ComplexUpset","pheatmap",
  "circlize","networkD3","htmlwidgets"
)

miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  message("Missing packages detected (names only): ", paste(miss, collapse = ", "))
  stop("Please install the missing packages, then rerun this script.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(Matrix)

  library(nichenetr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(babelgene)

  library(ggplot2)
  library(ggrepel)
  library(ggalluvial)
  library(ComplexUpset)
  library(pheatmap)
  library(circlize)
  library(networkD3)
  library(htmlwidgets)
})

#-----------------------------#
# 1) Argument parsing         #
#-----------------------------#
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 1 && length(args) >= hit + 1) return(args[hit + 1])
  default
}

path_chat <- get_arg("--chat_rds")
path_mgl  <- get_arg("--mgl_rds")
group_col <- get_arg("--group_col", default = "group")
out_dir   <- get_arg("--out", default = file.path("results", "NicheNet_ChAT_from_MGL"))

if (is.null(path_chat)) stop("Missing required argument: --chat_rds <path_to_receiver_seurat_rds>")
if (is.null(path_mgl))  stop("Missing required argument: --mgl_rds <path_to_sender_seurat_rds>")
if (!file.exists(path_chat)) stop("Receiver RDS not found (path not printed).")
if (!file.exists(path_mgl))  stop("Sender RDS not found (path not printed).")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Starting NicheNet pipeline (paths are not printed).")
message("Outputs will be written under the specified output directory.")

#-----------------------------#
# 2) User-configurable values #
#-----------------------------#
cmp_list <- list(
  KO = list(case = "KO-axotomy", ctrl = "KO-sham", tag = "KO_axotomy_vs_KO_sham"),
  WT = list(case = "WT-axotomy", ctrl = "WT-sham", tag = "WT_axotomy_vs_WT_sham")
)

# Optional: UMAP FeaturePlot targets (only for WT comparison in this script)
genes_umap_chat <- c(
  "Aplp2","Ddr1","Dysf","Egfr","F2r","Fgfr1","Igsf11",
  "Itga3","Itga6","Itga7","Itga9","Itgav","Itgb1","Itgb4","Itgb5",
  "Plat","Sdc2"
)
genes_umap_mgl  <- c(
  "Anxa2","F13a1","H2-D1","H2-K1","H2-Q4","H2-Q6","H2-Q7",
  "Mmp2","Tgfbi","Vsir"
)

#-----------------------------#
# 3) Utilities                #
#-----------------------------#
normalize_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[nchar(x) == 0] <- NA_character_
  x
}

extract_symbolish <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.[0-9]+$", "", x)  # strip trailing ".1" etc
  x <- gsub("(?i)^(gene|symbol|target|geneid|ensembl)[:_ -]*", "", x, perl = TRUE)
  x <- gsub("(?i)[:_ -]*(gene|symbol|target)$", "", x, perl = TRUE)
  x <- gsub("\\s+", "", x)
  x
}

expressed_genes <- function(seu, pct = 0.05, assay = NULL) {
  if (!is.null(assay)) DefaultAssay(seu) <- assay
  mat <- GetAssayData(seu, slot = "data")
  frac <- Matrix::rowMeans(mat > 0)
  names(frac)[frac >= pct]
}

# Convert NicheNet extdata (human symbols) to mouse symbols via babelgene orthologs
# Returns: list(ltm = ligand_target_matrix_mouse, lr = lr_network_mouse)
relabel_to_mouse <- function(ltm, lr) {
  mapL <- babelgene::orthologs(colnames(ltm), species = "mouse")
  mapR <- babelgene::orthologs(rownames(ltm), species = "mouse")

  mapL <- mapL[!is.na(mapL$mouse_symbol) & nzchar(mapL$mouse_symbol),
               c("human_symbol","mouse_symbol")] %>% distinct()
  mapR <- mapR[!is.na(mapR$mouse_symbol) & nzchar(mapR$mouse_symbol),
               c("human_symbol","mouse_symbol")] %>% distinct()

  col_new <- setNames(mapL$mouse_symbol, mapL$human_symbol)[colnames(ltm)]
  row_new <- setNames(mapR$mouse_symbol, mapR$human_symbol)[rownames(ltm)]

  ltm_m <- as.matrix(ltm)

  # Collapse duplicated ligands (columns)
  col_dup <- split(seq_along(col_new), col_new)
  ltm2 <- do.call(cbind, lapply(col_dup, function(idx) Matrix::rowSums(ltm_m[, idx, drop = FALSE])))
  # Apply receptor renaming
  rownames(ltm2) <- row_new

  # Collapse duplicated receptors (rows)
  row_dup <- split(seq_len(nrow(ltm2)), rownames(ltm2))
  ltm3 <- do.call(rbind, lapply(row_dup, function(idx) Matrix::colSums(ltm2[idx, , drop = FALSE])))

  # LR network renaming
  lr_m <- lr %>%
    mutate(ligand_h = ligand, receptor_h = receiver) %>%
    left_join(rename(mapL, ligand_h = human_symbol, ligand = mouse_symbol), by = "ligand_h") %>%
    left_join(rename(mapL, receptor_h = human_symbol, receiver = mouse_symbol), by = "receptor_h") %>%
    transmute(
      ligand  = coalesce(ligand, ligand_h),
      receiver = coalesce(receiver, receptor_h)
    ) %>%
    distinct()

  list(ltm = ltm3, lr = lr_m)
}

get_links_top <- function(ligands, ligand_target_matrix, n = 100) {
  if (missing(ligands) || !length(ligands)) {
    return(tibble(ligand = character(), target = character(), weight = numeric()))
  }
  lig_all <- colnames(ligand_target_matrix)
  map_col <- setNames(lig_all, normalize_symbol(lig_all))
  lig_use <- intersect(normalize_symbol(ligands), names(map_col))
  if (!length(lig_use)) {
    return(tibble(ligand = character(), target = character(), weight = numeric()))
  }

  cols <- unname(map_col[lig_use])
  bind_rows(lapply(cols, function(cx) {
    sc <- ligand_target_matrix[, cx, drop = FALSE]
    ord <- order(sc[, 1], decreasing = TRUE)
    take <- head(ord, n = min(n, length(ord)))
    tibble(
      ligand = colnames(sc)[1],
      target = rownames(sc)[take],
      weight = as.numeric(sc[take, 1])
    )
  }))
}

build_node_colors <- function(ligands, receptors) {
  col_lig <- "#4477AA"
  col_rec <- "#AA7733"
  nodes <- unique(c(ligands, receptors))
  side  <- ifelse(nodes %in% ligands, "ligand", "receptor")
  setNames(ifelse(side == "ligand", col_lig, col_rec), nodes)
}

save_dual <- function(p, base, w = 8, h = 5, dpi = 400) {
  ggsave(paste0(base, ".pdf"),  p, width = w, height = h, units = "in", device = "pdf", useDingbats = FALSE)
  ggsave(paste0(base, ".tiff"), p, width = w, height = h, units = "in", dpi = dpi, device = "tiff", compression = "lzw")
}

save_heatmap_dual <- function(mat, base, ann_col = NULL, w = 8, h = 10, dpi = 400) {
  grDevices::pdf(paste0(base, ".pdf"), width = w, height = h)
  pheatmap::pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                     annotation_col = ann_col, show_rownames = FALSE)
  grDevices::dev.off()

  grDevices::tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  pheatmap::pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                     annotation_col = ann_col, show_rownames = FALSE)
  grDevices::dev.off()
}

make_chord <- function(df_pairs, outbase, top_links = 60, dpi = 400, size_in = 10) {
  if (nrow(df_pairs) == 0) return(invisible(NULL))
  wcol <- if ("sumscore" %in% names(df_pairs)) "sumscore" else if ("count" %in% names(df_pairs)) "count" else names(df_pairs)[3]
  df_top <- df_pairs %>% arrange(desc(.data[[wcol]])) %>% slice_head(n = min(top_links, n()))

  ligands   <- unique(df_top$ligand)
  receptors <- unique(df_top$receptor)
  grid.col  <- build_node_colors(ligands, receptors)

  draw_chord <- function(device_fun) {
    circos.clear()
    circos.par(
      gap.after = c(rep(2, max(0, length(ligands) - 1)), 8, rep(2, max(0, length(receptors) - 1)), 8),
      start.degree = 90,
      track.margin = c(0.01, 0.01)
    )

    device_fun()
    chordDiagram(
      df_top %>% select(ligand, receptor, !!wcol),
      grid.col = grid.col,
      order = c(ligands, receptors),
      directional = 0,
      transparency = 0.25,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = 0.08)
    )
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index,
                  facing = "bending", niceFacing = TRUE, cex = 0.6)
    }, bg.border = NA)
    grDevices::dev.off()
    circos.clear()
  }

  draw_chord(function() grDevices::pdf(paste0(outbase, ".pdf"), width = size_in, height = size_in))
  draw_chord(function() grDevices::tiff(paste0(outbase, ".tiff"), width = size_in, height = size_in, units = "in",
                                        res = dpi, compression = "lzw"))
}

dotplot_receptors <- function(seu, receptors, base, group_col, order_groups = NULL, ordered = FALSE, dpi = 400) {
  receptors <- intersect(receptors, rownames(seu))
  if (!length(receptors)) return(invisible(NULL))

  if (!group_col %in% colnames(seu@meta.data)) return(invisible(NULL))

  if (!is.null(order_groups)) {
    seu[[group_col]] <- factor(seu[[group_col]][, 1], levels = order_groups)
  }
  Idents(seu) <- seu[[group_col]][, 1]

  feats <- receptors
  if (ordered) {
    avg <- AverageExpression(seu, features = feats, assays = DefaultAssay(seu), group.by = group_col)$RNA
    ord <- order(rowMeans(avg), decreasing = TRUE)
    feats <- rownames(avg)[ord]
  }

  p <- DotPlot(seu, features = feats, cols = c("#BBBBBB", "#4477AA")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  save_dual(p, base, w = 10, h = 6, dpi = dpi)
}

feature_umap <- function(seu, genes, out_dir, prefix, dpi = 400, tiff_only_for = NULL) {
  if (!"umap" %in% Reductions(seu)) {
    seu <- RunPCA(seu, verbose = FALSE)
    seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  }
  for (g in genes) {
    if (!g %in% rownames(seu)) next
    p <- FeaturePlot(seu, features = g, reduction = "umap") + theme_void() + ggtitle(g)
    pdf_path <- file.path(out_dir, sprintf("%s_%s.pdf", prefix, g))
    ggsave(pdf_path, p, width = 5, height = 4, units = "in", device = "pdf", useDingbats = FALSE)

    if (is.null(tiff_only_for) || g %in% tiff_only_for) {
      tiff_path <- file.path(out_dir, sprintf("%s_%s.tiff", prefix, g))
      ggsave(tiff_path, p, width = 5, height = 4, units = "in", dpi = dpi, device = "tiff", compression = "lzw")
    }
  }
}

sankey_alluvial <- function(df_pairs, base, top_links = 60, dpi = 400) {
  if (nrow(df_pairs) == 0) return(invisible(NULL))
  if (!"sumscore" %in% names(df_pairs)) return(invisible(NULL))

  df_top <- df_pairs %>% arrange(desc(sumscore)) %>% slice_head(n = min(top_links, n()))

  # Static (ggalluvial)
  p <- ggplot(df_top, aes(axis1 = ligand, axis2 = receptor, y = sumscore)) +
    scale_x_discrete(limits = c("Ligand", "Receptor"), expand = c(.1, .1)) +
    geom_alluvium(aes(fill = ligand), alpha = .6) +
    geom_stratum(width = 1/10) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")

  save_dual(p, base, w = 12, h = 6, dpi = dpi)

  # Interactive HTML (networkD3)
  try({
    ligs <- unique(df_top$ligand)
    recs <- unique(df_top$receptor)
    nodes <- data.frame(name = c(ligs, recs))

    links <- df_top %>%
      mutate(
        source = match(ligand, nodes$name) - 1,
        target = match(receptor, nodes$name) - 1,
        value  = sumscore
      ) %>%
      select(source, target, value)

    sn <- sankeyNetwork(
      Links = links, Nodes = nodes,
      Source = "source", Target = "target",
      Value = "value", NodeID = "name",
      fontSize = 12, nodeWidth = 20
    )

    htmlwidgets::saveWidget(sn, file = paste0(base, ".html"), selfcontained = TRUE)
  }, silent = TRUE)
}

#-----------------------------#
# 4) Load input data          #
#-----------------------------#
chat0 <- readRDS(path_chat)
mgl0  <- readRDS(path_mgl)

if (!inherits(chat0, "Seurat")) stop("Receiver object is not a Seurat object.")
if (!inherits(mgl0,  "Seurat")) stop("Sender object is not a Seurat object.")
if (!group_col %in% colnames(chat0@meta.data)) stop("group_col not found in receiver meta.data.")
if (!group_col %in% colnames(mgl0@meta.data))  stop("group_col not found in sender meta.data.")

#-----------------------------#
# 5) NicheNet networks (mouse)#
#-----------------------------#
ligand_target_matrix_h <- readRDS(system.file("extdata", "ligand_target_matrix.rds", package = "nichenetr"))
lr_network_h           <- readRDS(system.file("extdata", "lr_network.rds", package = "nichenetr"))

mm <- relabel_to_mouse(ligand_target_matrix_h, lr_network_h)
ligand_target_matrix_m <- mm$ltm
lr_network_m           <- mm$lr
rm(mm, ligand_target_matrix_h, lr_network_h)

#-----------------------------#
# 6) Core runner per compare  #
#-----------------------------#
run_cmp <- function(tag, case_label, ctrl_label,
                    chat0, mgl0, group_col, out_dir,
                    ligand_target_matrix_m, lr_network_m) {

  message("=== Running comparison: ", tag, " ===")

  # Receiver (ChAT): restrict to 2 groups
  chat <- subset(chat0, subset = !!as.name(group_col) %in% c(case_label, ctrl_label))
  Idents(chat) <- chat[[group_col]][, 1]

  # Sender (Microglia): restrict to same 2 groups for expression filters
  mgl <- subset(mgl0, subset = !!as.name(group_col) %in% c(case_label, ctrl_label))

  # DEG in receiver
  deg <- FindMarkers(
    chat, ident.1 = case_label, ident.2 = ctrl_label,
    logfc.threshold = 0.1, min.pct = 0.05, test.use = "wilcox",
    verbose = FALSE
  ) %>%
    rownames_to_column("gene") %>%
    mutate(gene = extract_symbolish(gene))

  write_csv(deg, file.path(out_dir, sprintf("DEG_receiver_%s.csv", tag)))

  geneset_up <- deg %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    pull(gene) %>%
    unique()

  background_expressed <- expressed_genes(chat, pct = 0.05)

  # Expression-constrained LR network
  sender_expressed   <- expressed_genes(mgl,  pct = 0.05)
  receiver_expressed <- expressed_genes(chat, pct = 0.05)

  lr_expressed <- lr_network_m %>%
    filter(ligand %in% sender_expressed, receiver %in% receiver_expressed) %>%
    distinct()

  potential_ligands <- intersect(lr_expressed$ligand, colnames(ligand_target_matrix_m))

  ligand_activities <- get_ligand_activities(
    geneset = geneset_up,
    background_expressed_genes = background_expressed,
    ligand_target_matrix = ligand_target_matrix_m,
    potential_ligands = potential_ligands
  )

  write_csv(ligand_activities, file.path(out_dir, sprintf("LigandActivity_mouse_%s.csv", tag)))

  # Top ligands
  topN <- 30
  top_ligands <- ligand_activities %>%
    arrange(desc(aupr_corrected)) %>%
    slice_head(n = min(topN, n())) %>%
    pull(test_ligand)

  summary10 <- ligand_activities %>%
    arrange(desc(aupr_corrected)) %>%
    mutate(rank = row_number()) %>%
    slice_head(n = min(10, n()))

  write_csv(summary10, file.path(out_dir, sprintf("Summary_top10_ligands_%s.csv", tag)))

  # Ligand→target links (top per ligand)
  links_top <- get_links_top(top_ligands, ligand_target_matrix_m, n = 100)
  write_csv(links_top, file.path(out_dir, sprintf("LigandTargetLinks_mouse_%s.csv", tag)))

  # Aggregate to LR candidates
  ligand_score <- links_top %>%
    group_by(ligand) %>%
    summarise(sumscore = sum(weight), count = dplyr::n(), .groups = "drop")

  lr_pairs <- lr_expressed %>%
    inner_join(ligand_score, by = "ligand") %>%
    transmute(ligand, receptor = receiver, sumscore, count) %>%
    arrange(desc(sumscore))

  write_csv(lr_pairs, file.path(out_dir, sprintf("LR_candidates_mouse_%s.csv", tag)))

  # Barplot: top ligands
  p_bar <- ligand_activities %>%
    arrange(desc(aupr_corrected)) %>%
    slice_head(n = 20) %>%
    mutate(test_ligand = factor(test_ligand, levels = rev(test_ligand))) %>%
    ggplot(aes(x = test_ligand, y = aupr_corrected)) +
    geom_col() +
    coord_flip() +
    labs(x = "Ligand", y = "AUPR (corrected)", title = tag) +
    theme_minimal(base_size = 11)

  save_dual(p_bar, file.path(out_dir, sprintf("Barplot_TopLigands_%s", tag)), w = 6, h = 6)

  # Heatmap: ligand-target (top links)
  top_targets <- links_top %>%
    arrange(desc(weight)) %>%
    slice_head(n = 200) %>%
    pull(target) %>%
    unique()

  mat_ht <- ligand_target_matrix_m[
    intersect(top_targets, rownames(ligand_target_matrix_m)),
    intersect(top_ligands, colnames(ligand_target_matrix_m)),
    drop = FALSE
  ]

  if (nrow(mat_ht) > 1 && ncol(mat_ht) > 1) {
    save_heatmap_dual(mat = mat_ht,
                      base = file.path(out_dir, sprintf("Heatmap_LigandTarget_mouse_%s", tag)),
                      w = 8, h = 10)
  }

  # Chord + Sankey for LR candidates
  make_chord(lr_pairs, outbase = file.path(out_dir, sprintf("Chord_LR_%s", tag)), top_links = 60)
  sankey_alluvial(lr_pairs, base = file.path(out_dir, sprintf("Sankey_LR_%s", tag)), top_links = 60)

  # DotPlot: receptor expression in receiver
  rec_all <- unique(lr_expressed$receiver)
  rec_top <- lr_pairs %>% arrange(desc(sumscore)) %>% pull(receptor) %>% unique() %>% head(50)

  dotplot_receptors(chat, rec_all,
                    base = file.path(out_dir, sprintf("DotPlot_receiver_receptors_%s", tag)),
                    group_col = group_col, order_groups = NULL, ordered = FALSE)

  dotplot_receptors(chat, rec_all,
                    base = file.path(out_dir, sprintf("DotPlot_receiver_receptors_%s_ordered", tag)),
                    group_col = group_col, order_groups = NULL, ordered = TRUE)

  # Optional canonical ordering (useful when plotting all groups together in a shared object)
  canonical_order <- c("WT-sham","WT-axotomy","KO-sham","KO-axotomy")
  dotplot_receptors(chat, rec_all,
                    base = file.path(out_dir, sprintf("DotPlot_receiver_receptors_canonical_order_%s", tag)),
                    group_col = group_col, order_groups = canonical_order, ordered = TRUE)

  dotplot_receptors(chat, rec_top,
                    base = file.path(out_dir, sprintf("DotPlot_receiver_receptors_topLigands_%s", tag)),
                    group_col = group_col, order_groups = NULL, ordered = TRUE)

  list(
    tag = tag,
    activities = ligand_activities %>%
      select(test_ligand, aupr_corrected) %>%
      rename(!!tag := aupr_corrected),
    top15 = ligand_activities %>% arrange(desc(aupr_corrected)) %>% slice_head(n = 15) %>% pull(test_ligand)
  )
}

#-----------------------------#
# 7) Run comparisons          #
#-----------------------------#
res_KO <- run_cmp(
  tag = cmp_list$KO$tag, case_label = cmp_list$KO$case, ctrl_label = cmp_list$KO$ctrl,
  chat0 = chat0, mgl0 = mgl0, group_col = group_col, out_dir = out_dir,
  ligand_target_matrix_m = ligand_target_matrix_m, lr_network_m = lr_network_m
)

res_WT <- run_cmp(
  tag = cmp_list$WT$tag, case_label = cmp_list$WT$case, ctrl_label = cmp_list$WT$ctrl,
  chat0 = chat0, mgl0 = mgl0, group_col = group_col, out_dir = out_dir,
  ligand_target_matrix_m = ligand_target_matrix_m, lr_network_m = lr_network_m
)

#-----------------------------#
# 8) WT vs KO summary         #
#-----------------------------#
act_join <- full_join(res_WT$activities, res_KO$activities, by = "test_ligand") %>%
  replace_na(setNames(as.list(rep(0, 2)), c(cmp_list$WT$tag, cmp_list$KO$tag)))

write_csv(act_join, file.path(out_dir, "LigandActivity_WT_vs_KO.csv"))

p_sc <- ggplot(act_join, aes(x = .data[[cmp_list$WT$tag]], y = .data[[cmp_list$KO$tag]], label = test_ligand)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  labs(
    x = paste0("AUPR: ", cmp_list$WT$tag),
    y = paste0("AUPR: ", cmp_list$KO$tag),
    title = "Ligand activities: WT vs KO"
  ) +
  theme_minimal(base_size = 11)

save_dual(p_sc, file.path(out_dir, "Scatter_LigandActivity_WT_vs_KO"), w = 6, h = 6)

# Top15 overlap
top15_WT <- res_WT$top15
top15_KO <- res_KO$top15
ovl      <- intersect(top15_WT, top15_KO)
only_WT  <- setdiff(top15_WT, top15_KO)
only_KO  <- setdiff(top15_KO, top15_WT)

write_csv(tibble(ligand = ovl),     file.path(out_dir, "TopLigandOverlap_WTvsKO_top15.csv"))
write_csv(tibble(ligand = only_WT), file.path(out_dir, "TopLigands_WTminusKO.csv"))
write_csv(tibble(ligand = only_KO), file.path(out_dir, "TopLigands_KOminusWT.csv"))

# Overlap barplot for shared ligands
if (length(ovl) > 0) {
  act_for_bar <- act_join %>%
    filter(test_ligand %in% ovl) %>%
    pivot_longer(cols = -test_ligand, names_to = "comparison", values_to = "AUPR")

  p_bar_ovl <- ggplot(act_for_bar, aes(x = reorder(test_ligand, AUPR, FUN = median), y = AUPR, fill = comparison)) +
    geom_col(position = "dodge") +
    coord_flip() +
    labs(x = "Ligands (overlap in top15)", y = "AUPR", title = "WT vs KO overlap (top15)") +
    theme_minimal(base_size = 11)

  save_dual(p_bar_ovl, file.path(out_dir, "Bar_OverlapTopLigands_WTvsKO_top15"), w = 7, h = 6)
}

# UpSet (ComplexUpset)
sets_df <- tibble(
  ligand = unique(c(top15_WT, top15_KO)),
  WT = ligand %in% top15_WT,
  KO = ligand %in% top15_KO
)

p_up <- ComplexUpset::upset(
  sets_df,
  sets = c("WT", "KO"),
  name = "Top15 ligands",
  base_annotations = list("Intersection size" = intersection_size())
) +
  ggtitle("Overlap of Top15 ligands (WT vs KO)")

save_dual(p_up, file.path(out_dir, "UpSet_TopLigands_WTvsKO_top15"), w = 8, h = 5)

#-----------------------------#
# 9) Optional UMAP plots (WT) #
#-----------------------------#
# Receiver (ChAT)
chat_WT <- subset(chat0, subset = !!as.name(group_col) %in% c("WT-sham", "WT-axotomy"))
feature_umap(chat_WT, genes_umap_chat, out_dir = out_dir,
             prefix = "UMAP_receiver_WT_axotomy_vs_WT_sham", tiff_only_for = c("Itga9"))

# Sender (Microglia)
mgl_WT <- subset(mgl0, subset = !!as.name(group_col) %in% c("WT-sham", "WT-axotomy"))
feature_umap(mgl_WT, genes_umap_mgl, out_dir = out_dir,
             prefix = "UMAP_sender_WT_axotomy_vs_WT_sham", tiff_only_for = c("F13a1"))

#-----------------------------#
# 10) Reproducibility         #
#-----------------------------#
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("All outputs completed.")
