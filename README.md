# single-nucleus-RNA-seq-code

This repository contains R scripts used for the single-nucleus RNA-seq analyses in:

Iba1 deficiency impairs microglial synaptic remodeling and neuronal survival after axonal injury.

## Contents

- `cluster_vs_rest_deg_volcano_composition.R`: cluster-vs-rest DEG analysis and volcano plots.
- `marker_dotplot_graph_clusters.R`: marker-based dot plot for graph-based clusters.
- `go_enrichment_cluster_degs.R`: GO enrichment analysis for cluster-specific DEGs.
- `microglia_reclustering_ad_state_annotation.R`: microglial reclustering and HM/IRM/DAM-like annotation.
- `tail_score_module.R`: axotomy-related microglial module score analysis.
- `slingshot_pseudotime_microglia.R`: Slingshot trajectory and pseudotime analysis.
- `nichenet_microglia_to_chat.R`: NicheNet analysis of microglia-to-ChAT neuron signaling.
- `violin_microglia_cytokines.R`: violin plots of selected cytokine, chemokine, growth-factor, and Aif1l expression.
- `volcano_highlight_validation_genes.R`: volcano plot highlighting selected validation genes.

## Data availability

Raw and processed datasets are not included in this repository. Additional datasets are available from the corresponding author upon reasonable request, as described in the article.

## Requirements

The scripts were developed in R and use packages including Seurat, ggplot2, dplyr, readr, ggrepel, clusterProfiler, org.Mm.eg.db, fgsea, slingshot, SingleCellExperiment, and NicheNet.

## Usage

Each script is designed to be run from the command line. Example:

```bash
Rscript cluster_vs_rest_deg_volcano_composition.R \
  --rds path/to/seurat_object.rds \
  --out results/cluster_deg_volcano
