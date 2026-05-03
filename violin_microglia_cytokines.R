################################################################################
# Faceted violin plots of gene expression 
#
# - No hard-coded absolute/local paths
# - Input/output are provided via command-line arguments
# - Assumes an input table with columns: Gene, Expression, Group
# - Produces a horizontal violin plot faceted by Category (rows) and Group (cols)
#
# Usage:
#   Rscript violin_facet_public_safe.R \
#     --in   path/to/expression_data_all.csv \
#     --out  results/violin_facet \
#     --file Supplementary_Fig3_microglia_violin.tiff
#
# Input requirements:
#   - Gene: gene symbol (character)
#   - Expression: numeric (or coercible to numeric)
#   - Group: one of WT-sham, WT-axotomy, KO-sham, KO-axotomy
#
# Notes:
# - Expression is transformed as log1p(Expression).
# - Gene categories and ordering are defined in this script (edit as needed).
################################################################################

suppressPackageStartupMessages({
  pkgs <- c("dplyr", "ggplot2", "readr")
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
  lapply(pkgs, library, character.only = TRUE)
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
out_dir <- get_arg("--out", default = file.path("results", "violin_facet"))
out_file_name <- get_arg("--file", default = "Supplementary_Fig3_microglia_violin.tiff")

dpi_out <- as.integer(get_arg("--dpi", default = "600"))
width_in <- as.numeric(get_arg("--width", default = "16"))
height_in <- as.numeric(get_arg("--height", default = "12"))

if (is.null(in_file)) stop("Missing required argument: --in <path_to_expression_csv>")
if (!file.exists(in_file)) stop("Input file not found (path not printed).")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, out_file_name)

message("Generating faceted violin plot (paths are not printed).")
message("Output will be written under the specified output directory.")

#-----------------------------#
# 1) Load and validate input  #
#-----------------------------#
expression_data_all <- readr::read_csv(in_file, show_col_types = FALSE)

req <- c("Gene", "Expression", "Group")
if (!all(req %in% names(expression_data_all))) {
  stop(
    paste0(
      "Input must contain columns: ", paste(req, collapse = ", "), "\n",
      "Found columns include: ", paste(head(names(expression_data_all), 25), collapse = ", "),
      ifelse(ncol(expression_data_all) > 25, ", ...", "")
    )
  )
}

#-----------------------------#
# 2) Preprocessing            #
#-----------------------------#
expression_data_all <- expression_data_all %>%
  mutate(
    Gene          = as.character(Gene),
    Expression    = suppressWarnings(as.numeric(Expression)),
    LogExpression = log1p(Expression),
    Group = factor(as.character(Group), levels = c("WT-sham", "WT-axotomy", "KO-sham", "KO-axotomy"))
  )

#-----------------------------#
# 3) Gene categories & order  #
#-----------------------------#
# Aif1l is intentionally placed last, with a blank facet label (" ").
gene_categories <- data.frame(
  Gene = c(
    # Pro-inflammatory
    "Il1b", "Il6", "Tnf", "Ifng", "Il12a", "Il12b",
    # Anti-inflammatory
    "Il10", "Tgfb1", "Il4",
    # Chemokines
    "Ccl2", "Cxcl2", "Cxcl10",
    # Growth Factors
    "Csf1", "Csf2", "Il33",
    # Last (no facet title)
    "Aif1l"
  ),
  Category = c(
    rep("Pro-inflammatory Cytokines", 6),
    rep("Anti-inflammatory Cytokines", 3),
    rep("Chemokines", 3),
    rep("Growth Factors", 3),
    " "  # blank label for the final facet row
  ),
  stringsAsFactors = FALSE
)

cat_levels <- c(
  "Pro-inflammatory Cytokines",
  "Anti-inflammatory Cytokines",
  "Chemokines",
  "Growth Factors",
  " "
)
gene_order <- gene_categories$Gene

# Attach category and enforce ordering
expression_data_plot <- expression_data_all
expression_data_plot$Category <- gene_categories$Category[match(expression_data_plot$Gene, gene_categories$Gene)]

expression_data_plot <- expression_data_plot %>%
  filter(!is.na(Category), is.finite(LogExpression))

expression_data_plot$Gene     <- factor(expression_data_plot$Gene, levels = gene_order)
expression_data_plot$Category <- factor(expression_data_plot$Category, levels = cat_levels)

#-----------------------------#
# 4) Plot                     #
#-----------------------------#
# User-specified palette:
# WT-sham = light green, WT-axotomy = dark green, KO-sham = magenta, KO-axotomy = deep pink
cols <- c(
  "WT-sham"    = "#90EE90",
  "WT-axotomy" = "#006400",
  "KO-sham"    = "#FF00FF",
  "KO-axotomy" = "#C71585"
)

p <- ggplot(expression_data_plot, aes(y = Gene, x = LogExpression, fill = Group)) +
  geom_violin(trim = TRUE, scale = "width") +
  # Optional jitter overlay:
  # geom_jitter(aes(color = Group), height = 0.15, size = 0.2, alpha = 0.25, show.legend = FALSE) +
  facet_grid(
    rows = vars(Category),
    cols = vars(Group),
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_manual(values = cols, drop = FALSE) +
  # If jitter is enabled, also add: scale_color_manual(values = cols, drop = FALSE)
  theme_minimal(base_size = 28, base_family = "Arial") +
  theme(
    strip.text.y = element_text(size = 20, face = "bold"),
    strip.text.x = element_text(size = 20, face = "bold"),
    axis.text.y  = element_text(size = 16),
    axis.text.x  = element_text(size = 16),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Gene expression across experimental groups (log-transformed)",
    x = "log(Expression + 1)",
    y = "Genes"
  )

print(p)

#-----------------------------#
# 5) Save (TIFF, LZW)         #
#-----------------------------#
ggsave(
  filename = out_path,
  plot = p,
  device = "tiff",
  dpi = dpi_out,
  width = width_in,
  height = height_in,
  units = "in",
  compression = "lzw"
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done.")
################################################################################
