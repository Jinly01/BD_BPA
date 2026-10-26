# ssGSEA immune infiltration.
# 28 cell type signatures from mmc3.xlsx, group comparison, cell-cell and
# core gene vs cell correlations, abundance plots.

# 0. Packages
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_packages <- c("rio", "dplyr", "tidyr", "tibble", "ggplot2", "ggpubr",
                   "corrplot", "pheatmap", "RColorBrewer", "stringr",
                   "rstatix", "grid")
bioc_packages <- c("GSVA", "limma")

for (p in cran_packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}
for (p in bioc_packages) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
  library(p, character.only = TRUE)
}


# -- Configuration
expr_file      <- "data/GSE17114_RMA_all_genes.csv"
meta_file      <- "data/metadata_label.csv"
geneset_file   <- "data/mmc3.xlsx"
core_gene_file <- "results/ml/final_core_genes.csv"
output_dir     <- "results/ssgsea"

group_column  <- "Condition"
case_group    <- "Patient"
control_group <- "Control"
case_color    <- "#990000"
control_color <- "#003366"

cor_method <- "spearman"    # not normal, so pearson is out
p_cutoff   <- 0.05

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
score_file <- file.path(output_dir, "ssgsea_enrichment_scores.csv")

# Short display names for the immune cell types.
cell_name_map <- c(
  "Effector memeory CD4 T cell"    = "CD4 Tem",
  "Effector memeory CD8 T cell"    = "CD8 Tem",
  "Activated B cell"               = "B cell (Act)",
  "Activated CD4 T cell"           = "CD4 T (Act)",
  "Activated CD8 T cell"           = "CD8 T (Act)",
  "Activated dendritic cell"       = "DC (Act)",
  "CD56bright natural killer cell" = "NK (bright)",
  "CD56dim natural killer cell"    = "NK (dim)",
  "Central memory CD4 T cell"      = "CD4 Tcm",
  "Central memory CD8 T cell"      = "CD8 Tcm",
  "Gamma delta T cell"             = "Tgd",
  "Immature B cell"                = "B cell (Imm)",
  "Immature dendritic cell"        = "DC (Imm)",
  "Memory B cell"                  = "B cell (Mem)",
  "Natural killer cell"            = "NK cell",
  "Natural killer T cell"          = "NKT cell",
  "Plasmacytoid dendritic cell"    = "pDC",
  "Regulatory T cell"              = "Treg",
  "T follicular helper cell"       = "Tfh",
  "Type 1 T helper cell"           = "Th1",
  "Type 17 T helper cell"          = "Th17",
  "Type 2 T helper cell"           = "Th2"
)


# --- 1. Read the expression matrix and metadata -------
cat("--- Step 1: reading the expression matrix and metadata ---\n")

# check.names = FALSE keeps the original sample IDs intact.
dat <- read.csv(expr_file, row.names = 1, check.names = FALSE)
colnames(dat) <- gsub("\\..*", "", colnames(dat))   # drop file-name suffixes

metadata <- read.csv(meta_file, check.names = FALSE)
colnames(metadata)[1] <- "SampleID"
rownames(metadata) <- metadata$SampleID


# -- 2. Align the samples
cat("--- Step 2: aligning samples between matrix and metadata ---\n")
common_samples   <- intersect(colnames(dat), rownames(metadata))
dat_aligned      <- dat[, common_samples]
metadata_aligned <- metadata[common_samples, ]

if (!all(colnames(dat_aligned) == rownames(metadata_aligned)))
  stop("Sample alignment failed between the matrix and the metadata.")
cat(sprintf("Aligned on %d shared samples.\n", length(common_samples)))


# 3. Immune cell signature gene sets
cat("--- Step 3: loading the immune signature gene sets ---\n")
geneset <- rio::import(geneset_file, skip = 1)

# GSVA expects a named list: cell type -> character vector of genes.
signature <- geneset %>%
  dplyr::select(Metagene, `Cell type`) %>%
  split(f = .$`Cell type`) %>%
  lapply(FUN = `[[`, "Metagene")


# -- 4. Run ssGSEA
cat("--- Step 4: running ssGSEA ---\n")
gsva_params <- ssgseaParam(exprData = as.matrix(dat_aligned),
                           geneSets = signature)
ssGSEA_matrix <- gsva(param = gsva_params, verbose = FALSE)

write.csv(ssGSEA_matrix, score_file)
cat("Enrichment scores saved to:", score_file, "\n")

ssGSEA_matrix <- read.csv(score_file, header = TRUE, row.names = 1,
                          check.names = FALSE)


# -- 5. Heatmap of the enrichment scores
cat("--- Step 5: enrichment score heatmap ---\n")

# Fixed breaks keep the legend comparable across runs.
my_color  <- colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100)
my_breaks <- seq(-2, 2, length.out = 100)

png(file.path(output_dir, "ssgsea_heatmap.png"),
    width = 10, height = 8, units = "in", res = 300)
pheatmap(ssGSEA_matrix,
         annotation_col = metadata_aligned[, group_column, drop = FALSE],
         scale          = "row",    # z-score per cell type
         show_colnames  = FALSE,
         cluster_cols   = TRUE,
         cluster_rows   = TRUE,
         color          = my_color,
         breaks         = my_breaks,
         legend_breaks  = c(-2, -1, 0, 1, 2),
         border_color   = "white",
         main           = "ssGSEA Enrichment Score Heatmap")
dev.off()


# 6. Long-format table for the group comparisons
df_long <- as.data.frame(t(ssGSEA_matrix)) %>%
  tibble::rownames_to_column("SampleID") %>%
  merge(metadata_aligned, by = "SampleID") %>%
  tidyr::pivot_longer(cols      = -colnames(metadata_aligned),
                      names_to  = "CellType",
                      values_to = "EnrichmentScore") %>%
  mutate(CellType = stringr::str_replace(CellType, ".*Type: ", ""),
         CellType = stringr::str_replace_all(CellType, cell_name_map))

comparison_groups <- list(c(control_group, case_group))


# 6a. All cell types, faceted boxplots
cat("--- Step 6: boxplots across all cell types ---\n")
label_y_pos <- max(df_long$EnrichmentScore) * 1.05

p_all_cells <- ggboxplot(df_long, x = group_column, y = "EnrichmentScore",
                         color = group_column, palette = "jco") +
  geom_jitter(aes(color = .data[[group_column]]), width = 0.2, size = 1.2,
              alpha = 0.8) +
  facet_wrap(~CellType, ncol = 6) +
  stat_compare_means(comparisons = comparison_groups, label = "p.signif",
                     label.y = label_y_pos) +
  labs(title = "Immune Cell Enrichment Between Groups") +
  theme(axis.text.y  = element_text(size = 6),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        strip.text   = element_text(size = 8, face = "bold"))

ggsave(file.path(output_dir, "immune_cell_boxplot_all.png"), p_all_cells,
       width = 14, height = 10, dpi = 300, bg = "white")


# --- 6b. Significantly different cell types -------
# rstatix needs a real formula, so build one from the configured column name.
score_formula <- as.formula(paste("EnrichmentScore ~", group_column))

stat_test <- df_long %>%
  group_by(CellType) %>%
  wilcox_test(score_formula) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

write.csv(stat_test, file.path(output_dir, "immune_cell_wilcox_tests.csv"),
          row.names = FALSE)

sig_cells <- stat_test %>% filter(p < p_cutoff) %>% pull(CellType)

if (length(sig_cells) == 0) {
  cat("No cell type reached p <", p_cutoff, "\n")
} else {
  cat("Cell types with p <", p_cutoff, ":\n")
  print(sig_cells)

  df_sig <- df_long %>% filter(CellType %in% sig_cells)
  n_cols <- min(length(sig_cells), 4)

  p_sig <- ggboxplot(df_sig, x = group_column, y = "EnrichmentScore",
                     color = group_column, palette = "jco", add = "jitter",
                     add.params = list(size = 1.5, alpha = 0.6)) +
    facet_wrap(~CellType, scales = "free", ncol = n_cols) +
    stat_compare_means(comparisons = comparison_groups, label = "p.format") +
    labs(title = sprintf("Significantly Different Immune Cells (P < %.2f)",
                         p_cutoff),
         y = "Enrichment Score", x = "") +
    theme_light() +
    theme(strip.text      = element_text(size = 12, face = "bold"),
          axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          legend.position = "bottom")

  ggsave(file.path(output_dir, "immune_cell_boxplot_significant.png"), p_sig,
         width = ifelse(length(sig_cells) > 2, 10, 6),
         height = ifelse(length(sig_cells) > 4, 8, 5), dpi = 300, bg = "white")
}


# --- 7. Correlation between immune cell types --------
cat("--- Step 7: cell-cell correlation ---\n")

ssGSEA_simplified <- ssGSEA_matrix
rownames(ssGSEA_simplified) <- str_replace_all(rownames(ssGSEA_simplified),
                                               cell_name_map)
cor_matrix <- cor(t(ssGSEA_simplified), method = cor_method)

# Pairwise p-values for the correlation matrix.
cor_pvalues <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      p.mat[i, j] <- p.mat[j, i] <- cor.test(mat[, i], mat[, j], ...)$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p_matrix <- cor_pvalues(t(ssGSEA_simplified))

png(file.path(output_dir, "immune_cell_correlation.png"),
    width = 10, height = 10, units = "in", res = 300)
corrplot.mixed(cor_matrix,
               upper      = "circle",
               upper.col  = rev(brewer.pal(n = 11, name = "RdYlBu")),
               lower      = "number",
               lower.col  = "grey40",
               tl.pos     = "lt",
               tl.col     = "black",
               tl.srt     = 45,
               tl.cex     = 0.7,
               p.mat      = p_matrix,
               insig      = "blank",       # hide non-significant cells
               sig.level  = p_cutoff,
               order      = "hclust",
               diag       = "n",
               mar        = c(0, 0, 1, 0),
               number.cex = 0.6)
dev.off()


# 8. Core gene vs immune cell correlation scatter plots
cat("--- Step 8: gene vs cell correlation scatter plots ---\n")

core_genes <- read.csv(core_gene_file, header = TRUE)[[1]]
valid_genes <- core_genes[core_genes %in% rownames(dat_aligned)]
all_cells <- rownames(ssGSEA_matrix)

#' Plot every gene-cell pair that passes the significance rule.
#'
#' @param mode "any"  keep the pair if either group is significant
#'             "both" keep it only if both groups are significant
plot_gene_cell_correlations <- function(mode = "any") {
  sub_dir <- file.path(output_dir,
                       ifelse(mode == "both",
                              "correlations_significant_both",
                              "correlations_significant_any"))
  if (!dir.exists(sub_dir)) dir.create(sub_dir, recursive = TRUE)

  plot_counter <- 0

  for (target_gene in valid_genes) {
    gene_expression <- as.data.frame(
      t(dat_aligned[target_gene, , drop = FALSE])) %>%
      rownames_to_column("SampleID")
    colnames(gene_expression)[2] <- "GeneExpression"

    for (target_cell in all_cells) {
      cell_scores <- as.data.frame(
        t(ssGSEA_matrix[target_cell, , drop = FALSE])) %>%
        rownames_to_column("SampleID")
      colnames(cell_scores)[2] <- "EnrichmentScore"

      df_scatter <- merge(metadata_aligned, cell_scores, by = "SampleID") %>%
        merge(gene_expression, by = "SampleID")

      cor_stats <- df_scatter %>%
        group_by(.data[[group_column]]) %>%
        rstatix::cor_test(vars = "GeneExpression", vars2 = "EnrichmentScore",
                          method = cor_method)

      keep <- if (mode == "both") {
        nrow(cor_stats) == 2 && all(cor_stats$p < p_cutoff)
      } else {
        any(cor_stats$p < p_cutoff, na.rm = TRUE)
      }
      if (!keep) next

      plot_counter <- plot_counter + 1

      p <- ggscatter(df_scatter, x = "GeneExpression", y = "EnrichmentScore",
                     color = group_column, add = "reg.line", conf.int = TRUE,
                     palette = "jco", size = 2, alpha = 0.6,
                     fullrange = FALSE) +
        stat_cor(aes(color = .data[[group_column]]), method = cor_method,
                 label.y.npc = "top", label.x.npc = "left", size = 3.5) +
        labs(x = paste(target_gene, "Expression (log2)"),
             y = paste(target_cell, "ssGSEA Enrichment Score")) +
        theme_bw() +
        theme(plot.title = element_text(size = 10, face = "bold"))

      clean_cell_name <- str_replace_all(target_cell, "[^[:alnum:]]", "_")
      ggsave(filename = file.path(
               sub_dir,
               paste0(target_gene, "_vs_", clean_cell_name, ".png")),
             plot = p, width = 5.5, height = 4.5, dpi = 300, bg = "white")
    }
  }

  cat(sprintf("Mode '%s': saved %d scatter plots to %s\n",
              mode, plot_counter, sub_dir))
}

plot_gene_cell_correlations(mode = "any")
plot_gene_cell_correlations(mode = "both")


# -- 9. Core gene vs immune cell correlation heatmap
cat("--- Step 9: gene vs cell correlation heatmap ---\n")

core_gene_expr    <- t(dat_aligned[valid_genes, , drop = FALSE])
immune_scores     <- t(ssGSEA_matrix)
colnames(immune_scores) <- stringr::str_replace_all(colnames(immune_scores),
                                                    cell_name_map)

cor_matrix_rect  <- cor(core_gene_expr, immune_scores, method = cor_method)
pval_matrix_rect <- matrix(
  NA, nrow = ncol(core_gene_expr), ncol = ncol(immune_scores),
  dimnames = list(colnames(core_gene_expr), colnames(immune_scores)))

for (gene in colnames(core_gene_expr)) {
  for (cell in colnames(immune_scores)) {
    pval_matrix_rect[gene, cell] <- cor.test(core_gene_expr[, gene],
                                             immune_scores[, cell],
                                             method = cor_method)$p.value
  }
}

sig_matrix_rect <- ifelse(pval_matrix_rect < 0.001, "***",
                   ifelse(pval_matrix_rect < 0.01, "**",
                   ifelse(pval_matrix_rect < 0.05, "*", "")))
sig_matrix_rect[is.na(pval_matrix_rect)] <- ""

write.csv(cor_matrix_rect,
          file.path(output_dir, "core_gene_immune_correlation.csv"))
write.csv(pval_matrix_rect,
          file.path(output_dir, "core_gene_immune_correlation_pvalues.csv"))

p_cor <- pheatmap(cor_matrix_rect,
                  color           = colorRampPalette(
                    rev(brewer.pal(n = 11, name = "RdBu")))(100),
                  border_color    = "grey90",
                  cluster_rows    = FALSE,
                  cluster_cols    = FALSE,
                  angle_col       = "45",
                  fontsize_row    = 10,
                  fontsize_col    = 9,
                  display_numbers = sig_matrix_rect,
                  fontsize_number = 8,
                  silent          = TRUE)

plot_width  <- ncol(cor_matrix_rect) * 0.25 + 3.5
plot_height <- nrow(cor_matrix_rect) * 0.25 + 2

# Drawing through a viewport keeps a margin around the heatmap.
png(file.path(output_dir, "core_gene_immune_correlation_heatmap.png"),
    width = plot_width, height = plot_height, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.95, height = 0.95))
grid.draw(p_cor$gtable)
dev.off()


# 10. Relative abundance stacked bar chart
cat("--- Step 10: relative abundance stacked bar chart ---\n")

abundance_data <- as.data.frame(t(ssGSEA_matrix))
abundance_data[abundance_data < 0] <- 0   # proportions require non-negatives

abundance_data <- abundance_data %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::left_join(metadata_aligned %>% dplyr::rename(sample = SampleID),
                   by = "sample") %>%
  filter(!is.na(.data[[group_column]])) %>%
  dplyr::rename(group = .data[[group_column]]) %>%
  dplyr::select(sample, group,
                any_of(colnames(as.data.frame(t(ssGSEA_matrix)))))

abundance_long <- abundance_data %>%
  tidyr::pivot_longer(cols = -c(sample, group), names_to = "cell",
                      values_to = "score") %>%
  mutate(cell = stringr::str_replace_all(cell, cell_name_map))

# Order the legend by mean abundance so the largest fractions come first.
cell_order <- abundance_long %>%
  group_by(cell) %>%
  summarise(mean_score = mean(score)) %>%
  arrange(desc(mean_score)) %>%
  pull(cell)
abundance_long$cell <- factor(abundance_long$cell, levels = cell_order)

color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(
  length(unique(abundance_long$cell)))

p_stacked_bar <- ggplot(abundance_long,
                        aes(x = sample, y = score, fill = cell)) +
  geom_col(position = "fill", width = 1.0) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = color_palette) +
  labs(x = NULL, y = "Relative Abundance", fill = "Immune Cell Type") +
  theme_classic() +
  theme(axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.title.y     = element_text(size = 14, face = "bold"),
        axis.text.y      = element_text(size = 12),
        legend.title     = element_text(size = 14, face = "bold"),
        legend.text      = element_text(size = 10),
        strip.background = element_rect(fill = "gray90", color = "black"),
        strip.text       = element_text(size = 14, face = "bold",
                                        color = "black"),
        axis.line        = element_line(colour = "black"),
        panel.spacing    = unit(0, "lines"))

ggsave(file.path(output_dir, "immune_abundance_stacked_bar.png"),
       p_stacked_bar, width = 12, height = 7, dpi = 300, bg = "white")


# --- 11. Differential boxplot with BH-adjusted p-values ---
cat("--- Step 11: differential boxplot ---\n")

plot_data <- abundance_long %>%
  filter(group %in% c(case_group, control_group))

stat_result <- plot_data %>%
  group_by(cell) %>%
  wilcox_test(score ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col     = "p.adj",
                   cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                   symbols   = c("***", "**", "*", "ns"))

# Most significant cell types first.
plot_data$cell <- factor(plot_data$cell,
                         levels = stat_result %>% arrange(p.adj) %>% pull(cell))
stat_result <- stat_result %>% add_xy_position(x = "cell", dodge = 0.8)

p_boxplot <- ggplot(plot_data, aes(x = cell, y = score)) +
  geom_boxplot(aes(fill = group), position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_jitter(aes(color = group),
              position = position_jitterdodge(0.8, jitter.width = 0.2),
              size = 1, alpha = 0.6) +
  stat_pvalue_manual(stat_result, label = "p.adj.signif", tip.length = 0.01,
                     hide.ns = TRUE) +
  scale_fill_manual(values = c(control_color, case_color)) +
  scale_color_manual(values = c(control_color, case_color)) +
  labs(x = NULL, y = "ssGSEA Enrichment Score", fill = "Group",
       color = "Group") +
  theme_classic() +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 12,
                                       color = "black"),
        axis.text.y     = element_text(size = 12, color = "black"),
        axis.title.y    = element_text(size = 14, face = "bold"),
        legend.position = "top",
        legend.title    = element_text(size = 14, face = "bold"),
        legend.text     = element_text(size = 12))

ggsave(file.path(output_dir, "immune_differential_boxplot.png"), p_boxplot,
       width = 14, height = 8, dpi = 300, bg = "white")

cat("Step 10 complete. Results in:", output_dir, "\n")
