# limma DEG + WGCNA on the GSE70403 matrix from step 02.
# DEG and WGCNA hub genes are pooled (union), then intersected with the
# BPA-Behcet target list to give the final candidate set.

# --- 0. Environment ------
rm(list = ls())
options(stringsAsFactors = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("limma", "WGCNA", "data.table", "dplyr", "tidyr", "ggplot2",
              "ggrepel", "pheatmap", "ggvenn", "ggpubr", "gplots", "UpSetR")
bioc_packages <- c("limma", "WGCNA")

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    if (p %in% bioc_packages) BiocManager::install(p, ask = FALSE)
    else install.packages(p)
  }
  library(p, character.only = TRUE)
}
allowWGCNAThreads()


# Configuration
expr_file  <- "results/GSE70403_expression_samples_x_genes.csv"
bpa_file   <- "data/bpa_behcet_targets.csv"
output_dir <- "results/DEG_WGCNA"

# DEG thresholds
adj_p_cutoff  <- 0.05
logfc_cutoff  <- 0.585      # ~1.5 fold change

# WGCNA parameters
n_genes_to_keep  <- 5000    # 8000 works too but TOM gets slow on a laptop
# n_genes_to_keep <- as.integer(ncol(datExpr) * 0.25)
min_module_size  <- 30
merge_cut_height <- 0.25
kme_cutoff       <- 0.9     # 0.85 gives a much bigger pool
fallback_power   <- 6
r2_target        <- 0.85

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat("Results will be written to:", file.path(getwd(), output_dir), "\n")


# --- 1. Load the expression matrix (samples x genes) ----
message("--- PART 1: loading the preprocessed matrix ---")
if (!file.exists(expr_file)) stop(paste("Expression file not found:", expr_file))

datExpr_df <- fread(expr_file) %>% as.data.frame()
rownames(datExpr_df) <- datExpr_df[, 1]          # first column holds sample IDs
datExpr <- datExpr_df[, -1]
datExpr <- as.data.frame(apply(datExpr, 2, as.numeric))
rownames(datExpr) <- rownames(datExpr_df)

cat(sprintf("Loaded %d samples and %d genes.\n", nrow(datExpr), ncol(datExpr)))


# --- 2. Differential expression with limma ----
message("--- PART 2: limma differential expression (single-group design) ---")

limma_input <- t(datExpr)                        # limma expects genes x samples
design      <- model.matrix(~ 1, data = datExpr)
fit         <- lmFit(limma_input, design)
fit         <- eBayes(fit)

all_deg_results      <- topTable(fit, coef = 1, number = Inf,
                                 adjust.method = "BH")
all_deg_results$GENE <- rownames(all_deg_results)
write.csv(all_deg_results, file.path(output_dir, "DEG_results_all_genes.csv"),
          row.names = FALSE)

deg_genes_df  <- all_deg_results %>%
  filter(adj.P.Val < adj_p_cutoff & abs(logFC) > logfc_cutoff)
deg_gene_list <- unique(deg_genes_df$GENE)

write.csv(data.frame(Gene = deg_gene_list),
          file.path(output_dir, "list_DEG_significant.csv"), row.names = FALSE)
cat("Significant DEGs:", length(deg_gene_list), "\n")


# 2b. Volcano plot
all_deg_results <- all_deg_results %>%
  mutate(Sig = case_when(
    adj.P.Val < adj_p_cutoff & logFC >  logfc_cutoff ~ "Up",
    adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff ~ "Down",
    TRUE ~ "Not"
  )) %>%
  mutate(Sig = factor(Sig, levels = c("Down", "Not", "Up")))

volcano_plot <- ggplot(all_deg_results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "green", "Not" = "grey")) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), lty = "dashed") +
  geom_hline(yintercept = -log10(adj_p_cutoff), lty = "dashed") +
  labs(x = "log2(Fold Change)", y = "-log10(Adjusted P-value)") +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank())

ggsave(file.path(output_dir, "DEG_volcano_plot.png"), volcano_plot,
       width = 8, height = 7, units = "in", dpi = 300)


# -- 2c. Heatmap of significant DEGs
if (length(deg_gene_list) < 2) {
  cat("Fewer than two significant DEGs; skipping the heatmap.\n")
} else {
  deg_expr_matrix <- datExpr[, deg_gene_list, drop = FALSE]
  matrix_to_plot  <- t(deg_expr_matrix)          # genes x samples

  cal_z_score <- function(x) {
    if (is.na(sd(x, na.rm = TRUE)) || sd(x, na.rm = TRUE) == 0)
      return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }

  plot_matrix_scaled <- t(apply(matrix_to_plot, 1, cal_z_score))
  colnames(plot_matrix_scaled) <- rownames(deg_expr_matrix)
  rownames(plot_matrix_scaled) <- colnames(deg_expr_matrix)

  # Clip extreme z-scores so the colour scale stays readable.
  z_limit <- 2.5
  plot_matrix_scaled[plot_matrix_scaled >  z_limit] <-  z_limit
  plot_matrix_scaled[plot_matrix_scaled < -z_limit] <- -z_limit

  my_colors <- colorRampPalette(c("steelblue", "white", "firebrick"))(100)
  my_breaks <- seq(-z_limit, z_limit, length.out = 101)

  pheatmap_plot <- pheatmap(
    plot_matrix_scaled,
    scale                    = "none",   # already scaled above
    color                    = my_colors,
    breaks                   = my_breaks,
    fontsize                 = 8,
    fontsize_row             = 6,
    fontsize_col             = 8,
    show_rownames            = (nrow(plot_matrix_scaled) < 50),
    show_colnames            = FALSE,
    border_color             = NA,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method        = "ward.D2"
  )
  ggsave(file.path(output_dir, "DEG_heatmap.png"), pheatmap_plot$gtable,
         width = 8, height = 8, units = "in", dpi = 300, bg = "white")
}


# -- 3a. Clean sample names
message("--- PART 3: WGCNA ---")

new_names <- gsub("_.*$", "", rownames(datExpr))  # keep the GSM accession only
if (any(duplicated(new_names)))
  warning("Cleaning the sample names produced duplicates; please check.")
rownames(datExpr) <- new_names


# --- 3b. Keep the most variable genes ---
# The full gene set makes the TOM computation very slow, so only the most
# variable genes are carried into the network.
gene_vars <- apply(datExpr, 2, var)
if (ncol(datExpr) > n_genes_to_keep) {
  message(sprintf("Filtering %d genes down to the %d most variable ...",
                  ncol(datExpr), n_genes_to_keep))
  genes_to_keep <- names(sort(gene_vars, decreasing = TRUE))[1:n_genes_to_keep]
  datExpr <- datExpr[, genes_to_keep]
}


# -- 3c. Sample and gene QC
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

sampleTree <- hclust(dist(datExpr), method = "average")
png(file.path(output_dir, "WGCNA_sample_clustering.png"),
    width = 10, height = 6, units = "in", res = 300)
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "", xlab = "", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.6)
dev.off()


# --- 3d. Soft threshold -----
powers     <- 1:20
sft        <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
r2vals     <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
sel_powers <- sft$fitIndices[, 1][which(r2vals >= r2_target)]
softPower  <- if (length(sel_powers) > 0) min(sel_powers) else fallback_power
if (length(sel_powers) == 0)
  message(sprintf("No power reached R^2 >= %.2f; falling back to %d.",
                  r2_target, fallback_power))
cat("Selected soft threshold power:", softPower, "\n")

png(file.path(output_dir, "WGCNA_soft_threshold.png"),
    width = 10, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))
cex1 <- 0.8
plot(sft$fitIndices[, 1], r2vals, type = "n",
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale independence")
text(sft$fitIndices[, 1], r2vals, labels = powers, cex = cex1, col = "red")
abline(h = r2_target, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n",
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers,
     cex = cex1, col = "red")
dev.off()


# 3e. Network construction and module detection
adjMat   <- adjacency(datExpr, power = softPower, type = "signed")
TOM      <- TOMsimilarity(adjMat, TOMType = "signed")
dissTOM  <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

dynamicMods   <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = min_module_size)
dynamicColors <- labels2colors(dynamicMods)

mergeRes    <- mergeCloseModules(datExpr, dynamicColors,
                                 cutHeight = merge_cut_height, verbose = 3)
finalColors <- mergeRes$colors
finalMEs    <- mergeRes$newMEs

png(file.path(output_dir, "WGCNA_dendrogram_modules.png"),
    width = 12, height = 8, units = "in", res = 300)
plotDendroAndColors(geneTree, finalColors, "Merged dynamic",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()


# -- 3f. Gene significance vs module membership
message("--- PART 3b: GS vs MM across modules ---")
logfc_col <- grep("logFC", colnames(all_deg_results),
                  ignore.case = TRUE, value = TRUE)[1]

if (!is.na(logfc_col)) {
  common_genes <- intersect(colnames(datExpr), rownames(all_deg_results))
  if (length(common_genes) == 0)
    stop("No overlap between datExpr columns and the DEG table rows.")

  GS_vector        <- abs(all_deg_results[common_genes, logfc_col])
  names(GS_vector) <- common_genes

  plot_data_all <- data.frame()
  all_modules   <- setdiff(unique(finalColors), "grey")   # grey = unassigned

  for (mod in all_modules) {
    mod_genes   <- colnames(datExpr)[finalColors == mod]
    valid_genes <- intersect(mod_genes, common_genes)

    if (length(valid_genes) > 10) {
      me_col   <- paste0("ME", mod)
      mm_vals  <- as.numeric(cor(datExpr[, valid_genes],
                                 finalMEs[, me_col], use = "p"))
      plot_data_all <- rbind(plot_data_all, data.frame(
        Gene   = valid_genes,
        GS     = GS_vector[valid_genes],
        MM     = abs(mm_vals),   # sign is irrelevant for membership strength
        Module = mod
      ))
    }
  }

  p_all <- ggplot(plot_data_all, aes(x = MM, y = GS)) +
    geom_point(aes(color = Module), alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", color = "black", se = TRUE, size = 0.5) +
    stat_cor(method = "pearson", label.x = 0.05,
             label.y.npc = "top", size = 2.8) +
    facet_wrap(~Module, scales = "free") +
    scale_color_identity() +
    theme_bw() +
    labs(title = "Gene Significance vs. Module Membership (All Modules)",
         x = "Module Membership (kME)",
         y = "Gene Significance (|logFC|)") +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          strip.background = element_rect(fill = "white"))

  ggsave(file.path(output_dir, "WGCNA_GS_vs_MM_all_modules.png"), p_all,
         width = 12, height = 10, dpi = 300)
} else {
  warning("No 'logFC' column found in the DEG table; skipping GS vs MM.")
}


# 3g. Module eigengene relationships
png(file.path(output_dir, "WGCNA_eigengene_network.png"),
    width = 8, height = 8, units = "in", res = 300)
par(cex = 1.0)
plotEigengeneNetworks(finalMEs, "", marDendro = c(0, 4, 1, 2),
                      marHeatmap = c(3, 4, 1, 2), cex.lab = 0.8,
                      xLabelsAngle = 90)
dev.off()

# Correlation heatmap of the module eigengenes (grey module excluded).
plotMEs <- if ("MEgrey" %in% colnames(finalMEs))
  finalMEs[, colnames(finalMEs) != "MEgrey"] else finalMEs
ME_cor <- cor(plotMEs, use = "p")

mod_colors      <- substring(colnames(ME_cor), 3)   # "MEblue" -> "blue"
annotation_row  <- data.frame(Module = mod_colors)
rownames(annotation_row) <- colnames(ME_cor)
ann_colors      <- list(Module = setNames(mod_colors, mod_colors))

png(file.path(output_dir, "WGCNA_eigengene_correlation.png"),
    width = 14, height = 8, units = "in", res = 300)
pheatmap(ME_cor,
         color             = colorRampPalette(c("navy", "white",
                                                "firebrick3"))(100),
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         border_color      = "grey60",
         number_color      = "black",
         fontsize          = 12,
         fontsize_number   = 10,
         display_numbers   = TRUE,
         number_format     = "%.2f",
         annotation_row    = annotation_row,
         annotation_col    = annotation_row,
         annotation_colors = ann_colors,
         annotation_legend = FALSE,
         main              = "")
dev.off()


# --- 3h. Hub genes (kME above the cutoff) --------
message(sprintf("--- Selecting WGCNA hub genes (kME > %.2f) ---", kme_cutoff))

module_info          <- data.frame(Gene = colnames(datExpr),
                                   ModuleColor = finalColors)
geneModuleMembership <- as.data.frame(cor(datExpr, finalMEs, use = "p"))
colnames(geneModuleMembership) <- paste0("kME_",
                                         substring(colnames(finalMEs), 3))

kME_table_long <- cbind(module_info, geneModuleMembership) %>%
  pivot_longer(cols = starts_with("kME_"),
               names_to = "kME_Module_Name", values_to = "kME_Value")

wgcna_hub_genes_list <- kME_table_long %>%
  filter(ModuleColor != "grey",
         paste0("kME_", ModuleColor) == kME_Module_Name,
         kME_Value > kme_cutoff) %>%
  pull(Gene) %>%
  unique()

write.csv(data.frame(Gene = wgcna_hub_genes_list),
          file.path(output_dir, "list_WGCNA_hub_genes.csv"), row.names = FALSE)
cat("WGCNA hub genes:", length(wgcna_hub_genes_list), "\n")


# 4. Candidate pool: union of DEGs and hub genes
message("--- PART 4: candidate gene pool ---")
union_genes_list <- union(deg_gene_list, wgcna_hub_genes_list)
cat("  DEGs:            ", length(deg_gene_list), "\n")
cat("  WGCNA hub genes: ", length(wgcna_hub_genes_list), "\n")
cat("  Candidate union: ", length(union_genes_list), "\n")

write.csv(data.frame(Gene = union_genes_list),
          file.path(output_dir, "list_candidate_genes_union.csv"),
          row.names = FALSE)

venn_plot <- ggvenn(list(DEG = deg_gene_list, WGCNA = wgcna_hub_genes_list),
                    fill_color = c("#FF8C69", "#87CEFA"), stroke_size = 0.5,
                    set_name_size = 6, text_size = 5)
ggsave(file.path(output_dir, "candidate_genes_venn.png"), venn_plot,
       width = 7, height = 7, units = "in", dpi = 300, bg = "white")


# 4b. UpSet plot: DEGs vs individual WGCNA modules
gene_color_map <- data.frame(Gene = colnames(datExpr), Module = finalColors)
hub_info       <- gene_color_map %>% filter(Gene %in% wgcna_hub_genes_list)
wgcna_sets     <- split(hub_info$Gene, hub_info$Module)
names(wgcna_sets) <- paste0("Hub_", names(wgcna_sets))

list_input <- c(list(DEG = deg_gene_list), wgcna_sets)
list_input <- list_input[sapply(list_input, length) >= 5]

# Colour each bar with its own module colour; DEGs get a distinct red.
get_color_for_set <- function(set_name) {
  if (set_name == "DEG") "#E41A1C" else gsub("Hub_", "", set_name)
}
color_queries <- lapply(names(list_input), function(set_name) {
  list(query  = intersects,
       params = list(set_name),
       color  = get_color_for_set(set_name),
       active = TRUE)
})

png(file.path(output_dir, "candidate_genes_upset.png"),
    width = 12, height = 8, units = "in", res = 300)
upset(fromList(list_input),
      nsets            = 20,
      order.by         = "freq",
      sets.bar.color   = "#4682B4",
      main.bar.color   = "black",
      matrix.color     = "#36454F",
      shade.color      = "#F5F5F5",
      text.scale       = c(1.4, 1.3, 1.4, 1.3, 1.3, 1.1),
      queries          = color_queries,
      mainbar.y.label  = "intersection size",
      sets.x.label     = "set size",
      query.legend     = "none")
dev.off()


# -- 5. Expression matrix of the candidate pool
message("--- PART 5: saving the candidate expression matrix ---")
genes_in_expr <- intersect(union_genes_list, colnames(datExpr))
if (length(genes_in_expr) == 0)
  stop("None of the candidate genes are present in datExpr.")

candidate_expr_matrix <- datExpr[, genes_in_expr, drop = FALSE]
final_df_to_save <- as.data.frame(candidate_expr_matrix) %>%
  mutate(SampleID = rownames(candidate_expr_matrix)) %>%
  select(SampleID, everything())

write.csv(final_df_to_save,
          file.path(output_dir, "candidate_genes_union_matrix.csv"),
          row.names = FALSE)


# -- 6. Intersection with the BPA-Behcet target list
message("--- PART 6: intersecting the candidate pool with the BPA targets ---")
if (!file.exists(bpa_file))
  stop(paste("BPA target file not found:", bpa_file))

candidate_genes <- unique(
  fread(file.path(output_dir, "list_candidate_genes_union.csv"))[[1]]
)
bpa_targets <- unique(fread(bpa_file)[[1]])   # first column holds the symbols

message(paste("  Candidate genes:   ", length(candidate_genes)))
message(paste("  BPA-Behcet targets:", length(bpa_targets)))

final_overlap_genes <- intersect(candidate_genes, bpa_targets)
message(paste("  Final intersection:", length(final_overlap_genes)))

if (length(final_overlap_genes) == 0) {
  warning("Empty intersection; check that both files use gene symbols.")
} else {
  print(head(final_overlap_genes))
}

write.csv(data.frame(Gene = final_overlap_genes),
          file.path(output_dir, "final_genes_bpa_intersection.csv"),
          row.names = FALSE)

venn_plot_final <- ggvenn(
  list("Clinical Signatures"   = candidate_genes,
       "Toxicological Targets" = bpa_targets),
  fill_color = c("#FF8C69", "#87CEFA"), stroke_size = 0.5,
  set_name_size = 5, text_size = 5
)
ggsave(file.path(output_dir, "final_bpa_intersection_venn.png"),
       venn_plot_final, width = 7, height = 7, units = "in",
       dpi = 300, bg = "white")

save.image(file.path(output_dir, "deg_wgcna_session.RData"))
message("Step 03 complete.")
