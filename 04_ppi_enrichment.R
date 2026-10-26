# KEGG / GO enrichment + STRING PPI for the final candidate genes.
# STRING can run offline from the downloaded flat files (default) or through
# the STRINGdb web service, see use_offline_string below.

# --- 0. Packages --------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_packages <- c("igraph", "ggraph", "dplyr", "tidyr", "ggplot2",
                   "ggrepel", "readr", "data.table", "scales")
bioc_packages <- c("clusterProfiler", "org.Hs.eg.db", "STRINGdb")

for (p in cran_packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}
for (p in bioc_packages) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
  library(p, character.only = TRUE)
}


# -- Configuration
gene_file   <- "results/DEG_WGCNA/final_genes_bpa_intersection.csv"
deg_file    <- "results/DEG_WGCNA/DEG_results_all_genes.csv"
output_dir  <- "results/enrichment_ppi"

# STRING settings. Set use_offline_string to TRUE to read the downloaded
# STRING flat files instead of querying the web service.
string_dir        <- "data/string_data"
string_version    <- "11.5"
string_species    <- 9606              # Homo sapiens
score_threshold   <- 400               # medium confidence
use_offline_string <- TRUE
alias_file <- file.path(string_dir,
                        paste0(string_species, ".protein.aliases.v",
                               string_version, ".txt.gz"))
links_file <- file.path(string_dir,
                        paste0(string_species, ".protein.links.v",
                               string_version, ".txt.gz"))

top_n_hubs <- 8                        # nodes placed on the inner ring

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat("Results will be written to:", output_dir, "\n")


# --- 1. Candidate genes ---
if (!file.exists(gene_file)) stop(paste("Gene list not found:", gene_file))
genes_df <- read.csv(gene_file, stringsAsFactors = FALSE)
genes    <- unique(genes_df$Gene)
cat("Candidate genes:", length(genes), "\n")


# 2. Up/down status from the DEG table
if (!file.exists(deg_file)) stop(paste("DEG table not found:", deg_file))
deg_results <- read.csv(deg_file, stringsAsFactors = FALSE)
if (!all(c("GENE", "logFC") %in% colnames(deg_results)))
  stop("The DEG table must contain 'GENE' and 'logFC' columns.")

gene_status_df <- deg_results %>%
  mutate(status = ifelse(logFC > 0, "Upregulated", "Downregulated")) %>%
  select(GENE, logFC, adj.P.Val, status) %>%
  distinct(GENE, .keep_all = TRUE)

cat(length(intersect(genes, gene_status_df$GENE)),
    "candidate genes matched a DEG status.\n")


# 3. Symbol -> Entrez ID
gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
cat("Genes with an Entrez ID:", nrow(gene_entrez), "\n")
write.csv(gene_entrez, file.path(output_dir, "gene_symbol_to_entrezid.csv"),
          row.names = FALSE)


# --- 4. KEGG enrichment ----
options(timeout = 600)
kegg_enrich <- enrichKEGG(gene          = gene_entrez$ENTREZID,
                          organism      = "hsa",
                          pAdjustMethod = "BH",
                          qvalueCutoff  = 0.05)

if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
  write.csv(kegg_enrich@result,
            file.path(output_dir, "kegg_enrichment_results.csv"),
            row.names = FALSE)

  p_kegg_dot <- dotplot(kegg_enrich, showCategory = 20, font.size = 10,
                        color = "p.adjust") +
    scale_color_gradientn(
      colors = c("#FF0000", "#FF4500", "#FFA500", "#FFD700", "#ADFF2F"),
      limits = c(0, 0.05), name = "p.adjust") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "GeneRatio")

  ggsave(file.path(output_dir, "kegg_dotplot.png"), p_kegg_dot,
         width = 10, height = 8, dpi = 300)
} else {
  cat("No significantly enriched KEGG pathway.\n")
}


# -- 5. GO enrichment (BP, MF, CC)
go_enrich <- enrichGO(gene          = gene_entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

if (!is.null(go_enrich) && nrow(go_enrich@result) > 0) {
  write.csv(go_enrich@result,
            file.path(output_dir, "go_enrichment_results.csv"),
            row.names = FALSE)

  # Top 10 terms per ontology, ranked by adjusted p-value.
  go_df_top10 <- as.data.frame(go_enrich) %>%
    filter(p.adjust < 0.05) %>%
    group_by(ONTOLOGY) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    ungroup()

  p_go_bar <- ggplot(go_df_top10,
                     aes(x = reorder(Description, Count), y = Count,
                         fill = p.adjust)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient2(
      low      = "#DC605C",
      mid      = "#666699",
      high     = "#9999CC",
      midpoint = 0.05,
      name     = "p.adjust",
      limits   = c(0, max(go_df_top10$p.adjust, 0.08)),
      oob      = scales::squish) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    theme_bw() +
    theme(plot.title       = element_text(hjust = 0.5),
          axis.title.y     = element_blank(),
          strip.background = element_rect(fill = "grey85", color = "black"),
          strip.text.y     = element_text(angle = 0)) +
    labs(x = "Count")

  ggsave(file.path(output_dir, "go_barplot.png"), p_go_bar,
         width = 12, height = 10, dpi = 300)
} else {
  cat("No significantly enriched GO term.\n")
}


# -- 6. PPI edge list from STRING
options(timeout = 3600)
ppi_edge_file <- file.path(output_dir, "ppi_edge_list.csv")

if (use_offline_string) {
  # Offline mode: read the downloaded STRING flat files directly.
  message("Building the PPI network from local STRING files ...")
  if (!file.exists(alias_file) || !file.exists(links_file))
    stop("STRING flat files not found; check 'string_dir'.")

  aliases_df <- fread(alias_file, select = c(1, 2),
                      col.names = c("STRING_id", "alias"), header = FALSE)

  mapped_df <- aliases_df %>%
    filter(alias %in% genes) %>%
    distinct(alias, .keep_all = TRUE) %>%   # one STRING id per gene
    rename(GENE = alias)
  cat("Mapped genes:", nrow(mapped_df), "\n")

  links_df <- fread(links_file,
                    col.names = c("protein1", "protein2", "combined_score"))

  target_ids <- mapped_df$STRING_id
  ppi_manual <- links_df %>%
    filter(protein1 %in% target_ids, protein2 %in% target_ids,
           combined_score >= score_threshold)
  cat("Interactions found:", nrow(ppi_manual), "\n")

  id_to_symbol <- setNames(mapped_df$GENE, mapped_df$STRING_id)
  ppi_readable <- ppi_manual %>%
    mutate(from = id_to_symbol[protein1],
           to   = id_to_symbol[protein2]) %>%
    select(from, to, combined_score)

} else {
  # Online mode: let STRINGdb handle the mapping and the API call.
  message("Querying STRING through STRINGdb ...")
  string_db <- STRINGdb$new(version         = string_version,
                            species         = string_species,
                            score_threshold = score_threshold,
                            input_directory = string_dir)

  map_step1 <- string_db$map(data.frame(GENE = genes), "GENE",
                             removeUnmappedRows = FALSE)

  unmapped <- map_step1$GENE[is.na(map_step1$STRING_id)]
  if (length(unmapped) > 0)
    cat("Unmapped symbols:", paste(unmapped, collapse = ", "), "\n")

  # Keep the first hit per gene and drop anything that failed to map.
  mapped_genes_final <- map_step1 %>%
    distinct(GENE, .keep_all = TRUE) %>%
    filter(!is.na(STRING_id))
  cat("Genes used for the network:", nrow(mapped_genes_final), "\n")

  ppi <- string_db$get_interactions(mapped_genes_final$STRING_id)
  cat("Interactions found:", nrow(ppi), "\n")

  ppi_readable <- ppi
  ppi_readable$from <- mapped_genes_final$GENE[
    match(ppi$from, mapped_genes_final$STRING_id)]
  ppi_readable$to <- mapped_genes_final$GENE[
    match(ppi$to, mapped_genes_final$STRING_id)]
}

write.csv(ppi_readable, ppi_edge_file, row.names = FALSE)
cat("PPI edge list saved to:", ppi_edge_file, "\n")


# -- 7. Network layout
ppi_data <- read.csv(ppi_edge_file, stringsAsFactors = FALSE)
g <- graph_from_data_frame(ppi_data, directed = FALSE)

V(g)$degree <- degree(g)
V(g)$status <- gene_status_df$status[match(V(g)$name, gene_status_df$GENE)]
V(g)$status[is.na(V(g)$status)] <- "Unknown"

# Concentric layout: the highest-degree nodes sit on a small inner circle,
# everything else on a wider outer circle.
hubs   <- names(sort(degree(g), decreasing = TRUE)[1:min(top_n_hubs, vcount(g))])
others <- setdiff(V(g)$name, hubs)
V(g)$node_type <- ifelse(V(g)$name %in% hubs, "Hub Genes", "Other Genes")

layout_matrix <- matrix(0, nrow = vcount(g), ncol = 2)
rownames(layout_matrix) <- V(g)$name
if (length(hubs) > 0)
  layout_matrix[hubs, ] <-
    layout_in_circle(make_empty_graph(length(hubs))) * 0.25
if (length(others) > 0)
  layout_matrix[others, ] <-
    layout_in_circle(make_empty_graph(length(others))) * 0.5
layout_final <- layout_matrix[V(g)$name, ]

# Shared styling for the three network figures below.
base_network_theme <- function(p, legend_title = "") {
  p +
    geom_edge_link(color = "grey70", width = 0.5, alpha = 0.6) +
    geom_node_text(aes(label = name), color = "black",
                   fontface = "bold", size = 3.5) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    coord_fixed(ratio = 1) +
    theme_void() +
    theme(legend.position = "right",
          legend.title    = element_text(face = "bold"),
          legend.text     = element_text(size = 10),
          plot.margin     = margin(10, 10, 10, 10)) +
    labs(color = legend_title)
}


# 7a. Nodes coloured by up/down regulation
custom_colors <- c("Upregulated"   = "#d35400",
                   "Downregulated" = "#2980b9",
                   "Unknown"       = "grey70")

p_status <- base_network_theme(
  ggraph(g, layout = layout_final) +
    geom_node_point(aes(color = status), size = 22)
) +
  scale_color_manual(values = custom_colors) +
  guides(color = guide_legend(override.aes = list(size = 6)))

ggsave(file.path(output_dir, "ppi_network_by_status.png"), p_status,
       width = 10, height = 8, dpi = 300, bg = "white")


# -- 7b. Nodes coloured by hub / non-hub
color_scheme <- c("Hub Genes" = "#d35400", "Other Genes" = "#f39c12")

p_hub <- base_network_theme(
  ggraph(g, layout = layout_final) +
    geom_node_point(aes(color = node_type), size = 22),
  legend_title = "Gene Type"
) +
  scale_color_manual(values = color_scheme) +
  guides(color = guide_legend(override.aes = list(size = 6)))

ggsave(file.path(output_dir, "ppi_network_by_hub.png"), p_hub,
       width = 10, height = 8, dpi = 300, bg = "white")


# --- 7c. Nodes coloured by a degree gradient --------
grad_colors   <- colorRampPalette(c("#ffeda0", "#f39c12", "#d35400"))(100)
V(g)$color    <- grad_colors[cut(V(g)$degree, 100, labels = FALSE)]

p_degree <- ggraph(g, layout = layout_final) +
  geom_edge_link(color = "grey70", width = 0.5, alpha = 0.6) +
  geom_node_point(aes(color = I(V(g)$color)), size = 22) +
  geom_node_text(aes(label = name), color = "black",
                 fontface = "bold", size = 3.5) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(legend.position  = "none",
        plot.margin      = margin(0, 0, 0, 0),
        panel.background = element_blank(),
        plot.background  = element_blank())

ggsave(file.path(output_dir, "ppi_network_by_degree.png"), p_degree,
       width = 10, height = 8, dpi = 300, bg = "white", limitsize = FALSE)

cat("Degree range:", min(V(g)$degree), "-", max(V(g)$degree), "\n")
cat("Number of nodes:", vcount(g), "\n")
cat("Top hub genes by degree:\n")
print(sort(setNames(V(g)$degree, V(g)$name),
           decreasing = TRUE)[1:min(top_n_hubs, vcount(g))])

message("Step 04 complete.")
