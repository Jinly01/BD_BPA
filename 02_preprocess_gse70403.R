# GSE70403: raw Agilent two-colour scans -> gene level M-value matrix.
# Needs the unpacked RAW.tar in data/GSE70403_RAW/ plus the series matrix
# (for the sample metadata) and the GPL4133 annotation.

# --- 0. Packages -----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma")
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")

library(GEOquery)
library(limma)
library(data.table)

# --- Configuration (edit these paths for your machine) ----
raw_files_directory <- "data/GSE70403_RAW/"
matrix_file_path    <- "data/GSE70403_series_matrix.txt.gz"
annotation_file     <- file.path(raw_files_directory,
                                 "GPL4133_old_annotations.txt.gz")
output_dir          <- "results"
output_file         <- file.path(output_dir,
                                 "GSE70403_expression_samples_x_genes.csv")

gene_symbol_col <- "GENE_SYMBOL"
probe_name_col  <- "NAME"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(raw_files_directory))
  stop("RAW directory not found; check 'raw_files_directory'.")
if (!file.exists(matrix_file_path))
  stop("Series matrix file not found; check 'matrix_file_path'.")


# -- 1. Sample metadata from the series matrix
message("Loading sample metadata from the series matrix ...")
gse_meta   <- getGEO(filename = matrix_file_path, GSEMatrix = TRUE,
                     getGPL = FALSE)
pheno_data <- pData(gse_meta)

# limma needs a targets frame: one row per array, pointing at the raw file.
targets <- data.frame(
  FileName = basename(as.character(pheno_data$supplementary_file)),
  Group    = "BD"   # all 41 arrays in this series are BD patients
)
rownames(targets) <- pheno_data$geo_accession
print(head(targets))


# --- 2. Load and normalise the raw intensities --
message("Reading raw Agilent files ...")
raw_data <- read.maimages(targets,
                          path   = raw_files_directory,
                          source = "agilent")

message("Background correction (normexp) ...")
y <- backgroundCorrect(raw_data, method = "normexp", offset = 16)

# Within-array loess normalisation corrects the red/green dye bias and
# produces M = log2(R/G) and A = 0.5 * log2(R*G).
message("Within-array normalisation (loess) ...")
MA <- normalizeWithinArrays(y, method = "loess")

# Between-array quantile normalisation makes M-values comparable across arrays.
message("Between-array normalisation (quantile) ...")
MA_norm <- normalizeBetweenArrays(MA, method = "quantile")


# -- 3. Platform annotation
message("Loading the platform annotation ...")
if (!file.exists(annotation_file))
  stop(paste("Annotation file not found:", annotation_file))

annot_file <- as.data.frame(fread(annotation_file))

if (!probe_name_col %in% colnames(annot_file))
  stop(paste("Probe column not found:", probe_name_col))
if (!gene_symbol_col %in% colnames(annot_file))
  stop(paste("Gene symbol column not found:", gene_symbol_col))


# --- 4. Map probes to gene symbols and aggregate --
message("Mapping probes to gene symbols ...")
probe_names <- MA_norm$genes$ProbeName

annot_clean <- annot_file[!is.na(annot_file[[gene_symbol_col]]) &
                            annot_file[[gene_symbol_col]] != "", ]

# Lookup table: probe name -> gene symbol.
lookup_map   <- setNames(annot_clean[[gene_symbol_col]],
                         annot_clean[[probe_name_col]])
gene_symbols <- lookup_map[probe_names]

valid_probes          <- !is.na(gene_symbols) & gene_symbols != ""
MA_norm_filtered      <- MA_norm[valid_probes, ]
gene_symbols_filtered <- gene_symbols[valid_probes]

message("Averaging replicate probes per gene ...")
MA_gene <- avereps(MA_norm_filtered, ID = gene_symbols_filtered)

message(paste("Probes before aggregation:", nrow(MA_norm_filtered)))
message(paste("Genes after aggregation: ", nrow(MA_gene)))
if (nrow(MA_gene) == 0) stop("Aggregation produced zero genes; check the mapping.")


# --- 5. Save the samples x genes matrix ------
# MA_gene$M is genes x samples, so transpose it for the downstream steps.
final_expression_matrix <- t(MA_gene$M)
message("Final matrix dimensions (samples x genes):")
print(dim(final_expression_matrix))

write.csv(final_expression_matrix, file = output_file, row.names = TRUE)
message(paste("Saved:", normalizePath(output_file)))
