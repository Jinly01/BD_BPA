"""
Series matrix -> gene symbol matrix.

Probe IDs are mapped through the GPL annotation and duplicate symbols are
averaged. Labels come from the series matrix metadata (or metadata.csv).
"""

import os
import re

import pandas as pd

# --- CONFIG
DATA_DIR = "data"
RESULT_DIR = "results"

SERIES_MATRIX_FILES = {
    "GSE209567": os.path.join(DATA_DIR, "GSE209567_series_matrix.txt"),
    "GSE17114": os.path.join(DATA_DIR, "GSE17114_series_matrix.txt"),
}

# Platform annotation used to map probe IDs to gene symbols (GPL570 here).
PLATFORM_ANNOTATION = os.path.join(DATA_DIR, "GPL570-55999.txt")
PROBE_MAP_FILE = os.path.join(RESULT_DIR, "probe_to_gene_symbol.csv")

# Metadata file holding the sample IDs and the case/control column.
METADATA_FILE = os.path.join(DATA_DIR, "metadata.csv")
SAMPLE_ID_COLUMN = "SampleID_CEL"
LABEL_COLUMN = "Condition"
LABEL_ENCODING = {"Patient": 1, "Control": 0}

os.makedirs(RESULT_DIR, exist_ok=True)


# 1. Series matrix -> CSV
def series_matrix_to_csv(txt_path, csv_path):
    """Read a GEO series matrix file, dropping the '!' metadata lines."""
    df = pd.read_csv(txt_path, sep="\t", comment="!", low_memory=False)
    df.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path} with shape {df.shape}")
    return df


# --- 2. Platform annotation -> probe/gene map
def build_probe_map(annotation_path, output_path,
                    probe_col="ID", symbol_col="Gene Symbol"):
    """Keep only the probe ID and gene symbol columns of a GPL annotation."""
    annotation = pd.read_csv(annotation_path, sep="\t", comment="#",
                             low_memory=False)
    probe_map = annotation[[probe_col, symbol_col]].dropna(subset=[symbol_col])
    probe_map.to_csv(output_path, index=False)
    print(f"Wrote {output_path} with {len(probe_map)} probes")
    return probe_map


# ---- 3. Probe-level matrix -> gene-symbol matrix -------
def collapse_to_gene_symbols(expression_file, probe_map_file, output_file,
                             probe_col="ID", symbol_col="Gene Symbol"):
    """Merge with the probe map, drop unmapped probes, average duplicates."""
    expression = pd.read_csv(expression_file)
    probe_map = pd.read_csv(probe_map_file)

    expression = expression.rename(columns={"ID_REF": probe_col})
    merged = pd.merge(expression, probe_map, on=probe_col, how="inner")
    merged = merged.drop(columns=[probe_col])
    merged = merged.dropna(subset=[symbol_col])
    merged = merged.set_index(symbol_col)
    merged = merged.groupby(merged.index).mean()

    merged.to_csv(output_file)
    print(f"Wrote {output_file} with shape {merged.shape}")
    return merged


# 4. Sample labels from the series matrix metadata
def extract_sample_labels(matrix_file, output_file):
    """Pull GSM accessions and the disease-status characteristic from metadata."""
    with open(matrix_file, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    sample_ids = []
    for line in lines:
        if line.startswith("!Sample_geo_accession"):
            sample_ids = line.strip().split("\t")[1:]
            break

    disease_labels = []
    for line in lines:
        if line.startswith("!Sample_characteristics_ch1") and re.search(
            r"disease status|group|condition", line, re.IGNORECASE
        ):
            disease_labels = line.strip().split("\t")[1:]
            break

    # "disease status: Behcet" -> "Behcet"
    simplified = [
        label.split(":")[-1].strip() if ":" in label else label.strip()
        for label in disease_labels
    ]

    labels = pd.DataFrame({"sample": sample_ids, "label": simplified})
    labels.to_csv(output_file, index=False)
    print(f"Wrote {output_file} with {len(labels)} samples")
    return labels


# --- 5. Utility: strip stray quotes left by some exports
def strip_quotes(csv_path, output_path):
    """Some exports leave stray quotes in the text columns. Strip them."""
    df = pd.read_csv(csv_path)
    for col in df.columns:
        if df[col].dtype == "object":
            df[col] = (df[col]
                       .str.replace('"', "", regex=False)
                       .str.replace("'", "", regex=False)
                       .str.strip())
    df.to_csv(output_path, index=False)
    print(f"Cleaned {csv_path} -> {output_path}")
    return df


# 6. Attach labels to the expression matrix
def attach_labels(expression_file, metadata_file, output_file):
    """Transpose to samples x genes and join the condition column."""
    expression = pd.read_csv(expression_file, index_col=0)
    print(f"Expression matrix shape: {expression.shape}")

    if expression.index.duplicated().any():
        expression = expression.groupby(expression.index).mean()
        print(f"Collapsed duplicate symbols -> {expression.shape}")

    metadata = pd.read_csv(metadata_file).set_index(SAMPLE_ID_COLUMN)

    merged = expression.T.join(metadata[[LABEL_COLUMN]])

    missing = merged[merged[LABEL_COLUMN].isnull()].index.tolist()
    if missing:
        print(f"Warning: {len(missing)} samples have no label, e.g. {missing[:10]}")

    merged[LABEL_COLUMN] = merged[LABEL_COLUMN].replace(LABEL_ENCODING)
    merged.to_csv(output_file)
    print(f"Wrote {output_file} with shape {merged.shape}")
    return merged


# --- Pipeline
if __name__ == "__main__":
    # 1. Series matrix -> CSV
    for accession, path in SERIES_MATRIX_FILES.items():
        if os.path.exists(path):
            series_matrix_to_csv(
                path, os.path.join(RESULT_DIR, f"{accession}_series_matrix.csv")
            )

    # 2. Probe -> gene symbol map
    if os.path.exists(PLATFORM_ANNOTATION):
        build_probe_map(PLATFORM_ANNOTATION, PROBE_MAP_FILE)

    # 3. Collapse probes to gene symbols
    probe_level = os.path.join(RESULT_DIR, "GSE209567_series_matrix.csv")
    gene_level = os.path.join(RESULT_DIR, "GSE209567_expression_genesymbol.csv")
    if os.path.exists(probe_level) and os.path.exists(PROBE_MAP_FILE):
        collapse_to_gene_symbols(probe_level, PROBE_MAP_FILE, gene_level)

    # 4. Sample labels
    for accession, path in SERIES_MATRIX_FILES.items():
        if os.path.exists(path):
            extract_sample_labels(
                path, os.path.join(RESULT_DIR, f"{accession}_labels.csv")
            )

    # 5. Labelled matrix (samples x genes + condition)
    if os.path.exists(gene_level) and os.path.exists(METADATA_FILE):
        attach_labels(
            gene_level,
            METADATA_FILE,
            os.path.join(RESULT_DIR, "label_expression_matrix.csv"),
        )
