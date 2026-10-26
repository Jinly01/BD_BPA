# BD-BPA biomarker pipeline

Code for identifying and validating core genes linking bisphenol A exposure to
Behçet disease, using three public GEO microarray cohorts.

| Accession | Platform | Role |
|---|---|---|
| GSE17114 | GPL570 | training / internal validation |
| GSE209567 | GPL570 | external validation |
| GSE70403 | GPL4133, two-colour | candidate discovery and an independent check |

## Scripts

Run in order from the repository root (paths are relative to it).

```
01_prepare_expression_matrix.py   series matrix -> gene symbol matrix, labels
02_preprocess_gse70403.R          Agilent raw scans -> M-value matrix
03_deg_wgcna.R                    limma DEG + WGCNA, intersect with BPA targets
04_ppi_enrichment.R               KEGG / GO enrichment, STRING PPI
05_feature_selection_loocv.py     LOOCV stability selection -> core genes
06_model_validation.py            logistic regression, ROC / PR / DCA
07_model_interpretation.py        SHAP, nomogram
08_expression_validation.py       heatmaps, boxplots, PCA, radar
09_gse70403_validation.py         one-sample tests on the M-value cohort
10_ssgsea_immune.R                ssGSEA immune infiltration
```

Each script starts with a configuration block holding file paths, cutoffs and
model parameters. Change those rather than the analysis code.

## Setup

```bash
pip install -r requirements.txt
Rscript install_r_packages.R
```

Then put the input files in `data/` — see `data/README.md` for what is needed
and where to get it. Neither `data/` nor `results/` is tracked in git.

## Notes

GSE70403 has no control arrays. Every array is a patient sample hybridised
against a common reference, so the values are M-values and the limma model in
step 03 is a single-group intercept model testing deviation from zero, not a
case-control contrast. Step 09 validates the core genes the same way, with
one-sample t-tests.

The scaling used for `test_scaled.csv` has to come from the training cohort,
not from the test cohort itself.

Seeds are fixed (2025 in Python, 123 in R), but results still move a little
with package versions, WGCNA and GSVA in particular.

## License

MIT, see `LICENSE`.
