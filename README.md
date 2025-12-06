# Pompe Transcriptomics (GSE38680)

Reproducible R workflow for analysing the Infantile-Onset Pompe Disease (IOPD) dataset **GSE38680** (Affymetrix U133 Plus 2.0): normalisation, differential expression, functional enrichment, QC plots, and optional multiMiR queries.

## Quick start

```bash
# 1) set up environment
conda env create -f environment/environment.yml
conda activate pompe

# 2) obtain data (LFS)
git lfs install
git lfs pull  # pulls CEL files and pheno.csv into metadata/

# 3) run analysis (multiMiR disabled in CI/offline)
SKIP_MULTIMIR=1 Rscript microarray_analysis_pompe.R
```

Results and plots are written to `results/` (override with `RESULTS_DIR=/path/to/out`).

## Project structure
- `microarray_analysis_pompe.R` – main pipeline script
- `config/` – `samplesheet.csv` and `params.yaml` (DEG thresholds, plot settings, multiMiR flag)
- `metadata/` – CEL files and `pheno.csv` (tracked via Git LFS), provenance and variable dictionary
- `environment/` – Conda environment specification
- `docs/` – additional documentation (methods, reproducibility notes)

## Data and ethics
All inputs come from public GEO accession **GSE38680**. No personal health information is included.

## Citation
Please cite:

> Tosun F., Bayrak H. 2025. *Transcriptomic and Functional Analysis of Pompe Disease*.

See `CITATION.cff` for machine-readable metadata.
