# Reproducibility

1. Create the software environment:
   ```bash
   conda env create -f environment/environment.yml
   conda activate pompe
   ```
2. Pull required metadata tracked with Git LFS (CEL files and `pheno.csv`):
   ```bash
   git lfs install
   git lfs pull
   ```
3. Verify `config/samplesheet.csv` matches your CEL files and groups.
4. Validate configuration (recommended):
   ```bash
   python tools/validate_config.py
   ```
5. Run preflight check (optional but recommended):
   ```bash
   python tools/preflight_check.py
   ```
6. Run the analysis (disable network-heavy multiMiR when offline/CI):
   ```bash
   SKIP_MULTIMIR=1 Rscript microarray_analysis_pompe.R
   ```
7. Outputs (CSV, PDF, session info) are written to `results/` by default. Override with `RESULTS_DIR=/path/to/out`.

Session information is captured in `results/sessionInfo.txt` for transparency.
