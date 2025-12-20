# Contributing

Thanks for your interest in improving the Pompe microarray analysis pipeline!

## Getting Started
- Fork the repository and create your branch from `main`.
- Ensure Git LFS is installed (`git lfs install`).
- Set up the conda environment: `conda env create -f environment/environment.yml && conda activate pompe`.
- Install pre-commit hooks: `pre-commit install` (optional but recommended).
- Validate configuration before making changes: `python tools/validate_config.py`.
- Run the analysis locally to confirm changes: `Rscript microarray_analysis_pompe.R`.

## Coding Standards
- Use R style compatible with `styler` and `lintr`.
- Place data under `metadata/` and results under `results/`.
- Validate configuration files using `python tools/validate_config.py` before committing.
- Run preflight checks before testing: `python tools/preflight_check.py`.
- Update documentation and examples when changing interfaces.

## Pull Requests
- Reference related issues in commit messages.
- Include a summary of changes and testing steps.
- Ensure validation tools pass: `python tools/validate_config.py`.
- Ensure CI workflows pass.

## Reporting Issues
+Please open an issue using the templates in `.github/ISSUE_TEMPLATE` and provide
+as much detail as possible.
 
