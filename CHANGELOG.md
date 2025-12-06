# Changelog
All notable changes to this project will be documented in this file.

## [0.3.0] - 2025-12-06
### Changed
- Rebuilt R pipeline to consume `config/samplesheet.csv`, respect `params.yaml`, and write to `results/`.
- Added offline/CI guardrails (optional multiMiR via `SKIP_MULTIMIR`, results dir override, session info output).
- Hardened error handling for missing data and annotation.
### Added
- Cleaned metadata docs and corrected sample sheet/variable dictionary.
- Updated Conda environment and CI workflow to create the environment and skip network-heavy steps when data are absent.

## [0.2.0] - 2025-05-04
### Added
- Governance documents: CODE_OF_CONDUCT, CONTRIBUTING, SECURITY, LICENSE, CITATION.cff
- Environment definition and configuration templates
- Basic `.pre-commit-config.yaml` and `.editorconfig`
### Fixed
- Replaced incorrect `.gitignore` content with proper ignore rules

## [0.1.0] - 2024-12-01
### Added
- Initial R script for Pompe microarray analysis
- LFS-tracked CEL files and phenotype metadata
