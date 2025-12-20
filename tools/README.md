# Tools

Python-based validation and automation tools for the Pompe microarray analysis pipeline.

## Overview

This directory contains Python utilities that follow the repository's Python-first directive for tooling, automation, and validation tasks.

## Available Tools

### validate_config.py

Validates configuration files without executing the analysis pipeline.

**Purpose**: Static validation of YAML and CSV configuration files
**Usage**:
```bash
python tools/validate_config.py
```

**What it checks**:
- `config/params.yaml`: Schema, data types, value ranges
- `config/samplesheet.csv`: Required columns, data integrity, duplicates
- `environment/environment.yml`: Structure and required fields

**Exit codes**:
- `0`: All validations passed
- `1`: Validation errors found

**Example output**:
```
Pompe Microarray Analysis - Configuration Validator
============================================================

[1/3] Validating params.yaml...
  ✓ params.yaml is valid
    - logFC threshold: 1.0
    - padj threshold: 0.05
    - multiMiR enabled: True

[2/3] Validating samplesheet.csv...
  ✓ samplesheet.csv is valid
    - Total samples: 19
    - Control: 10 samples
    - Pompe: 9 samples

[3/3] Validating environment.yml...
  ✓ environment.yml is valid

============================================================
✓ All validations PASSED
```

### preflight_check.py

Comprehensive pre-analysis validation tool that checks system readiness.

**Purpose**: Verify all prerequisites before running the R pipeline
**Usage**:
```bash
python tools/preflight_check.py
```

**What it checks**:
1. R installation and version
2. Configuration file validity (via `validate_config`)
3. CEL file availability
4. Output directory creation
5. Environment variables
6. Pipeline script existence

**Exit codes**:
- `0`: All checks passed, ready to run analysis
- `1`: Critical errors, cannot proceed
- `2`: Warnings present, can proceed with caution

**Example output**:
```
Pompe Microarray Analysis - Preflight Check
============================================================
[1/6] Checking R installation...
  ✓ R is installed: R scripting front-end version 4.3.0

[2/6] Validating configuration files...
  ✓ params.yaml is valid
  ✓ samplesheet.csv is valid
  ✓ environment.yml is valid

[3/6] Checking CEL files...
  ✓ All 19 CEL files found

[4/6] Checking output directory...
  ✓ Output directory exists: results

[5/6] Checking environment variables...
  ℹ SKIP_MULTIMIR=1
  ℹ RESULTS_DIR=results

[6/6] Checking pipeline script...
  ✓ Pipeline script found: microarray_analysis_pompe.R

============================================================
✓ ALL CHECKS PASSED: Ready to run analysis

To run the pipeline:
  Rscript microarray_analysis_pompe.R
```

## Dependencies

All tools require only Python 3.6+ and PyYAML, which is already included in the conda environment (`environment/environment.yml`).

## Integration

### Pre-commit Hooks

These tools can be integrated into `.pre-commit-config.yaml`:

```yaml
- repo: local
  hooks:
    - id: validate-config
      name: Validate configuration files
      entry: python tools/validate_config.py
      language: system
      pass_filenames: false
      files: ^config/.*\.(yaml|csv)$
```

### CI/CD Workflow

Add to `.github/workflows/pipelineA.yml` before running the R script:

```yaml
- name: Validate configuration
  shell: bash -l {0}
  run: |
    python tools/validate_config.py

- name: Preflight check
  shell: bash -l {0}
  run: |
    python tools/preflight_check.py
```

### Manual Workflow

Recommended workflow for local development:

```bash
# 1. Validate configuration
python tools/validate_config.py

# 2. Run comprehensive preflight check
python tools/preflight_check.py

# 3. If all checks pass, run the analysis
SKIP_MULTIMIR=1 Rscript microarray_analysis_pompe.R
```

## Design Principles

These tools follow the repository's engineering directives:

1. **Python-first**: All tooling is Python-based as per directive I.1
2. **Static analysis**: No code execution, only validation (directive I.4)
3. **Safety by construction**: Defensive validation, clear error messages
4. **Single responsibility**: Each tool has one clear purpose
5. **Zero side effects**: Tools are read-only except for creating output directories

## Troubleshooting

### "PyYAML not found"

Install PyYAML (already in conda environment):
```bash
conda activate pompe  # or pip install pyyaml
```

### "Rscript not found"

Activate the conda environment:
```bash
conda env create -f environment/environment.yml
conda activate pompe
```

### "CEL files not found"

Pull data files using Git LFS:
```bash
git lfs install
git lfs pull
```

## Future Enhancements

Potential additions following the same principles:

- Schema validation for results files
- Automated data quality reports
- Configuration migration tools
- Result comparison utilities
- Automated documentation generation

All future tools should follow Python-first directive and maintain atomicity.
