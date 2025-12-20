#!/usr/bin/env python3
"""
Config and samplesheet validation tool for Pompe microarray analysis.

This tool validates:
- YAML parameter files against expected schema
- Sample sheet structure and content
- CEL file references and existence
- Configuration consistency

Usage:
    python tools/validate_config.py

Exit codes:
    0: All validations passed
    1: Validation errors found
"""

import sys
import csv
from pathlib import Path
from typing import Dict, List, Any

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not found. Install via: pip install pyyaml", file=sys.stderr)
    sys.exit(1)


class ValidationError(Exception):
    """Custom exception for validation failures."""
    pass


def validate_params_yaml(params_path: Path) -> Dict[str, Any]:
    """
    Validate params.yaml structure and required fields.
    
    Args:
        params_path: Path to params.yaml file
        
    Returns:
        Parsed parameters dictionary
        
    Raises:
        ValidationError: If validation fails
    """
    if not params_path.exists():
        raise ValidationError(f"params.yaml not found at {params_path}")
    
    try:
        with open(params_path, 'r') as f:
            params = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ValidationError(f"Invalid YAML syntax in params.yaml: {e}")
    
    if params is None:
        raise ValidationError("params.yaml is empty")
    
    # Validate required sections
    required_sections = ['deg', 'plots', 'enrichment', 'multimir']
    for section in required_sections:
        if section not in params:
            raise ValidationError(f"Missing required section '{section}' in params.yaml")
    
    # Validate deg parameters
    deg = params['deg']
    if 'logFC' not in deg or not isinstance(deg['logFC'], (int, float)):
        raise ValidationError("deg.logFC must be a numeric value")
    if deg['logFC'] < 0:
        raise ValidationError("deg.logFC must be non-negative")
    
    if 'padj' not in deg or not isinstance(deg['padj'], (int, float)):
        raise ValidationError("deg.padj must be a numeric value")
    if not (0 < deg['padj'] <= 1):
        raise ValidationError("deg.padj must be between 0 and 1")
    
    # Validate plots parameters
    plots = params['plots']
    if 'volcano_topn' not in plots or not isinstance(plots['volcano_topn'], int):
        raise ValidationError("plots.volcano_topn must be an integer")
    if plots['volcano_topn'] < 1:
        raise ValidationError("plots.volcano_topn must be positive")
    
    # Validate multimir parameters
    multimir = params['multimir']
    if 'enabled' not in multimir or not isinstance(multimir['enabled'], bool):
        raise ValidationError("multimir.enabled must be a boolean")
    
    return params


def validate_samplesheet(samplesheet_path: Path, metadata_dir: Path) -> List[Dict[str, str]]:
    """
    Validate samplesheet.csv structure and content.
    
    Args:
        samplesheet_path: Path to samplesheet.csv
        metadata_dir: Path to metadata directory containing CEL files
        
    Returns:
        List of sample records
        
    Raises:
        ValidationError: If validation fails
    """
    if not samplesheet_path.exists():
        raise ValidationError(f"samplesheet.csv not found at {samplesheet_path}")
    
    try:
        with open(samplesheet_path, 'r') as f:
            # Skip empty lines
            lines = [line for line in f if line.strip()]
            reader = csv.DictReader(lines)
            samples = list(reader)
    except Exception as e:
        raise ValidationError(f"Failed to read samplesheet.csv: {e}")
    
    if not samples:
        raise ValidationError("samplesheet.csv is empty")
    
    # Validate required columns
    required_cols = ['sample_id', 'group', 'filename']
    first_row_keys = set(samples[0].keys())
    missing_cols = set(required_cols) - first_row_keys
    if missing_cols:
        raise ValidationError(
            f"Missing required columns in samplesheet.csv: {', '.join(sorted(missing_cols))}"
        )
    
    # Validate sample data
    sample_ids_seen = set()
    filenames_seen = set()
    
    for idx, sample in enumerate(samples, start=2):
        # Check for empty required fields
        for col in required_cols:
            if not sample.get(col, '').strip():
                raise ValidationError(
                    f"Row {idx}: Required column '{col}' is empty"
                )
        
        # Check for duplicate sample IDs
        sample_id = sample['sample_id']
        if sample_id in sample_ids_seen:
            raise ValidationError(
                f"Row {idx}: Duplicate sample_id '{sample_id}'"
            )
        sample_ids_seen.add(sample_id)
        
        # Check for duplicate filenames
        filename = sample['filename']
        if filename in filenames_seen:
            raise ValidationError(
                f"Row {idx}: Duplicate filename '{filename}'"
            )
        filenames_seen.add(filename)
        
        # Validate group values
        group = sample['group']
        if group not in ['Control', 'Pompe']:
            raise ValidationError(
                f"Row {idx}: Invalid group '{group}'. Must be 'Control' or 'Pompe'"
            )
        
        # Check CEL file existence (if metadata_dir exists)
        if metadata_dir.exists():
            cel_path = metadata_dir / filename
            if not cel_path.exists():
                print(
                    f"WARNING: Row {idx}: CEL file not found: {filename}",
                    file=sys.stderr
                )
    
    return samples


def validate_environment_yml(env_path: Path) -> None:
    """
    Validate environment.yml structure.
    
    Args:
        env_path: Path to environment.yml
        
    Raises:
        ValidationError: If validation fails
    """
    if not env_path.exists():
        raise ValidationError(f"environment.yml not found at {env_path}")
    
    try:
        with open(env_path, 'r') as f:
            env = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ValidationError(f"Invalid YAML syntax in environment.yml: {e}")
    
    if env is None:
        raise ValidationError("environment.yml is empty")
    
    # Validate required fields
    if 'name' not in env:
        raise ValidationError("environment.yml missing 'name' field")
    
    if 'dependencies' not in env or not env['dependencies']:
        raise ValidationError("environment.yml missing or empty 'dependencies'")


def main() -> int:
    """
    Main validation entry point.
    
    Returns:
        Exit code (0 for success, 1 for failure)
    """
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent
    
    print("Pompe Microarray Analysis - Configuration Validator")
    print("=" * 60)
    
    validation_passed = True
    
    # Validate params.yaml
    print("\n[1/3] Validating params.yaml...")
    try:
        params_path = repo_root / "config" / "params.yaml"
        params = validate_params_yaml(params_path)
        print("  ✓ params.yaml is valid")
        print(f"    - logFC threshold: {params['deg']['logFC']}")
        print(f"    - padj threshold: {params['deg']['padj']}")
        print(f"    - multiMiR enabled: {params['multimir']['enabled']}")
    except ValidationError as e:
        print(f"  ✗ FAILED: {e}", file=sys.stderr)
        validation_passed = False
    
    # Validate samplesheet.csv
    print("\n[2/3] Validating samplesheet.csv...")
    try:
        samplesheet_path = repo_root / "config" / "samplesheet.csv"
        metadata_dir = repo_root / "metadata"
        samples = validate_samplesheet(samplesheet_path, metadata_dir)
        print("  ✓ samplesheet.csv is valid")
        print(f"    - Total samples: {len(samples)}")
        
        groups = {}
        for sample in samples:
            group = sample['group']
            groups[group] = groups.get(group, 0) + 1
        
        for group, count in sorted(groups.items()):
            print(f"    - {group}: {count} samples")
    except ValidationError as e:
        print(f"  ✗ FAILED: {e}", file=sys.stderr)
        validation_passed = False
    
    # Validate environment.yml
    print("\n[3/3] Validating environment.yml...")
    try:
        env_path = repo_root / "environment" / "environment.yml"
        validate_environment_yml(env_path)
        print("  ✓ environment.yml is valid")
    except ValidationError as e:
        print(f"  ✗ FAILED: {e}", file=sys.stderr)
        validation_passed = False
    
    # Summary
    print("\n" + "=" * 60)
    if validation_passed:
        print("✓ All validations PASSED")
        return 0
    else:
        print("✗ Some validations FAILED", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
