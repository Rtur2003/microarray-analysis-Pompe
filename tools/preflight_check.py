#!/usr/bin/env python3
"""
Pre-flight validation tool for Pompe microarray analysis pipeline.

This tool performs all necessary checks before running the R analysis:
- Environment readiness
- Configuration validity
- Data file availability
- Dependencies verification

Usage:
    python tools/preflight_check.py

Exit codes:
    0: All checks passed, ready to run analysis
    1: Critical errors found, cannot proceed
    2: Warnings present, can proceed with caution
"""

import sys
import os
import subprocess
from pathlib import Path
from typing import List, Tuple

# Import validation functions from validate_config
sys.path.insert(0, str(Path(__file__).parent))
try:
    from validate_config import (
        validate_params_yaml,
        validate_samplesheet,
        validate_environment_yml,
        ValidationError
    )
except ImportError:
    print("ERROR: Cannot import validate_config module", file=sys.stderr)
    sys.exit(1)


class PreflightCheck:
    """Manages pre-flight validation checks."""
    
    def __init__(self, repo_root: Path):
        self.repo_root = repo_root
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.info: List[str] = []
    
    def add_error(self, message: str) -> None:
        """Add a critical error."""
        self.errors.append(message)
    
    def add_warning(self, message: str) -> None:
        """Add a warning."""
        self.warnings.append(message)
    
    def add_info(self, message: str) -> None:
        """Add informational message."""
        self.info.append(message)
    
    def check_r_installation(self) -> bool:
        """Check if R is installed and accessible."""
        print("[1/6] Checking R installation...")
        try:
            result = subprocess.run(
                ['Rscript', '--version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                version_line = result.stderr.split('\n')[0] if result.stderr else "R version unknown"
                print(f"  ✓ R is installed: {version_line}")
                self.add_info(f"R: {version_line}")
                return True
            else:
                self.add_error("R is installed but not responding correctly")
                print("  ✗ R installation check failed", file=sys.stderr)
                return False
        except FileNotFoundError:
            self.add_error("Rscript not found in PATH. Install R or activate conda environment.")
            print("  ✗ Rscript not found", file=sys.stderr)
            return False
        except Exception as e:
            self.add_error(f"Failed to check R: {e}")
            print(f"  ✗ Error checking R: {e}", file=sys.stderr)
            return False
    
    def check_configuration(self) -> bool:
        """Validate all configuration files."""
        print("\n[2/6] Validating configuration files...")
        config_valid = True
        
        # Validate params.yaml
        try:
            params_path = self.repo_root / "config" / "params.yaml"
            params = validate_params_yaml(params_path)
            print("  ✓ params.yaml is valid")
            self.add_info(f"logFC={params['deg']['logFC']}, padj={params['deg']['padj']}")
        except ValidationError as e:
            self.add_error(f"params.yaml validation failed: {e}")
            print(f"  ✗ params.yaml: {e}", file=sys.stderr)
            config_valid = False
        
        # Validate samplesheet.csv
        try:
            samplesheet_path = self.repo_root / "config" / "samplesheet.csv"
            metadata_dir = self.repo_root / "metadata"
            samples = validate_samplesheet(samplesheet_path, metadata_dir)
            print("  ✓ samplesheet.csv is valid")
            self.add_info(f"Total samples: {len(samples)}")
        except ValidationError as e:
            self.add_error(f"samplesheet.csv validation failed: {e}")
            print(f"  ✗ samplesheet.csv: {e}", file=sys.stderr)
            config_valid = False
        
        # Validate environment.yml
        try:
            env_path = self.repo_root / "environment" / "environment.yml"
            validate_environment_yml(env_path)
            print("  ✓ environment.yml is valid")
        except ValidationError as e:
            self.add_error(f"environment.yml validation failed: {e}")
            print(f"  ✗ environment.yml: {e}", file=sys.stderr)
            config_valid = False
        
        return config_valid
    
    def check_cel_files(self) -> bool:
        """Check availability of CEL files."""
        print("\n[3/6] Checking CEL files...")
        
        samplesheet_path = self.repo_root / "config" / "samplesheet.csv"
        metadata_dir = self.repo_root / "metadata"
        
        if not samplesheet_path.exists():
            self.add_error("samplesheet.csv not found")
            return False
        
        try:
            import csv
            with open(samplesheet_path, 'r') as f:
                lines = [line for line in f if line.strip()]
                reader = csv.DictReader(lines)
                samples = list(reader)
            
            missing_files = []
            found_files = []
            
            for sample in samples:
                filename = sample.get('filename', '')
                if filename:
                    cel_path = metadata_dir / filename
                    if cel_path.exists():
                        found_files.append(filename)
                    else:
                        missing_files.append(filename)
            
            if missing_files:
                self.add_warning(f"{len(missing_files)} CEL files not found in metadata/")
                print(f"  ⚠ {len(missing_files)} CEL files missing (Git LFS may be needed)")
                if found_files:
                    print(f"  ✓ {len(found_files)} CEL files found")
                return len(found_files) > 0
            else:
                print(f"  ✓ All {len(found_files)} CEL files found")
                return True
                
        except Exception as e:
            self.add_error(f"Failed to check CEL files: {e}")
            print(f"  ✗ Error: {e}", file=sys.stderr)
            return False
    
    def check_output_directory(self) -> bool:
        """Check or create output directory."""
        print("\n[4/6] Checking output directory...")
        
        results_dir = os.environ.get('RESULTS_DIR', 'results')
        results_path = self.repo_root / results_dir
        
        try:
            if results_path.exists():
                if results_path.is_dir():
                    print(f"  ✓ Output directory exists: {results_path}")
                    self.add_info(f"Results dir: {results_path}")
                    return True
                else:
                    self.add_error(f"Results path exists but is not a directory: {results_path}")
                    print(f"  ✗ Not a directory: {results_path}", file=sys.stderr)
                    return False
            else:
                results_path.mkdir(parents=True, exist_ok=True)
                print(f"  ✓ Created output directory: {results_path}")
                self.add_info(f"Results dir created: {results_path}")
                return True
        except Exception as e:
            self.add_error(f"Cannot create output directory: {e}")
            print(f"  ✗ Error: {e}", file=sys.stderr)
            return False
    
    def check_environment_variables(self) -> bool:
        """Check relevant environment variables."""
        print("\n[5/6] Checking environment variables...")
        
        skip_multimir = os.environ.get('SKIP_MULTIMIR', '0')
        results_dir = os.environ.get('RESULTS_DIR', 'results')
        
        print(f"  ℹ SKIP_MULTIMIR={skip_multimir}")
        print(f"  ℹ RESULTS_DIR={results_dir}")
        
        if skip_multimir == '1':
            self.add_info("multiMiR queries will be skipped")
        
        return True
    
    def check_pipeline_script(self) -> bool:
        """Verify main pipeline script exists."""
        print("\n[6/6] Checking pipeline script...")
        
        script_path = self.repo_root / "microarray_analysis_pompe.R"
        
        if script_path.exists():
            print(f"  ✓ Pipeline script found: {script_path.name}")
            return True
        else:
            self.add_error("Pipeline script not found: microarray_analysis_pompe.R")
            print("  ✗ Pipeline script not found", file=sys.stderr)
            return False
    
    def print_summary(self) -> None:
        """Print check summary."""
        print("\n" + "=" * 60)
        print("PREFLIGHT CHECK SUMMARY")
        print("=" * 60)
        
        if self.errors:
            print("\n✗ CRITICAL ERRORS:")
            for error in self.errors:
                print(f"  - {error}")
        
        if self.warnings:
            print("\n⚠ WARNINGS:")
            for warning in self.warnings:
                print(f"  - {warning}")
        
        if self.info:
            print("\nℹ INFO:")
            for info in self.info:
                print(f"  - {info}")
        
        print("\n" + "=" * 60)
        if self.errors:
            print("✗ FAILED: Cannot proceed with analysis")
            print("Fix critical errors before running the pipeline.")
            return 1
        elif self.warnings:
            print("⚠ WARNINGS PRESENT: Proceed with caution")
            print("You can run the pipeline, but be aware of warnings.")
            return 2
        else:
            print("✓ ALL CHECKS PASSED: Ready to run analysis")
            print("\nTo run the pipeline:")
            print("  Rscript microarray_analysis_pompe.R")
            return 0


def main() -> int:
    """
    Main preflight check entry point.
    
    Returns:
        Exit code (0=pass, 1=errors, 2=warnings)
    """
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent
    
    print("Pompe Microarray Analysis - Preflight Check")
    print("=" * 60)
    
    checker = PreflightCheck(repo_root)
    
    # Run all checks
    checker.check_r_installation()
    checker.check_configuration()
    checker.check_cel_files()
    checker.check_output_directory()
    checker.check_environment_variables()
    checker.check_pipeline_script()
    
    # Print summary and return exit code
    return checker.print_summary()


if __name__ == '__main__':
    sys.exit(main())
