#!/usr/bin/env python3
"""
SBML model validation.
Validates SBML models using libSBML library for structure and compliance checking.

Validation checks:
  - Internal consistency: Model structure and mathematical consistency
  - SBML L3v2 compatibility: Compliance with SBML Level 3 Version 2 specification
"""

import os
import argparse
import sys
import yaml
import libsbml


def run_validation_check(sbml_document, check_function, check_name, verbose=True):
    """
    Execute specific libSBML validation check and collect results.
    
    Args:
        sbml_document: libSBML document object
        check_function: Name of libSBML check function to call
        check_name: Human-readable check name
        verbose: Print progress messages
        
    Returns:
        Dictionary containing check name and any errors found
    """
    results = {"check": check_name, "errors": []}
    
    if verbose:
        print(f"[INFO] Running {check_name}")
    
    # Execute validation check
    error_count = getattr(sbml_document, check_function)()
    
    if error_count == 0:
        if verbose:
            print(f"[OK] {check_name}: No errors")
    else:
        if verbose:
            print(f"[WARNING] {check_name}: Found {error_count} error(s)")
        
        # Collect error details
        for i in range(error_count):
            error = sbml_document.getError(i)
            error_msg = f"[{i + 1}] Level {error.getSeverity()} - {error.getMessage()}"
            
            if verbose:
                print(f"  {error_msg}")
            results["errors"].append(error_msg)
    
    return results


def validate_model(model_path, outdir="outputs/validation_logs", 
                   force_overwrite=False, verbose=True):
    """
    Validate SBML model and generate report.
    
    Args:
        model_path: Path to SBML model file
        outdir: Output directory for validation logs
        force_overwrite: Overwrite existing log files
        verbose: Print detailed progress
        
    Returns:
        Path to generated validation log file
        
    Raises:
        FileNotFoundError: If model file does not exist
        FileExistsError: If output exists and force_overwrite is False
    """
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")

    if verbose:
        print(f"[INFO] Loading model: {model_path}")
    
    # Read SBML document
    reader = libsbml.SBMLReader()
    document = reader.readSBML(model_path)

    # Check for initial read errors
    if document.getNumErrors() > 0:
        print("[ERROR] Critical errors reading SBML document:", file=sys.stderr)
        document.printErrors()
        sys.exit(1)
    else:
        if verbose:
            print("[OK] SBML document loaded successfully")

    # Define validation checks
    checks = [
        ("checkInternalConsistency", "Internal consistency"),
        ("checkL3v2Compatibility", "SBML L3v2 compatibility"),
    ]

    # Execute validation checks
    results = []
    for check_func, check_desc in checks:
        result = run_validation_check(document, check_func, check_desc, verbose=verbose)
        results.append(result)

    # Prepare output path
    os.makedirs(outdir, exist_ok=True)
    model_basename = os.path.splitext(os.path.basename(model_path))[0]
    log_path = os.path.join(outdir, f"{model_basename}_validation.txt")

    # Check for existing file
    if os.path.exists(log_path) and not force_overwrite:
        raise FileExistsError(
            f"Validation log already exists: {log_path}. Use --force to overwrite."
        )

    # Write validation report
    with open(log_path, "w") as log_file:
        log_file.write(f"Validation report for: {model_path}\n")
        log_file.write("=" * 60 + "\n\n")
        
        for result in results:
            status = "PASS" if not result["errors"] else f"FAIL ({len(result['errors'])} errors)"
            log_file.write(f"Check: {result['check']} - {status}\n")
            
            for error in result["errors"]:
                log_file.write(f"  {error}\n")
            log_file.write("\n")

    if verbose:
        print(f"[DONE] Validation report saved to: {log_path}")
    
    return log_path


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Validate SBML models using libSBML",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--config", type=str, default="config.yaml",
                       help="Configuration YAML file")
    parser.add_argument("--model", "--models", nargs="+", dest="models",
                       help="Model file path(s) - supports multiple models (overrides config)")
    parser.add_argument("--outdir", type=str,
                       help="Output directory for validation logs (overrides config)")
    parser.add_argument("--verbose", action="store_true",
                       help="Show detailed progress messages")
    parser.add_argument("--force", action="store_true",
                       help="Overwrite existing validation logs")
    
    args = parser.parse_args()

    # Load configuration
    config = {}
    if os.path.exists(args.config):
        with open(args.config, "r") as f:
            config = yaml.safe_load(f) or {}
    else:
        print(f"[WARNING] Config file not found, using CLI arguments only")

    # Get validation section from config
    prep_config = config.get("preparation", {})
    val_config = prep_config.get("validation") or {}

    # Priority: CLI > config
    model_paths = args.models if args.models else prep_config.get("models") or prep_config.get("model_path")
    
    # Convert single model to list for uniform processing
    if isinstance(model_paths, str):
        model_paths = [model_paths]
    
    outdir = args.outdir if args.outdir else val_config.get("output_dir", "outputs/validation_logs")
    verbose = args.verbose or val_config.get("verbose", False)
    force_overwrite = args.force or val_config.get("force", False)

    if not model_paths:
        parser.error("No models specified. Use --model or define in config.")

    # Process each model
    output_paths = []
    for model_path in model_paths:
        try:
            log_path = validate_model(
                model_path=model_path,
                outdir=outdir,
                force_overwrite=force_overwrite,
                verbose=verbose
            )
            output_paths.append(log_path)
            
        except Exception as e:
            print(f"[ERROR] Failed to process {model_path}: {e}", file=sys.stderr)
            # Continue with other models instead of failing completely
    
    if verbose:
        print(f"\n[INFO] Validation completed for {len(output_paths)}/{len(model_paths)} model(s)")
    
    return 0 if len(output_paths) == len(model_paths) else 1


if __name__ == "__main__":
    exit(main())