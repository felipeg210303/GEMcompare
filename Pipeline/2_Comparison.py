#!/usr/bin/env python3
"""
Master comparison script for GEM models.
Part 2 of the pipeline: COMPARISON

Orchestrates three analysis modules:
  1. Structural: Structural comparison (reactions, metabolites, genes)
  2. Functional: Functional comparison via FBA and flux analysis
  3. CarbonSources: Growth analysis across different carbon sources

Usage modes:
  - Pipeline mode: Automatic execution after Preparation
  - Standalone mode: python "2. Comparison.py" --models model1.xml model2.xml
"""

import argparse
import os
import sys
import yaml
import subprocess
from pathlib import Path


def load_config(config_path):
    """
    Load and validate YAML configuration file.
    
    Args:
        config_path: Path to configuration file
        
    Returns:
        Configuration dictionary
        
    Raises:
        FileNotFoundError: If configuration file does not exist
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, "r") as config_file:
        return yaml.safe_load(config_file) or {}


def detect_models(config=None, cli_models=None):
    """
    Detect models using priority system.

    Priority:
    1. Models produced by Preparation
    2. YAML configuration
    3. Command-line input
    """
    models = []
    source = None

    # Priority 1 — Preparation output directory
    prep_dir = None
    if config:
        prep_dir = config.get("preparation", {}).get("output_models",
                    "outputs/prepared_models")

    if prep_dir and os.path.exists(prep_dir):
        prep_models = list(Path(prep_dir).rglob("*.xml"))
        if prep_models:
            models = [str(m) for m in prep_models]
            source = f"Preparation ({len(models)} models)"

    # Priority 2 — YAML config
    if not models and config:
        yaml_models = config.get("comparison", {}).get("models", [])
        if yaml_models:
            models = yaml_models
            source = f"config.yaml ({len(models)} models)"

    # Priority 3 — CLI
    if not models and cli_models:
        models = [str(m) for m in cli_models]
        source = f"CLI ({len(models)} models)"

    if source:
        print(f"[INFO] Models detected from {source}")

    if not models:
        print("[ERROR] No models found for comparison.", file=sys.stderr)
        sys.exit(1)

    # Validate files exist
    valid = [m for m in models if os.path.exists(m)]
    if not valid:
        print("[ERROR] None of the specified models exist.", file=sys.stderr)
        sys.exit(1)

    return valid


def run_subscript(script_path, args, verbose=True):
    """
    Execute subscript with robust error handling.
    
    Args:
        script_path: Path to script to execute
        args: List of arguments for the script
        verbose: Show detailed output
        
    Returns:
        True if execution was successful, False otherwise
    """
    if not os.path.exists(script_path):
        print(f"[ERROR] Script not found: {script_path}", file=sys.stderr)
        return False
        
    try:
        if verbose:
            print(f"[INFO] Executing: {os.path.basename(script_path)}")
            
        result = subprocess.run(
            ["python", script_path] + args,
            check=True,
            capture_output=True,
            text=True
        )
        
        if verbose:
            print(f"[OK] {os.path.basename(script_path)} completed successfully")
            if result.stdout.strip():
                print(f"  Output: {result.stdout.strip()}")
                
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed: {os.path.basename(script_path)}", file=sys.stderr)
        if e.stderr:
            print(f"  {e.stderr.strip()}", file=sys.stderr)
        return False


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="GEM models comparison pipeline - structural, functional, essentiality, and carbon source analyses",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Automatic execution (uses Preparation outputs)
  python "2. Comparison.py" --verbose
  
  # Specific models
  python "2. Comparison.py" --models model1.xml model2.xml
  
  # Custom configuration
  python "2. Comparison.py" --config my_config.yaml
        """
    )

    parser.add_argument("--models", "-m", nargs="+",
                        help="List of model files to compare (overrides YAML)")
    parser.add_argument("--config", default="config.yaml",
                        help="Configuration YAML file (default: config.yaml)")
    parser.add_argument("--skip", nargs="+", choices=["structural", "functional", "carbon_sources", "essentiality", "biolog"],
                        help="Analyses to skip (e.g., --skip structural functional)")
    parser.add_argument("--outdir", type=str,
                        help="Base output directory (overrides YAML)")
    parser.add_argument("--verbose", action="store_true",
                        help="Show detailed progress messages")

    args = parser.parse_args()

    # Load configuration
    try:
        config = load_config(args.config)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    # Base outdir priority: CLI > YAML > default
    base_outdir = args.outdir or config.get("comparison", {}).get("output_dir", "outputs/comparison")
    os.makedirs(base_outdir, exist_ok=True)

    # Detect models
    models = detect_models(config=config, cli_models=args.models)

    verbose = args.verbose or config.get("comparison", {}).get("verbose", False)

    if verbose:
        print("\n" + "="*60)
        print("COMPARISON PIPELINE - INITIALIZATION")
        print("="*60)
        print(f"Models detected: {len(models)}")
        for m in models:
            print(f"  - {os.path.basename(m)}")
        print(f"Base output directory: {base_outdir}")
        print(f"Analyses skipped: {args.skip or 'none'}")
        print("="*60)

    # Analysis script configuration
    scripts_config = {
        "structural": {
            "script": "2_Comparison/2.1_Structural.py",
            "config_section": "structural"
        },
        "functional": {
            "script": "2_Comparison/2.2_Functional.py",
            "config_section": "functional"
        },
        "carbon_sources": {
            "script": "2_Comparison/2.3_Diff_CarbonSources.py",
            "config_section": "carbon_sources"
        },
        "essentiality": {
            "script": "2_Comparison/2.4_GeneEssentiality.py",
            "config_section": "gene_essentiality"
        },
        "biolog": {
            "script": "2_Comparison/2.5_CarbonUtilization.py",
            "config_section": "biolog"
        }
    }

    base_dir = Path(__file__).parent
    exit_code = 0
    executed_scripts = 0
    skipped_scripts = 0

    # === MAIN PIPELINE LOOP ===
    for analysis_name, info in scripts_config.items():
        script_rel = info["script"]
        script_path = base_dir / script_rel
        config_section = info["config_section"]

        # Verbose header
        if verbose:
            print("\n" + "-"*60)
            print(f"[ANALYSIS] {analysis_name.upper()}")
            print("-"*60)

        # Skip via CLI
        if args.skip and analysis_name in args.skip:
            print(f"[SKIP] {analysis_name} (--skip flag)")
            skipped_scripts += 1
            continue

        # Skip via config
        section_conf = config.get("comparison", {}).get(config_section, {})
        if section_conf.get("run") is False:
            print(f"[SKIP] {analysis_name} (run: false in config)")
            skipped_scripts += 1
            continue

        # Script missing
        if not script_path.exists():
            print(f"[WARNING] Script not found: {script_rel}, skipping")
            skipped_scripts += 1
            continue

        # Prepare args
        analysis_outdir = os.path.join(base_outdir, analysis_name)
        script_args = ["--models"] + models + ["--outdir", analysis_outdir]

        # pass config file
        if args.config:
            script_args += ["--config", args.config]

        if verbose:
            script_args.append("--verbose")
            print(f"[INFO] Running script: {script_rel}")
            print(f"[INFO] Output directory: {analysis_outdir}")

        # Run script
        success = run_subscript(str(script_path), script_args, verbose)

        if success:
            executed_scripts += 1
        else:
            exit_code = 1
            print(f"[ERROR] {analysis_name} failed — continuing with remaining analyses.")

    # === FINAL SUMMARY ===
    print("\n" + "="*60)
    if exit_code == 0:
        print("COMPARISON PIPELINE - COMPLETED SUCCESSFULLY")
    else:
        print("COMPARISON PIPELINE - COMPLETED WITH ERRORS")
    print(f"Analyses executed: {executed_scripts}")
    print(f"Analyses skipped: {skipped_scripts}")
    print(f"Results directory: {base_outdir}")
    print("="*60)

    sys.exit(exit_code)


if __name__ == "__main__":
    main()