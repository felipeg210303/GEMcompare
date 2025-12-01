#!/usr/bin/env python3
"""
Model quality control.
Comprehensive quality assessment for GEM models using memote validation suite.

Output includes:
  - Model statistics (reactions, metabolites, genes)
  - Basic quality assessment report
  - Instructions for detailed HTML report generation
"""

import os
import argparse
import sys
import yaml
import cobra
import subprocess
from memote.suite.api import test_model

def quality_assessment(model_path, outdir="outputs/quality_reports",
                      outname=None, force_overwrite=False, verbose=True):

    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")

    if verbose:
        print(f"[INFO] Validating model: {model_path}")

    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Determine output file path
    if outname:
        report_path = os.path.join(outdir, outname)
    else:
        base = os.path.splitext(os.path.basename(model_path))[0]
        report_path = os.path.join(outdir, f"{base}_memote_report.html")

    if os.path.exists(report_path) and not force_overwrite:
        raise FileExistsError(
            f"Report already exists: {report_path}. Use force_overwrite=True."
        )

    if verbose:
        print("[INFO] Running Memote snapshot report (CLI)...")

    # Run Memote CLI snapshot report
    cmd = [
        "memote", "report", "snapshot",
        model_path,
        "--filename", report_path
    ]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Memote failed: {e}")

    if verbose:
        print(f"[DONE] Memote HTML report saved to: {report_path}")

    return report_path

def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Quality control for GEM models using memote",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--config", type=str, default="config.yaml",
                       help="Configuration YAML file")
    parser.add_argument("--model", "--models", nargs="+", dest="models",
                       help="Model file path(s) - supports multiple models (overrides config)")
    parser.add_argument("--outdir", type=str,
                       help="Output directory for quality reports (overrides config)")
    parser.add_argument("--outname", type=str,
                       help="Custom output filename")
    parser.add_argument("--verbose", action="store_true",
                       help="Show detailed progress messages")
    parser.add_argument("--force", action="store_true",
                       help="Overwrite existing quality reports")
    
    args = parser.parse_args()

    # Load configuration
    config = {}
    if os.path.exists(args.config):
        with open(args.config, "r") as f:
            config = yaml.safe_load(f) or {}
    else:
        print(f"[WARNING] Config file not found, using CLI arguments only")

    # Get quality section from config
    prep_config = config.get("preparation", {})
    quality_config = prep_config.get("quality") or {}

    # Priority: CLI > config
    model_paths = args.models if args.models else prep_config.get("models") or prep_config.get("model_path")
    
    # Convert single model to list for uniform processing
    if isinstance(model_paths, str):
        model_paths = [model_paths]
    
    outdir = args.outdir if args.outdir else quality_config.get("output_dir", "outputs/quality_reports")
    outname = args.outname if args.outname else quality_config.get("output_name")
    verbose = args.verbose or quality_config.get("verbose", False)
    force_overwrite = args.force or quality_config.get("force", False)

    if not model_paths:
        parser.error("No models specified. Use --model or define in config.")

    # Process each model
    output_paths = []
    for model_path in model_paths:
        try:
            # Generate unique output name for batch processing
            if len(model_paths) > 1 and outname:
                model_basename = os.path.splitext(os.path.basename(model_path))[0]
                current_outname = f"{outname}_{model_basename}.txt"
            else:
                current_outname = outname
            
            report_path = quality_assessment(
                model_path=model_path,
                outdir=outdir,
                outname=current_outname,
                force_overwrite=force_overwrite,
                verbose=verbose
            )
            output_paths.append(report_path)
            
        except Exception as e:
            print(f"[ERROR] Failed to process {model_path}: {e}", file=sys.stderr)
            # Continue with other models instead of failing completely
    
    if verbose:
        print(f"\n[INFO] Quality assessment completed for {len(output_paths)}/{len(model_paths)} model(s)")
    
    return 0 if len(output_paths) == len(model_paths) else 1


if __name__ == "__main__":
    exit(main())