#!/usr/bin/env python3
"""
Culture medium application to models.
Apply culture medium constraints to metabolic models based on Excel definitions.

Process:
  1. Close all exchange reactions
  2. Apply medium constraints from Excel file
  3. Save modified model

Excel format requirements:
  - Columns: "ID", "LOWER_BOUND", "UPPER_BOUND"
  - Sheet name specifies which medium to apply
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import cobra
from cobra.io import read_sbml_model, write_sbml_model


def apply_medium(model, excel_file, sheet_name, verbose=False):
    """
    Apply culture medium constraints to metabolic model.
    
    Args:
        model: cobra.Model object to modify
        excel_file: Path to Excel file with medium definitions
        sheet_name: Sheet name containing medium composition
        verbose: Print detailed information
        
    Returns:
        Modified cobra.Model with applied medium constraints
        
    Raises:
        FileNotFoundError: If Excel file does not exist
        ValueError: If required columns are missing
    """
    # Close all exchange reactions
    for rxn in model.exchanges:
        rxn.bounds = (0, 0)

    # Read medium composition
    medium_df = pd.read_excel(excel_file, sheet_name=sheet_name)

    # Apply medium constraints
    for _, row in medium_df.iterrows():
        rxn_id = row["ID"]
        
        # Handle numeric formatting (comma decimal separators)
        lower = float(str(row["LOWER_BOUND"]).replace(",", "."))
        upper = float(str(row["UPPER_BOUND"]).replace(",", "."))
        
        # Apply bounds if reaction exists
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.lower_bound = lower
            rxn.upper_bound = upper
        else:
            if verbose:
                print(f"[WARNING] Reaction {rxn_id} not found in model, skipping")

    return model


def set_media(model_path, excel_file, sheet_name, outdir, outname=None, 
              force_overwrite=False, verbose=True):
    """
    Load model, apply medium, and save result.
    
    Args:
        model_path: Path to input model file
        excel_file: Path to Excel file with medium definitions
        sheet_name: Sheet name with medium composition
        outdir: Output directory
        outname: Custom output filename (optional)
        force_overwrite: Overwrite existing files
        verbose: Print progress messages
        
    Returns:
        Path to generated output file
        
    Raises:
        FileNotFoundError: If model or Excel file does not exist
        FileExistsError: If output exists and force_overwrite is False
    """
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")
    
    if not os.path.exists(excel_file):
        raise FileNotFoundError(f"Excel file not found: {excel_file}")

    if verbose:
        print(f"[INFO] Loading model: {model_path}")

    # Load model
    model = read_sbml_model(model_path)

    if verbose:
        print(f"[INFO] Applying medium from sheet: {sheet_name}")

    # Apply medium
    modified_model = apply_medium(model, excel_file, sheet_name, verbose)

    # Generate output path
    os.makedirs(outdir, exist_ok=True)
    model_basename = os.path.splitext(os.path.basename(model_path))[0]
    
    if outname:
        output_filename = outname
    else:
        output_filename = f"{model_basename}_{sheet_name}.xml"
    
    output_path = os.path.join(outdir, output_filename)

    # Check for existing file
    if os.path.exists(output_path) and not force_overwrite:
        raise FileExistsError(
            f"Output already exists: {output_path}. Use --force to overwrite."
        )

    # Save modified model
    write_sbml_model(modified_model, output_path)

    if verbose:
        print(f"[OK] Model with applied medium saved to: {output_path}")

    return output_path


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Apply culture medium constraints to metabolic models",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--config", type=str, default="config.yaml",
                       help="Configuration YAML file")
    parser.add_argument("--model", "--models", nargs="+", dest="models",
                       help="Model file path(s) - supports multiple models (overrides config)")
    parser.add_argument("--excel", type=str,
                       help="Excel file with medium definitions (overrides config)")
    parser.add_argument("--sheet", type=str,
                       help="Excel sheet name with medium composition (overrides config)")
    parser.add_argument("--outdir", type=str,
                       help="Output directory (overrides config)")
    parser.add_argument("--outname", type=str,
                       help="Custom output filename")
    parser.add_argument("--verbose", action="store_true",
                       help="Show detailed progress messages")
    parser.add_argument("--force", action="store_true",
                       help="Overwrite existing output files")
    
    args = parser.parse_args()

    # Load configuration
    config = {}
    if os.path.exists(args.config):
        with open(args.config, "r") as f:
            config = yaml.safe_load(f) or {}
    else:
        print(f"[WARNING] Config file not found, using CLI arguments only")

    # Get media section from config
    prep_config = config.get("preparation", {})
    media_config = prep_config.get("media") or {}

    # Priority: CLI > config
    model_paths = args.models if args.models else prep_config.get("models") or prep_config.get("model_path")
    
    # Convert single model to list for uniform processing
    if isinstance(model_paths, str):
        model_paths = [model_paths]
    
    excel_file = args.excel if args.excel else media_config.get("excel_file")
    sheet_name = args.sheet if args.sheet else media_config.get("sheet_name")
    outdir = args.outdir if args.outdir else media_config.get("output_dir", "outputs/custom_media_models")
    outname = args.outname if args.outname else media_config.get("output_name")
    verbose = args.verbose or media_config.get("verbose", False)
    force_overwrite = args.force or media_config.get("force", False)

    # Validate required parameters
    if not model_paths:
        parser.error("No models specified. Use --model or define in config.")
    if not excel_file:
        parser.error("No Excel file specified. Use --excel or define in config.")
    if not sheet_name:
        parser.error("No sheet name specified. Use --sheet or define in config.")

    # Process each model
    output_paths = []
    for model_path in model_paths:
        try:
            # Generate unique output name for batch processing
            if len(model_paths) > 1 and outname:
                model_basename = os.path.splitext(os.path.basename(model_path))[0]
                current_outname = f"{outname}_{model_basename}.xml"
            else:
                current_outname = outname
            
            output_path = set_media(
                model_path=model_path,
                excel_file=excel_file,
                sheet_name=sheet_name,
                outdir=outdir,
                outname=current_outname,
                force_overwrite=force_overwrite,
                verbose=verbose
            )
            output_paths.append(output_path)
            
        except Exception as e:
            print(f"[ERROR] Failed to process {model_path}: {e}", file=sys.stderr)
            # Continue with other models instead of failing completely
    
    if verbose:
        print(f"\n[INFO] Media application completed for {len(output_paths)}/{len(model_paths)} model(s)")
    
    return 0 if len(output_paths) == len(model_paths) else 1


if __name__ == "__main__":
    exit(main())