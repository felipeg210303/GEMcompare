#!/usr/bin/env python3
"""
Master preparation script for GEM models.
Part 1 of the pipeline: PREPARATION

Orchestrates four sequential processing steps:
  1. Translation: Convert model between namespace formats (BiGG/SEED/KEGG)
  2. Validation: Validate SBML structure and compliance
  3. Quality Control: Comprehensive model assessment using memote
  4. Media Setup: Apply culture medium constraints to model

Usage modes:
  - Configuration-driven: Uses config.yaml for step execution and parameters
  - Batch processing: Supports multiple models in single execution
  - Flexible skipping: Individual steps can be enabled/disabled
"""

import subprocess
import argparse
import yaml
import os
import sys
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
        description="GEM models preparation pipeline - translation, validation, quality control, and media setup",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single model with verbose output
  python "1. Preparation.py" --models model.xml --verbose
  
  # Batch processing
  python "1. Preparation.py" --models model1.xml model2.xml
  
  # Custom configuration
  python "1. Preparation.py" --config my_config.yaml
        """
    )
    
    parser.add_argument("--models", "-m", nargs="+",
                       help="Model file(s) to process (overrides preparation.models in config)")
    parser.add_argument("--config", type=str, default="config.yaml",
                       help="Configuration YAML file (default: config.yaml)")
    parser.add_argument("--verbose", action="store_true",
                       help="Show detailed progress messages")
    
    args = parser.parse_args()
    
    # Load configuration
    try:
        config = load_config(args.config)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
    
    # Extract preparation configuration
    prep_config = config.get("preparation", {})
    
    # Determine models to process (CLI takes precedence)
    models = args.models if args.models else prep_config.get("models", [])
    if not models:
        print("[ERROR] No models specified.", file=sys.stderr)
        print("  Use --models argument or define 'preparation.models' in config.", file=sys.stderr)
        sys.exit(1)
    
    # Determine base output directory
    base_outdir = prep_config.get("output_dir", "outputs/preparation")
    os.makedirs(base_outdir, exist_ok=True)
    
    verbose = args.verbose or prep_config.get("verbose", False)
    
    # Display configuration
    if verbose:
        print("\n" + "="*60)
        print("PREPARATION PIPELINE - INITIALIZATION")
        print("="*60)
        print(f"Models: {len(models)}")
        for model in models:
            print(f"  - {os.path.basename(model)}")
        print(f"Base output directory: {base_outdir}")
        print("="*60)
    
    # Pipeline step configuration
    scripts_config = {
        "translation": {
            "script": "1_Preparation/1.1_Translation.py",
            "config_section": "translation"
        },
        "validation": {
            "script": "1_Preparation/1.2_Validation.py",
            "config_section": "validation"
        },
        "quality": {
            "script": "1_Preparation/1.3_Quality_control.py",
            "config_section": "quality"
        },
        "media": {
            "script": "1_Preparation/1.4_Set_media.py",
            "config_section": "media"
        }
    }
    
    base_dir = Path(__file__).parent
    exit_code = 0
    executed_steps = 0
    skipped_steps = 0

    
    # Process each model through pipeline
    total_models = len(models)
    successful_models = 0
    
    for model_path in models:
        current_model = model_path
        model_failed = False
        
        if verbose:
            print(f"\n{'='*50}")
            print(f"PROCESSING MODEL: {os.path.basename(model_path)}")
            print("="*50)
        
        # Execute each preparation step
        for step_name, step_info in scripts_config.items():
            script_file = step_info["script"]
            config_section = step_info["config_section"]
            script_path = base_dir / script_file
            
            # Check if step should be skipped
            step_config = prep_config.get(config_section) or {}
            if not step_config.get("run", False):
                if verbose:
                    print(f"[SKIP] {step_name} (run: false in config)")
                skipped_steps += 1
                continue
            
            # Check if script exists
            if not script_path.exists():
                print(f"[WARNING] Script not found: {script_file}, skipping")
                skipped_steps += 1
                continue
            
            if verbose:
                print(f"\n[STEP] {step_name.title()}")
            
            # Prepare step-specific arguments
            step_outdir = step_config.get("output_dir", os.path.join(base_outdir, step_name))
            script_args = ["--model", current_model, "--outdir", step_outdir, "--config", args.config]
            
            if verbose:
                script_args.append("--verbose")
            
            # Step-specific argument handling
            if step_name == "translation":
                target_db = step_config.get("translation_db", "bigg")
                script_args += ["--to", target_db]
                outname = step_config.get("output_name")
                if outname:
                    script_args += ["--outname", outname]
            
            elif step_name == "quality":
                outname = step_config.get("output_name")
                if outname:
                    script_args += ["--output_name", outname]
            
            elif step_name == "media":
                excel_file = step_config.get("excel_file")
                sheet_name = step_config.get("sheet_name")
                outname = step_config.get("output_name")
                if excel_file:
                    script_args += ["--excel", excel_file]
                if sheet_name:
                    script_args += ["--sheet", sheet_name]
                if outname:
                    script_args += ["--outname", outname]
            
            # Execute script
            success = run_subscript(str(script_path), script_args, verbose)
            if not success:
                model_failed = True
                exit_code = 1
                print(f"[WARNING] {step_name} failed for {os.path.basename(model_path)}, continuing anyway .")
                
            
            executed_steps += 1
            
            # Update model path for translation output
            if step_name == "translation" and success:
                model_basename = os.path.splitext(os.path.basename(current_model))[0]
                outname = step_config.get("output_name")
                translated_name = outname if outname else f"{model_basename}_translated.xml"
                current_model = os.path.join(step_outdir, translated_name)
        
        # Track successful models
        if not model_failed:
            successful_models += 1
    
    # Pipeline completion
    print("\n" + "="*60)
    if exit_code == 0:
        print("PREPARATION PIPELINE - EXECUTION COMPLETE")
    else:
        print("PREPARATION PIPELINE - COMPLETED WITH ERRORS")
    
    print(f"Execution summary:")
    print(f"  - Models processed: {successful_models}/{total_models}")
    print(f"  - Steps executed: {executed_steps}")
    print(f"  - Steps skipped: {skipped_steps}")
    print(f"Results directory: {base_outdir}")
    print("="*60)
    
    sys.exit(exit_code)


if __name__ == "__main__":
    main()