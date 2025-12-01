#!/usr/bin/env python3
"""
Structural comparison of GEM models.
Counts genes, reactions, and metabolites in each model and generates
a CSV file and comparative bar plot.
"""

import os
import argparse
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from cobra.io import read_sbml_model


def load_model(model_path):
    """
    Load SBML model using COBRApy.
    
    Args:
        model_path: Path to SBML model file
        
    Returns:
        Loaded cobra.Model object
        
    Raises:
        FileNotFoundError: If model file does not exist
        RuntimeError: If model cannot be loaded
    """
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")
    try:
        return read_sbml_model(model_path)
    except Exception as e:
        raise RuntimeError(f"Error loading model {model_path}: {e}")


def structural_comparison(models, outdir, output_name=None, verbose=True):
    """
    Perform structural comparison of multiple GEM models.
    
    Args:
        models: List of paths to SBML model files
        outdir: Output directory for results
        output_name: Base name for output files (default: "Structural")
        verbose: Print progress messages
        
    Returns:
        DataFrame containing comparison results
        
    Raises:
        ValueError: If no models could be processed successfully
    """
    results = []
    
    # Normalize and create output directory
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    
    if verbose:
        print(f"[INFO] Saving results to: {outdir}")

    # Process each model
    for model_path in models:
        try:
            model_path = os.path.abspath(model_path)
            model = load_model(model_path)
            data = {
                "Model": os.path.basename(model_path),
                "Genes": len(model.genes),
                "Reactions": len(model.reactions),
                "Metabolites": len(model.metabolites)
            }
            results.append(data)
            if verbose:
                print(f"[OK] {data['Model']}: {data['Genes']} genes, "
                      f"{data['Reactions']} reactions, {data['Metabolites']} metabolites")
        except Exception as e:
            print(f"[ERROR] Could not process {model_path}: {e}")
            continue

    if not results:
        raise ValueError("No models could be processed successfully")

    df = pd.DataFrame(results)

    # Determine output file names
    if not output_name:
        output_name = "Structural"
    csv_path = os.path.join(outdir, f"{output_name}.csv")
    plot_path = os.path.join(outdir, f"{output_name}_graph.png")

    # Save CSV results
    df.to_csv(csv_path, index=False)

    # Create bar plot
    df.set_index("Model")[["Genes", "Reactions", "Metabolites"]].plot(
        kind="bar", figsize=(10, 6), rot=45, color=["#3498db", "#e74c3c", "#2ecc71"]
    )
    plt.title("Structural Comparison of Models", fontsize=14, fontweight='bold')
    plt.ylabel("Count", fontsize=12)
    plt.xlabel("Model", fontsize=12)
    plt.legend(title="Component", fontsize=10)
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.close()

    if verbose:
        print(f"[DONE] Results saved to:")
        print(f"  - CSV: {csv_path}")
        print(f"  - Plot: {plot_path}")

    return df


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Structural comparison of GEM models",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--config", type=str, default="config.yaml", 
                       help="Configuration YAML file")
    parser.add_argument("--models", "-m", nargs="+", 
                       help="List of model files to compare (overrides YAML)")
    parser.add_argument("--outdir", type=str, 
                       help="Output directory (overrides YAML)")
    parser.add_argument("--output_name", type=str, 
                       help="Base name for output files")
    parser.add_argument("--verbose", action="store_true", 
                       help="Show detailed progress messages")

    args = parser.parse_args()

    # Load configuration
    config = {}
    if os.path.exists(args.config):
        with open(args.config, "r") as f:
            config = yaml.safe_load(f) or {}
    else:
        print(f"[WARNING] Config file not found, using CLI arguments only")

    comp_config = config.get("comparison", {})
    cfg_section = comp_config.get("structural") or {}
    
    # Priority: CLI > config
    models = args.models if args.models else comp_config.get("models")
    # Priority: CLI > config.yaml)
    outdir = args.outdir if args.outdir else cfg_section.get("output_dir", "outputs/comparison/structural")
    output_name = args.output_name if args.output_name else cfg_section.get("output_name", None)
    verbose = args.verbose or cfg_section.get("verbose", False)

    if not models:
        parser.error("No models specified. Use --models or define 'comparison.models' in config.yaml")

    try:
        structural_comparison(models, outdir, output_name, verbose)
    except Exception as e:
        print(f"[ERROR] {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())