#!/usr/bin/env python3
"""
Functional comparison of GEM models.
Compares predicted FBA fluxes across models, focusing on biomass flux
and central metabolic pathways (EMP, TCA, PPP).
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from cobra.io import read_sbml_model
import scienceplots


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


def safe_get_flux(solution, rxn_id):
    """
    Safely extract flux value from optimization solution.
    
    Args:
        solution: COBRApy solution object
        rxn_id: Reaction identifier
        
    Returns:
        Flux value or 0.0 if not available
    """
    if solution is None:
        return 0.0
    try:
        flux_value = solution.fluxes.get(rxn_id)
        if flux_value is None or np.isnan(flux_value):
            return 0.0
        return float(flux_value)
    except Exception:
        return 0.0

def functional_comparison(models, pathways, outdir, output_name="Functional", 
                         verbose=True, force_overwrite=False):
    """
    Perform functional comparison of multiple GEM models via FBA.
    
    Args:
        models: List of paths to SBML model files
        pathways: Dictionary mapping pathway names to reaction ID lists
        outdir: Output directory for results
        output_name: Base name for output files
        verbose: Print progress messages
        force_overwrite: Overwrite existing output files
        
    Returns:
        Dictionary containing paths to generated files and DataFrame
        
    Raises:
        FileExistsError: If output exists and force_overwrite is False
    """
    os.makedirs(outdir, exist_ok=True)
    records = []

    # Process each model
    for model_path in models:
        model_name = os.path.basename(model_path)
        if verbose:
            print(f"\n[INFO] Loading {model_name}")

        try:
            model = load_model(model_path)
        except Exception as e:
            print(f"[ERROR] Could not load {model_path}: {e}", file=sys.stderr)
            continue

        # Run FBA optimization
        try:
            solution = model.optimize()
        except Exception as e:
            print(f"[ERROR] Optimization failed for {model_name}: {e}", file=sys.stderr)
            continue

        biomass_flux = solution.objective_value if solution else float("nan")

        # Extract pathway fluxes
        for pathway_name, rxn_list in pathways.items():
            for rxn_id in rxn_list:
                flux = safe_get_flux(solution, rxn_id)
                records.append({
                    "model": model_name,
                    "pathway": pathway_name,
                    "reaction": rxn_id,
                    "flux": flux,
                    "biomass_flux": biomass_flux
                })

        if verbose:
            print(f"[OK] Biomass flux: {biomass_flux:.4f}")

    df = pd.DataFrame(records)

    # Save CSV results
    csv_path = os.path.join(outdir, f"{output_name}.csv")
    if os.path.exists(csv_path) and not force_overwrite:
        raise FileExistsError(f"CSV already exists: {csv_path}. Use --force to overwrite.")
    df.to_csv(csv_path, index=False)

    # Create pathway comparison subplots
    subplots_path = None
    try:
        subplots_path = create_pathway_subplots(df, pathways, outdir, output_name, verbose)
    except Exception as e:
        print(f"[WARNING] Could not create pathway subplots: {e}")

    if verbose:
        print(f"\n[DONE] Results saved to:")
        print(f"  - CSV: {csv_path}")
        if subplots_path:
            print(f"  - Pathway subplots: {subplots_path}")

    return {
        "csv": csv_path, 
        "pathway_subplots": subplots_path, 
        "data": df
    }

def create_pathway_subplots(flux_dataframe, pathways, output_directory, output_name, verbose=True):
    """
    Create comparative subplots for pathway flux analysis across multiple models.
    
    Args:
        flux_dataframe: DataFrame containing flux data from functional_comparison
        pathways: Dictionary mapping pathway names to reaction ID lists
        output_directory: Directory to save the generated plot
        output_name: Base name for output file
        verbose: Print progress messages
        
    Returns:
        str: Path to the saved plot file
    """
 
    # Prepare data for plotting
    model_names = flux_dataframe['model'].unique()
    num_models = len(model_names)
    
    if num_models == 0:
        raise ValueError("No model data found for plotting")
    
    # Define color palette for multiple models
    if num_models <= 8:
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#C73E1D', '#3E8914', '#8A1C7C', '#1B5299', '#F0C808']
    else:
        # Generate colors if more than 8 models
        colors = plt.cm.Set3(np.linspace(0, 1, num_models))
    
    plt.style.use(['science','grid' , 'no-latex'])
    # Create figure with subplots
    fig, axes = plt.subplots(1, len(pathways), figsize=(5 * len(pathways), 5))
    if len(pathways) == 1:
        axes = [axes]  # Ensure axes is iterable for single pathway
    
    for idx, (pathway_name, reactions) in enumerate(pathways.items()):
        ax = axes[idx]
        
        # Extract flux values for each reaction in current pathway for all models
        model_fluxes = {model: [] for model in model_names}
        reaction_labels = []
        
        for reaction in reactions:
            for model_name in model_names:
                flux_values = flux_dataframe[
                    (flux_dataframe['model'] == model_name) & 
                    (flux_dataframe['reaction'] == reaction) & 
                    (flux_dataframe['pathway'] == pathway_name)
                ]['flux'].values
                
                # Handle None or NaN values safely
                flux_value = 0.0
                if len(flux_values) > 0:
                    flux_value = flux_values[0] if flux_values[0] is not None else 0.0
                    if np.isnan(flux_value):
                        flux_value = 0.0
                
                model_fluxes[model_name].append(flux_value)
            
            reaction_labels.append(reaction)
        
        # Configure bar positions
        x_positions = np.arange(len(reactions))
        total_width = 0.8  # Total width allocated for all bars
        bar_width = total_width / num_models
        
        # Create bars for each model
        all_bars = []
        for i, model_name in enumerate(model_names):
            bar_offset = (i - (num_models - 1) / 2) * bar_width
            bars = ax.bar(x_positions + bar_offset, model_fluxes[model_name], bar_width,
                         label=model_name, color=colors[i], alpha=0.8,
                         edgecolor='black', linewidth=0.5)
            all_bars.append(bars)
        
        # Customize subplot appearance
        ax.set_title(f'{pathway_name}', fontsize=14, fontweight='bold', pad=10)
        ax.set_xlabel('Reacciones', fontsize=10)
        ax.set_ylabel('Flujo (mmol/gDW/h)', fontsize=10)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(reaction_labels, rotation=45, ha='right', fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add flux values on top of bars - SAFELY handle None/NaN values
        all_flux_values = [flux for fluxes in model_fluxes.values() for flux in fluxes]
        valid_flux_values = [f for f in all_flux_values if f is not None and not np.isnan(f)]
        max_flux = max(valid_flux_values) if valid_flux_values else 1.0
        label_offset = 0.02 * max_flux if max_flux > 0 else 0.02
        
        for model_bars, model_color in zip(all_bars, colors):
            for bar in model_bars:
                height = bar.get_height()
                # Safe formatting - check for None and NaN before formatting
                if height is not None and not np.isnan(height) and abs(height) > 1e-6:
                    try:
                        formatted_value = f'{height:.2f}'
                        ax.text(bar.get_x() + bar.get_width()/2., height + label_offset,
                               formatted_value, ha='center', va='bottom', fontsize=7,
                               fontweight='bold', color=model_color)
                    except (TypeError, ValueError):
                        # Fallback if formatting fails
                        ax.text(bar.get_x() + bar.get_width()/2., height + label_offset,
                               'N/A', ha='center', va='bottom', fontsize=7,
                               fontweight='bold', color=model_color)
    
    # Add unified legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.02),
               ncol=min(4, num_models), frameon=True, fancybox=True, shadow=True, fontsize=10)
    
    plt.tight_layout()
    
    # Save the plot
    plot_filename = f"{output_name}_pathway_subplots.png"
    plot_filepath = os.path.join(output_directory, plot_filename)
    plt.savefig(plot_filepath, dpi=400, bbox_inches='tight')
    plt.close()
    
    if verbose:
        print(f"[INFO] Pathway comparison subplots saved to: {plot_filepath}")
        print(f"[INFO] Compared {num_models} models across {len(pathways)} pathways")
    
    return plot_filepath

def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Functional comparison of GEM models via FBA",
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
    parser.add_argument("--force", action="store_true",
                       help="Overwrite existing output files")
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

    cfg_section = config.get("comparison", {}).get("functional") or {}
    cfg_section2=config.get("comparison") or {}
    # Priority: CLI > config.yaml
    models = args.models if args.models else cfg_section2.get("models", [])
    if not models:
        parser.error("No models specified. Use --models or define 'comparison.functional.models' in config.yaml")

    # Get pathway definitions from config (with defaults)
    pathways = cfg_section.get("pathways", {
        "EMP": ["GLCpts", "PGI", "PFK","FBA", "TPI","GAPD","PGK","PGM","ENO", "PYK"],
        "TCA": ["CS", "ACONT", "ICDHyr","AKGDH","SUCOAS", "SUCD1", "FUM", "MDH","PDH", "PC"],
        "PPP": ["G6PDH2r","PGL","GND","RPI","RPE","TKT1","TALA","TKT2"]
    })

    outdir = args.outdir if args.outdir else cfg_section.get("output_dir", "outputs/comparison/functional")
    output_name = args.output_name if args.output_name else cfg_section.get("output_name", "Functional")
    verbose = args.verbose or cfg_section.get("verbose", False)
    force_overwrite = args.force or cfg_section.get("force", False)

    if verbose:
        print("[INFO] Functional comparison configuration:")
        print(f"  Models: {len(models)}")
        print(f"  Pathways: {list(pathways.keys())}")
        print(f"  Output directory: {outdir}")
        print(f"  Output name: {output_name}")

    try:
        functional_comparison(models, pathways, outdir, output_name, 
                            verbose=verbose, force_overwrite=force_overwrite)
    except Exception as e:
        print(f"[ERROR] {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())