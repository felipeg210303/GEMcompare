#!/usr/bin/env python3 
"""
Gene essentiality comparison for GEM models.
Compares predicted essential genes from FBA against experimental dataset.
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from cobra.io import read_sbml_model
import cobra.flux_analysis
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


def apply_standard_medium(model):
    """
    Apply standard minimal medium to model (salts and inorganic nutrients).
    
    Args:
        model: COBRApy model object
        
    Returns:
        Modified model (in-place)
    """
    # Close all exchanges first
    for rxn in model.exchanges:
        rxn.lower_bound = 0.0
    
    # Define standard minimal medium (no organic carbon)
    standard_medium = {
        'EX_o2_e':  -18,        # Oxygen
        'EX_h2o_e': -1000,      # Water
        'EX_pi_e':  -5,         # Inorganic phosphate
        'EX_so4_e': -5,         # Sulfate
        'EX_nh4_e': -5,         # Ammonium
        'EX_k_e':   -1000,       # Potassium
        'EX_na1_e': -1000,     # Sodium
        'EX_ca2_e': -1000,     # Calcium
        'EX_mg2_e': -1000,     # Magnesium
        'EX_fe3_e': -1000,     # Ferric iron
        'EX_fe2_e': -1000,     # Iron
        'EX_mn2_e': -1000,     # Manganese
        'EX_h_e':   -1000,       # Hydrogen
        'EX_co2_e': -1000,     # Carbon dioxide
        'EX_glc__D_e': -10      # Glucose
    }
    
    # Apply medium constraints
    for rxn_id, bound in standard_medium.items():
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = bound
    
    return model


def normalize_locus_tag(locus_tag):
    """
    Normalize locus tags for flexible matching.
    Removes underscores and converts to uppercase.
    
    Args:
        locus_tag: Locus tag string (e.g., "BSU_00010" or "BSU00010")
        
    Returns:
        Normalized locus tag (e.g., "BSU00010")
    """
    if pd.isna(locus_tag):
        return ""
    return str(locus_tag).replace("_", "").replace("-", "").strip().upper()


def parse_essentiality(value):
    """
    Parse essentiality value from dataset.
    
    Args:
        value: Value from 'essential' column
        
    Returns:
        Boolean: True if essential, False otherwise
    """
    if pd.isna(value):
        return False
    
    value_str = str(value).lower().strip()
    
    # Consider as essential
    if value_str in ['yes', 'y', 'true', '1', 'essential']:
        return True
    
    # Everything else (including 'no', empty, etc.) is non-essential
    return False


def calculate_metrics(tp, fp, tn, fn, total_genes_model, total_genes_dataset):
    """
    Calculate classification metrics from confusion matrix.
    
    Args:
        tp: True positives
        fp: False positives
        tn: True negatives
        fn: False negatives
        total_genes_model: Total number of genes in the model
        total_genes_dataset: Total number of genes in the dataset
        
    Returns:
        Dictionary with calculated metrics
    """
    metrics = {}
    
    # TPR - True Positive Rate (Sensitivity, Recall)
    if (tp + fn) > 0:
        metrics['tpr'] = tp / (tp + fn)
        metrics['recall'] = tp / (tp + fn)  # Same as TPR
    else:
        metrics['tpr'] = 0.0
        metrics['recall'] = 0.0
    
    # TNR - True Negative Rate (Specificity)
    if (tn + fp) > 0:
        metrics['tnr'] = tn / (tn + fp)
        metrics['specificity'] = tn / (tn + fp)  # Same as TNR
    else:
        metrics['tnr'] = 0.0
        metrics['specificity'] = 0.0
    
    # FDR - False Discovery Rate
    if (fp + tp) > 0:
        metrics['fdr'] = fp / (fp + tp)
    else:
        metrics['fdr'] = 0.0
    
    # Precision (Positive Predictive Value)
    if (tp + fp) > 0:
        metrics['precision'] = tp / (tp + fp)
    else:
        metrics['precision'] = 0.0
    
    # F1 Score
    if metrics['precision'] + metrics['recall'] > 0:
        metrics['f1_score'] = 2 * (metrics['precision'] * metrics['recall']) / \
                              (metrics['precision'] + metrics['recall'])
    else:
        metrics['f1_score'] = 0.0
    
    # Accuracy
    total = tp + fp + tn + fn
    if total > 0:
        metrics['accuracy'] = (tp + tn) / total
    else:
        metrics['accuracy'] = 0.0
    
    # MCC - Matthews Correlation Coefficient
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator > 0:
        metrics['mcc'] = (tp * tn - fp * fn) / denominator
    else:
        metrics['mcc'] = 0.0
    
    # Coverage - Number of genes in model / Total number of genes in database
    if total_genes_dataset > 0:
        metrics['coverage'] = total_genes_model / total_genes_dataset
    else:
        metrics['coverage'] = 0.0
    
    return metrics


def plot_metrics_comparison(df_summary, outdir, output_name, verbose=True):
    """
    Create grouped bar plot comparing metrics across models.
    
    Args:
        df_summary: DataFrame with summary metrics
        outdir: Output directory for plot
        output_name: Base name for output file
        verbose: Print progress messages
        
    Returns:
        Path to saved plot file
    """
    if df_summary.empty:
        if verbose:
            print("[WARNING] No data to plot")
        return None
    
    # Select metrics to plot (matching the paper figure)
    metrics = ['Precision', 'FDR', 'Recall', 'Specificity', 'MCC', 'Coverage']
    
    # Prepare data for plotting
    models = df_summary['model'].tolist()
    x = np.arange(len(metrics))
    width = 0.8 / len(models)  # Bar width adjusted by number of models
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Define colors for each model (cycle through if more than 4 models)
    plt.style.use(['science','grid' , 'no-latex'])
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#C73E1D', '#3E8914', '#8A1C7C', '#1B5299', '#F0C808']
    
    # Plot bars for each model
    for i, (idx, row) in enumerate(df_summary.iterrows()):
        model_name = row['model']
        values = [row[metric] for metric in metrics]
        offset = (i - len(models)/2 + 0.5) * width
        ax.bar(x + offset, values, width, label=model_name, color=colors[i % len(colors)])
    
    # Customize plot
    ax.set_ylabel('Prediction score', fontsize=12)
    ax.set_xlabel('')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=11)
    ax.set_ylim(0, 1.0)
    ax.legend(loc='upper right', frameon=True, fontsize=10)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add tight layout
    plt.tight_layout()
    
    # Save plot
    plot_path = os.path.join(outdir, f"{output_name}_metrics_comparison.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    if verbose:
        print(f"[INFO] Metrics comparison plot saved to: {plot_path}")
    
    return plot_path

def gene_essentiality_comparison(
            models, dataset_path, locus_column, essential_column,
            outdir, output_name, processes=4, verbose=True,
            force_overwrite=False
        ):
    """
    Compare predicted essential genes against experimental dataset.
    
    Args:
        models: List of paths to SBML model files
        dataset_path: Path to Excel file with experimental essentiality data
        locus_column: Column name for locus tags in dataset
        essential_column: Column name for essentiality values
        outdir: Output directory for results
        output_name: Base name for output files
        processes: Number of parallel processes for essentiality analysis
        verbose: Print detailed progress messages
        force_overwrite: Overwrite existing output files
        
    Returns:
        Dictionary with results and metrics
        
    Raises:
        FileNotFoundError: If dataset file does not exist
        ValueError: If required columns are missing
        FileExistsError: If output exists and force_overwrite is False
    """
    # Load experimental dataset
    if verbose:
        print(f"\n[INFO] Loading experimental dataset: {dataset_path}")
    
    if not os.path.exists(dataset_path):
        raise FileNotFoundError(f"Dataset file not found: {dataset_path}")
    
    try:
        df_exp = pd.read_excel(dataset_path)
    except Exception as e:
        raise RuntimeError(f"Error reading dataset {dataset_path}: {e}")
    
    # Validate required columns
    if locus_column not in df_exp.columns:
        raise ValueError(f"Dataset must contain '{locus_column}' column")
    if essential_column not in df_exp.columns:
        raise ValueError(f"Dataset must contain '{essential_column}' column")
    
    # Parse experimental data
    experimental_data = {}
    for _, row in df_exp.iterrows():
        locus = row[locus_column]
        essential_val = row[essential_column]
        
        # Normalize locus tag
        normalized = normalize_locus_tag(locus)
        if not normalized:
            continue
        
        # Parse essentiality
        is_essential = parse_essentiality(essential_val)
        experimental_data[normalized] = is_essential
    
    if verbose:
        exp_essential_count = sum(experimental_data.values())
        exp_nonessential_count = len(experimental_data) - exp_essential_count
        print(f"[INFO] Experimental data loaded:")
        print(f"  Total genes: {len(experimental_data)}")
        print(f"  Essential: {exp_essential_count}")
        print(f"  Non-essential: {exp_nonessential_count}")
    
    # Check output
    os.makedirs(outdir, exist_ok=True)
    output_file = os.path.join(outdir, f"{output_name}.xlsx")
    if os.path.exists(output_file) and not force_overwrite:
        raise FileExistsError(f"Output file already exists: {output_file}. Use --force to overwrite.")
    
    # Results storage
    all_results = []
    summary_data = []
    
    # Process each model
    for model_path in models:
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        
        if verbose:
            print(f"\n[INFO] Processing model: {model_name}")
        
        # Load model
        try:
            model = load_model(model_path)
        except Exception as e:
            print(f"[ERROR] Could not load model {model_path}: {e}", file=sys.stderr)
            continue
        
        # Apply minimal medium
        apply_standard_medium(model)
        
        if verbose:
            print(f"[INFO] Total genes in model: {len(model.genes)}")
        
        # Find essential genes using COBRApy
        try:
            if verbose:
                print(f"[INFO] Running essentiality analysis (processes={processes})...")
            
            essential_genes_obj = cobra.flux_analysis.find_essential_genes(
                model, processes=processes
            )
            predicted_essential_ids = {gene.id for gene in essential_genes_obj}
            
            if verbose:
                print(f"[INFO] Predicted essential genes: {len(predicted_essential_ids)}")
                
        except Exception as e:
            print(f"[ERROR] Essentiality analysis failed for {model_name}: {e}", file=sys.stderr)
            continue
        
        # Build mapping: normalized → original gene ID
        model_genes_map = {}
        for gene in model.genes:
            normalized = normalize_locus_tag(gene.id)
            if normalized:
                model_genes_map[normalized] = gene.id
        
        # Normalize predicted essential IDs
        predicted_essential_normalized = {normalize_locus_tag(gid) for gid in predicted_essential_ids}
        
        # Find common genes between dataset and model
        common_genes = set(experimental_data.keys()) & set(model_genes_map.keys())
        
        if verbose:
            print(f"[INFO] Common genes (dataset and model): {len(common_genes)}")
        
        if len(common_genes) == 0:
            print(f"[WARNING] No common genes found between dataset and {model_name}", 
                  file=sys.stderr)
            continue
        
        # Calculate confusion matrix
        tp, fp, tn, fn = 0, 0, 0, 0
        
        for gene_normalized in common_genes:
            exp_essential = experimental_data[gene_normalized]
            pred_essential = gene_normalized in predicted_essential_normalized
            
            # Determine classification
            if exp_essential and pred_essential:
                tp += 1
            elif not exp_essential and pred_essential:
                fp += 1
            elif not exp_essential and not pred_essential:
                tn += 1
            elif exp_essential and not pred_essential:
                fn += 1
        
        # Calculate metrics
        metrics = calculate_metrics(tp, fp, tn, fn, len(model.genes), len(experimental_data))
        
        if verbose:
            print(f"[INFO] Confusion Matrix:")
            print(f"  TP: {tp}, FP: {fp}")
            print(f"  FN: {fn}, TN: {tn}")
            print(f"[INFO] Metrics:")
            print(f"  TPR (Sensitivity/Recall): {metrics['tpr']:.4f}")
            print(f"  TNR (Specificity): {metrics['tnr']:.4f}")
            print(f"  FDR: {metrics['fdr']:.4f}")
            print(f"  Precision: {metrics['precision']:.4f}")
            print(f"  F1-Score: {metrics['f1_score']:.4f}")
            print(f"  Accuracy: {metrics['accuracy']:.4f}")
            print(f"  MCC: {metrics['mcc']:.4f}")
            print(f"  Coverage: {metrics['coverage']:.4f}")
        
        # Store summary
        summary_data.append({
            'model': model_name,
            'genes_in_model': len(model.genes),
            'genes_in_dataset': len(experimental_data),
            'common_genes': len(common_genes),
            'TP': tp,
            'FP': fp,
            'TN': tn,
            'FN': fn,
            'TPR': round(metrics['tpr'], 4),
            'TNR': round(metrics['tnr'], 4),
            'FDR': round(metrics['fdr'], 4),
            'MCC': round(metrics['mcc'], 4),
            'Precision': round(metrics['precision'], 4),
            'Recall': round(metrics['recall'], 4),
            'Specificity': round(metrics['specificity'], 4),
            'F1_Score': round(metrics['f1_score'], 4),
            'Accuracy': round(metrics['accuracy'], 4),
            'Coverage': round(metrics['coverage'], 4)
        })
    
    # Create output DataFrame
    df_summary = pd.DataFrame(summary_data)
    
    # Save to Excel - only Summary sheet
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df_summary.to_excel(writer, sheet_name='Summary', index=False)
    
    # GENERAR EL BARPLOT - Esta es la línea que faltaba
    if not df_summary.empty:
        plot_path = plot_metrics_comparison(df_summary, outdir, output_name, verbose)
    else:
        plot_path = None
    
    if verbose:
        print(f"\n[DONE] Results saved to: {output_file}")
        print(f"[INFO] Sheet created: Summary")
        if plot_path:
            print(f"[INFO] Plot saved to: {plot_path}")
    
    return {
        'summary': df_summary,
        'output_file': output_file,
        'plot_file': plot_path
    }


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Compare predicted essential genes against experimental dataset",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--config", type=str, default="config.yaml",
                       help="Configuration YAML file")
    parser.add_argument("--models", "-m", nargs="+",
                       help="List of model files to test (overrides YAML)")
    parser.add_argument("--dataset", type=str,
                       help="Excel file with experimental essentiality data (overrides YAML)")
    parser.add_argument("--locus-column", type=str,
                       help="Column name for locus tags (overrides YAML)")
    parser.add_argument("--essential-column", type=str,
                       help="Column name for essentiality values (overrides YAML)")
    parser.add_argument("--outdir", type=str,
                       help="Output directory (overrides YAML)")
    parser.add_argument("--output-name", type=str,
                       help="Base name for output files (overrides YAML)")
    parser.add_argument("--processes", "-p", type=int,
                       help="Number of parallel processes (overrides YAML)")
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
        if args.verbose:
            print(f"[WARNING] Config file {args.config} not found, using CLI arguments only")
    
    cfg_section = config.get("comparison", {}).get("gene_essentiality") or {}
    cfg_section2=config.get("comparison") or {}
    # Priority: CLI > config.yaml    
    
    # Priority: CLI > config.yaml > defaults
    models = args.models if args.models else cfg_section2.get("models", [])
    dataset_path = args.dataset if args.dataset else cfg_section.get(
        "dataset", 
        r"C:\Users\felip\Desktop\python\GEMcompare\Pipeline\2_Comparison\Essentiality_dataset.xlsx"
    )
    locus_column = args.locus_column if args.locus_column else cfg_section.get(
        "locus_column", "genomic_annotation.locus_tag"
    )
    essential_column = args.essential_column if args.essential_column else cfg_section.get(
        "essential_column", "essential"
    )
    outdir = args.outdir if args.outdir else cfg_section.get(
        "outdir", "outputs/comparison/gene_essentiality"
    )
    output_name = args.output_name if args.output_name else cfg_section.get(
        "output_name", "GeneEssentiality"
    )
    processes = args.processes if args.processes else cfg_section.get("processes", 4)
    verbose = args.verbose or cfg_section.get("verbose", False)
    force_overwrite = args.force or cfg_section.get("force", False)
    
    # Validate required parameters
    if not models:
        parser.error("No models specified. Use --models or define 'comparison.gene_essentiality.models' in config.yaml")
    
    if verbose:
        print("[INFO] Gene essentiality comparison configuration:")
        print(f"  Models: {len(models)}")
        for m in models:
            print(f"    - {m}")
        print(f"  Dataset: {dataset_path}")
        print(f"  Locus column: {locus_column}")
        print(f"  Essential column: {essential_column}")
        print(f"  Output directory: {outdir}")
        print(f"  Processes: {processes}")
    
    # Run comparison
    try:
        results = gene_essentiality_comparison(
            models, dataset_path, locus_column, essential_column,
            outdir, output_name, processes=processes,
            verbose=verbose, force_overwrite=force_overwrite
        )
        
        if verbose and results.get('plot_file'):
            print(f"[SUCCESS] Analysis completed. Plot generated: {results['plot_file']}")
            
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())