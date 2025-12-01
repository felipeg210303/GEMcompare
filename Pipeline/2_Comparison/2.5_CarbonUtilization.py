#!/usr/bin/env python3
"""
Biolog carbon source utilization comparison for GEM models.
Compares predicted growth against experimental Biolog data.
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
from multiprocessing import Pool


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


def parse_experimental_result(value):
    """
    Parse experimental result from Biolog data.
    
    Args:
        value: Value from 'result' column
        
    Returns:
        - True: if "+"
        - False: if "-" or empty/NaN
        - None: if "x" (ambiguous, excluded from analysis)
    """
    if pd.isna(value):
        return False
    
    value_str = str(value).strip()
    
    if value_str == "+":
        return True
    elif value_str == "-":
        return False
    elif value_str == "x":
        return None
    else:
        return False


def test_single_carbon_source(args_tuple):
    """
    Test growth on a single carbon source (designed for parallel execution).
    
    Args:
        args_tuple: Tuple of (model_path, carbon_id, exp_growth, carbon_exchanges, uptake_rate, threshold)
        
    Returns:
        Tuple of (carbon_id, exp_growth, pred_growth, success)
    """
    model_path, carbon_id, exp_growth, carbon_exchanges, uptake_rate, threshold = args_tuple
    
    try:
        # Load fresh model copy
        model = read_sbml_model(model_path)
        
        # Close all carbon sources
        for carbon_rxn_id in carbon_exchanges:
            if carbon_rxn_id in model.reactions:
                model.reactions.get_by_id(carbon_rxn_id).lower_bound = 0.0
        
        # Check if exchange reaction exists
        if carbon_id not in model.reactions:
            return (carbon_id, exp_growth, False, False)
        
        # Get exchange reaction
        carbon_rxn = model.reactions.get_by_id(carbon_id)
        
        # Verify it's an exchange reaction
        if carbon_rxn not in model.exchanges:
            return (carbon_id, exp_growth, False, False)
        
        # Open current carbon source
        carbon_rxn.lower_bound = -abs(uptake_rate)
        
        # Run FBA
        solution = model.optimize()
        
        if solution.status == "optimal":
            growth_rate = solution.objective_value
            pred_growth = growth_rate >= threshold
        else:
            pred_growth = False
        
        return (carbon_id, exp_growth, pred_growth, True)
        
    except Exception:
        return (carbon_id, exp_growth, False, False)


def calculate_metrics(tp, fp, tn, fn, total_sources_model, total_sources_dataset):
    """
    Calculate classification metrics from confusion matrix.
    
    Args:
        tp: True positives
        fp: False positives
        tn: True negatives
        fn: False negatives
        total_sources_model: Total carbon sources available in model
        total_sources_dataset: Total carbon sources in experimental dataset
        
    Returns:
        Dictionary with calculated metrics
    """
    metrics = {}
    
    # TPR - True Positive Rate (Sensitivity, Recall)
    if (tp + fn) > 0:
        metrics['tpr'] = tp / (tp + fn)
        metrics['recall'] = tp / (tp + fn)
    else:
        metrics['tpr'] = 0.0
        metrics['recall'] = 0.0
    
    # TNR - True Negative Rate (Specificity)
    if (tn + fp) > 0:
        metrics['tnr'] = tn / (tn + fp)
        metrics['specificity'] = tn / (tn + fp)
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
    
    # Coverage
    if total_sources_dataset > 0:
        metrics['coverage'] = total_sources_model / total_sources_dataset
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
    
    # Select metrics to plot
    metrics = ['Precision', 'FDR', 'TPR', 'TNR', 'MCC', 'Coverage']
    
    # Prepare data for plotting
    models = df_summary['model'].tolist()
    x = np.arange(len(metrics))
    width = 0.8 / len(models)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Define colors for each model
    colors = ['#4472C4', '#ED7D31', '#A5A5A5', '#FFC000', '#5B9BD5', '#70AD47']
    
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
    
    plt.tight_layout()
    
    # Save plot
    plot_path = os.path.join(outdir, f"{output_name}_metrics_comparison.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    if verbose:
        print(f"[INFO] Metrics comparison plot saved to: {plot_path}")
    
    return plot_path


def biolog_comparison(models, carbon_file, outdir, output_name, uptake_rate,
                     threshold, processes=4, verbose=True, force_overwrite=False):
    """
    Compare predicted carbon source utilization against experimental Biolog data.
    
    Args:
        models: List of paths to SBML model files
        carbon_file: Path to Excel file with experimental Biolog data
        outdir: Output directory for results
        output_name: Base name for output files
        uptake_rate: Maximum carbon uptake rate (positive value)
        threshold: Minimum growth rate to classify as positive
        processes: Number of parallel processes
        verbose: Print detailed progress messages
        force_overwrite: Overwrite existing output files
        
    Returns:
        Dictionary with results and metrics
        
    Raises:
        FileNotFoundError: If carbon_file does not exist
        ValueError: If required columns are missing
        FileExistsError: If output exists and force_overwrite is False
    """
    # Load experimental data
    if verbose:
        print(f"\n[INFO] Loading experimental Biolog data: {carbon_file}")
    
    if not os.path.exists(carbon_file):
        raise FileNotFoundError(f"Biolog data file not found: {carbon_file}")
    
    try:
        df_exp = pd.read_excel(carbon_file, sheet_name='Names')
    except ValueError as e:
        raise ValueError(f"Excel file must contain a 'Names' sheet: {e}")
    except Exception as e:
        raise RuntimeError(f"Error reading Excel file {carbon_file}: {e}")
    
    # Validate required columns
    if 'ID' not in df_exp.columns:
        raise ValueError("'Names' sheet must contain 'ID' column")
    if 'result' not in df_exp.columns:
        raise ValueError("'Names' sheet must contain 'result' column")
    
    # Parse experimental data
    experimental_data = {}
    ambiguous_sources = []
    
    for _, row in df_exp.iterrows():
        carbon_id = row['ID']
        exp_result = parse_experimental_result(row['result'])
        
        if exp_result is None:
            ambiguous_sources.append(carbon_id)
        else:
            experimental_data[carbon_id] = exp_result
    
    if verbose:
        exp_positive = sum(experimental_data.values())
        exp_negative = len(experimental_data) - exp_positive
        print(f"[INFO] Experimental data loaded:")
        print(f"  Total carbon sources: {len(experimental_data)}")
        print(f"  Positive (+): {exp_positive}")
        print(f"  Negative (-): {exp_negative}")
        print(f"  Ambiguous (x): {len(ambiguous_sources)} (excluded from analysis)")
    
    # Check output
    os.makedirs(outdir, exist_ok=True)
    output_file = os.path.join(outdir, f"{output_name}.xlsx")
    if os.path.exists(output_file) and not force_overwrite:
        raise FileExistsError(f"Output file already exists: {output_file}. Use --force to overwrite.")
    
    # Results storage
    summary_data = []
    
    # Process each model
    for model_path in models:
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        
        if verbose:
            print(f"\n[INFO] Processing model: {model_name}")
        
        # Load model
        try:
            base_model = load_model(model_path)
        except Exception as e:
            print(f"[ERROR] Could not load model {model_path}: {e}", file=sys.stderr)
            continue
        
        # Identify all carbon exchange reactions to close
        carbon_exchanges = []
        for rxn in base_model.exchanges:
            if rxn.metabolites:
                met = list(rxn.metabolites.keys())[0]
                if met.elements.get('C', 0) > 0 and met.compartment == 'e':
                    carbon_exchanges.append(rxn.id)
        
        if verbose:
            print(f"[INFO] Identified {len(carbon_exchanges)} carbon exchange reactions to manage")
        
        # Prepare arguments for parallel processing
        test_args = [
            (model_path, carbon_id, exp_growth, carbon_exchanges, uptake_rate, threshold)
            for carbon_id, exp_growth in experimental_data.items()
        ]
        
        if verbose:
            print(f"[INFO] Testing {len(test_args)} carbon sources using {processes} processes...")
        
        # Run parallel FBA tests
        with Pool(processes=processes) as pool:
            results = pool.map(test_single_carbon_source, test_args)
        
        # Calculate confusion matrix from results
        tp, fp, tn, fn = 0, 0, 0, 0
        tested_sources = 0
        
        for carbon_id, exp_growth, pred_growth, success in results:
            if not success:
                if verbose:
                    print(f"[WARNING] {carbon_id} could not be tested")
                continue
            
            # Update confusion matrix
            if exp_growth and pred_growth:
                tp += 1
            elif not exp_growth and pred_growth:
                fp += 1
            elif not exp_growth and not pred_growth:
                tn += 1
            elif exp_growth and not pred_growth:
                fn += 1
            
            tested_sources += 1
        
        if verbose:
            print(f"[INFO] Tested sources: {tested_sources}/{len(experimental_data)}")
            print(f"[INFO] Confusion Matrix:")
            print(f"  TP: {tp}, FP: {fp}")
            print(f"  FN: {fn}, TN: {tn}")
        
        # Calculate metrics
        total_model_sources = len([rxn for rxn in base_model.exchanges])
        metrics = calculate_metrics(tp, fp, tn, fn, total_model_sources, len(experimental_data))
        
        if verbose:
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
            'sources_in_model': total_model_sources,
            'sources_in_dataset': len(experimental_data),
            'sources_tested': tested_sources,
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
    
    # Save to Excel
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df_summary.to_excel(writer, sheet_name='Summary', index=False)
    
    if verbose:
        print(f"\n[DONE] Results saved to: {output_file}")
        print(f"[INFO] Sheet created: Summary")
    
    # Generate metrics comparison plot
    plot_path = plot_metrics_comparison(df_summary, outdir, output_name, verbose=verbose)
    
    return {
        'summary': df_summary,
        'output_file': output_file,
        'plot': plot_path
    }


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Compare predicted carbon source utilization against experimental Biolog data",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--config", type=str, default="config.yaml",
                       help="Configuration YAML file")
    parser.add_argument("--models", "-m", nargs="+",
                       help="List of model files to test (overrides YAML)")
    parser.add_argument("--carbon-file", type=str,
                       help="Excel file with experimental Biolog data (overrides YAML)")
    parser.add_argument("--outdir", type=str,
                       help="Output directory (overrides YAML)")
    parser.add_argument("--output-name", type=str,
                       help="Base name for output files (overrides YAML)")
    parser.add_argument("--uptake", "-u", type=float,
                       help="Maximum carbon uptake rate (overrides YAML)")
    parser.add_argument("--threshold", "-t", type=float,
                       help="Minimum growth rate for positive result (overrides YAML)")
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
    
    cfg_section = config.get("comparison", {}).get("biolog") or {}
    cfg_section2 = config.get("comparison") or {}
    
    # Priority: CLI > config.yaml > defaults
    models = args.models if args.models else cfg_section2.get("models", [])
    carbon_file = args.carbon_file if args.carbon_file else cfg_section.get(
        "carbon_file", "CarbonUtilization.xlsx"
    )
    outdir = args.outdir if args.outdir else cfg_section.get(
        "outdir", "outputs/comparison/biolog"
    )
    output_name = args.output_name if args.output_name else cfg_section.get(
        "output_name", "BiologComparison"
    )
    uptake_rate = args.uptake if args.uptake else cfg_section.get("uptake", 10.0)
    threshold = args.threshold if args.threshold else cfg_section.get("threshold", 1e-3)
    processes = args.processes if args.processes else cfg_section.get("processes", 4)
    verbose = args.verbose or cfg_section.get("verbose", False)
    force_overwrite = args.force or cfg_section.get("force", False)
    
    # Validate required parameters
    if not models:
        parser.error("No models specified. Use --models or define 'comparison.biolog.models' in config.yaml")
    
    if verbose:
        print("[INFO] Biolog comparison configuration:")
        print(f"  Models: {len(models)}")
        for m in models:
            print(f"    - {m}")
        print(f"  Biolog data file: {carbon_file}")
        print(f"  Output directory: {outdir}")
        print(f"  Uptake rate: {uptake_rate}")
        print(f"  Growth threshold: {threshold}")
        print(f"  Processes: {processes}")
    
    # Run comparison
    try:
        biolog_comparison(models, carbon_file, outdir, output_name,
                         uptake_rate, threshold, processes=processes,
                         verbose=verbose, force_overwrite=force_overwrite)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())