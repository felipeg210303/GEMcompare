#!/usr/bin/env python3
"""
2.3_Diff_CarbonSources.py

Compare predicted growth (FBA) across different carbon sources for one or more models.

Usage Modes:
  - Pipeline mode: Import as module and call run_carbon_sources() function
  - Standalone mode: python "2.3_Diff_CarbonSources.py" --config config.yaml --model model1.xml --model model2.xml

Core Functionality:
  - Performs FBA optimization for multiple models across specified carbon sources
  - Generates CSV results and comparative bar plots
  - Supports custom carbon source definitions via configuration

Dependencies:
  - cobrapy: For SBML model reading and FBA optimization
  - pandas: For data manipulation and CSV export
  - matplotlib: For visualization and plot generation
  - pyyaml: For configuration file parsing
"""

import os
import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt

try:
    import yaml
except ImportError:
    yaml = None

try:
    from cobra.io import read_sbml_model
except ImportError:
    read_sbml_model = None


# -------------------------------
# Configuration Constants
# -------------------------------
DEFAULT_SOURCES = [
    "glucose", "pyruvate", "lactate", "fructose", "malate", 
    "sucrose", "gluconate", "glycerol", "succinate", "ethanol"
]

FRIENDLY_SOURCES = {
    # Azúcares simples
    "glucose":    "EX_glc__D_e",
    "fructose":   "EX_fru_e",
    "galactose":  "EX_gal_e",
    "mannose":    "EX_man_e",
    "ribose":     "EX_rib__D_e",
    "xylose":     "EX_xyl__D_e",
    "arabinose":  "EX_arab__L_e",
    
    # Azúcares complejos
    "sucrose":    "EX_sucr_e",  # Fixed typo: was EX_suc_e in original
    "maltose":    "EX_malt_e",
    "lactose":    "EX_lcts_e",
    "trehalose":  "EX_tre_e",
    "cellobiose": "EX_cellb_e",
    
    # Alcoholes
    "glycerol":   "EX_glyc_e",  # Fixed: was EX_glycerol_e in original
    "ethanol":    "EX_etoh_e",
    "mannitol":   "EX_mnl_e",
    "sorbitol":   "EX_sbt__D_e",
    
    # Ácidos orgánicos
    "pyruvate":   "EX_pyr_e",
    "lactate":    "EX_lac__L_e",
    "malate":     "EX_mal__L_e",
    "succinate":  "EX_succ_e",
    "acetate":    "EX_ac_e",
    "citrate":    "EX_cit_e",
    "fumarate":   "EX_fum_e",
    "alpha_ketoglutarate": "EX_akg_e",
    "formate":    "EX_for_e",
    
    # Otros compuestos
    "gluconate":  "EX_glcn__D_e",
    "glutamate":  "EX_glu__L_e",
    "glutamine":  "EX_gln__L_e",
    "aspartate":  "EX_asp__L_e",
    "alanine":    "EX_ala__L_e",
    
    # Fuentes complejas
    "starch":     "EX_starch_e",
    
    # Nucleótidos y derivados
    "uracil":     "EX_ura_e",  # Fixed: was EX_urac_e in original
    "adenine":    "EX_ade_e",
    "guanine":    "EX_gua_e",
    
    # Amino azúcares
    "glucosamine": "EX_gam_e",
    "N_acetylglucosamine": "EX_acgam_e"  # Fixed: was EX_glcNac_e in original
}


# -------------------------------
# Configuration Loading
# -------------------------------
def load_config(config_path: str) -> dict:
    """
    Load and parse YAML configuration file.
    
    Args:
        config_path (str): Path to the YAML configuration file
        
    Returns:
        dict: Parsed configuration dictionary
        
    Raises:
        FileNotFoundError: If config file does not exist
        RuntimeError: If PyYAML is not installed
    """
    if yaml is None:
        raise RuntimeError("PyYAML is not installed. Install it with: pip install pyyaml")
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
        
    with open(config_path, "r") as config_file:
        config_data = yaml.safe_load(config_file) or {}
    return config_data


def parse_custom_sources(custom_sources_list: list) -> dict:
    """
    Parse custom carbon source definitions from command line arguments.
    
    Args:
        custom_sources_list (list): List of 'name=token' strings
        
    Returns:
        dict: Mapping of friendly names to exchange reaction tokens
        
    Raises:
        ValueError: If any entry doesn't contain '=' separator
    """
    mapping = {}
    if not custom_sources_list:
        return mapping
        
    for entry in custom_sources_list:
        if "=" not in entry:
            raise ValueError(f"Invalid custom source entry (expected name=token): {entry}")
        name, token = entry.split("=", 1)
        mapping[name.strip()] = token.strip()
    return mapping


def list_available_sources() -> str:
    """
    Generate formatted list of available carbon sources.
    
    Returns:
        str: Formatted string listing all available sources
    """
    lines = ["Available carbon sources (friendly name -> exchange reaction):"]
    for name, token in FRIENDLY_SOURCES.items():
        default_tag = " (default)" if name in DEFAULT_SOURCES else ""
        lines.append(f"  {name:25s} -> {token}{default_tag}")
    return "\n".join(lines)


# -------------------------------
# Model Processing Functions
# -------------------------------
def find_exchange_reaction(model, source_token: str):
    """
    Find the appropriate exchange reaction for a given carbon source.
    
    Args:
        model: COBRA model object
        source_token (str): Token identifying the carbon source
        
    Returns:
        Reaction or None: Exchange reaction if found, None otherwise
    """
    token = source_token.lower()

    # Direct reaction ID match
    if source_token in model.reactions:
        rxn = model.reactions.get_by_id(source_token)
        if rxn in model.exchanges:
            return rxn

    # Metabolite ID match
    matching_metabolites = [met for met in model.metabolites if token == met.id.lower()]
    if matching_metabolites:
        for met in matching_metabolites:
            for rxn in model.exchanges:
                if met in rxn.metabolites:
                    return rxn

    # Partial matches in reaction ID, name, or metabolite names
    for rxn in model.exchanges:
        if token in rxn.id.lower():
            return rxn
            
    for rxn in model.exchanges:
        if rxn.name and token in rxn.name.lower():
            return rxn
            
    for rxn in model.exchanges:
        for met in rxn.metabolites:
            if token in met.id.lower() or (met.name and token in met.name.lower()):
                return rxn
                
    return None


def configure_uptake(reaction, uptake_value: float):
    """
    Configure uptake bounds for an exchange reaction.
    
    Args:
        reaction: COBRA reaction object
        uptake_value (float): Uptake rate value (absolute value used)
    """
    uptake = float(abs(uptake_value))
    reaction.lower_bound = -uptake


# -------------------------------
# Core Workflow Function
# -------------------------------
def run_carbon_sources(models: list, sources_mapping: dict, uptake: float, 
                      output_dir: str, output_name: str, verbose: bool = True, 
                      force: bool = False) -> dict:
    """
    Main workflow: test carbon sources across models and generate results.
    
    Args:
        models (list): List of paths to SBML model files
        sources_mapping (dict): Mapping of friendly names to exchange tokens
        uptake (float): Carbon uptake rate for FBA
        output_dir (str): Directory for output files
        output_name (str): Base name for output files
        verbose (bool): Enable detailed progress output
        force (bool): Overwrite existing output files
        
    Returns:
        dict: Results dictionary containing paths and dataframe
        
    Raises:
        FileExistsError: If output files exist and force is False
    """
    if read_sbml_model is None:
        raise RuntimeError("cobrapy is not installed. Install it with: pip install cobra")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    records = []

    for model_path in models:
        model_abspath = os.path.abspath(model_path)
        if verbose:
            print(f"\n[model] Loading {model_abspath} ...")
            
        try:
            base_model = read_sbml_model(model_abspath)
        except Exception as e:
            print(f"[ERROR] Could not load model {model_path}: {e}", file=sys.stderr)
            # Record failure for all sources for this model
            for friendly_name in sources_mapping.keys():
                records.append({
                    "model": os.path.splitext(os.path.basename(model_path))[0],
                    "carbon_source": friendly_name,
                    "growth_rate": float("nan"),
                    "carbon_uptake_flux": float("nan"),
                    "note": "model_load_failed"
                })
            continue

        # Test each carbon source
        for friendly_name, exchange_token in sources_mapping.items():
            if verbose:
                print(f"[source] Testing '{friendly_name}' (token: '{exchange_token}')")

            # Create model copy for isolation
            model = base_model.copy()

            # Disable glucose uptake if present
            if "EX_glc__D_e" in model.reactions:
                model.reactions.get_by_id("EX_glc__D_e").lower_bound = 0.0

            # Find exchange reaction
            exchange_rxn = find_exchange_reaction(model, exchange_token)
            if exchange_rxn is None:
                if verbose:
                    print(f"[WARN] Exchange for token '{exchange_token}' not found")
                records.append({
                    "model": os.path.splitext(os.path.basename(model_path))[0],
                    "carbon_source": friendly_name,
                    "growth_rate": float("nan"),
                    "carbon_uptake_flux": float("nan"),
                    "note": "exchange_not_found"
                })
                continue

            # Configure uptake and optimize
            configure_uptake(exchange_rxn, uptake)
            
            try:
                solution = model.optimize()
                growth = solution.objective_value if solution else float("nan")
                flux_val = float(solution.fluxes.get(exchange_rxn.id, float("nan"))) if solution else float("nan")
                
                records.append({
                    "model": os.path.splitext(os.path.basename(model_path))[0],
                    "carbon_source": friendly_name,
                    "growth_rate": growth,
                    "carbon_uptake_flux": flux_val,
                    "note": ""
                })
                
                if verbose:
                    print(f"  -> growth: {growth:.4f}, {exchange_rxn.id} flux: {flux_val:.4f}")
                    
            except Exception as e:
                print(f"[ERROR] Optimization failed for {model_path} on source {friendly_name}: {e}", file=sys.stderr)
                records.append({
                    "model": os.path.splitext(os.path.basename(model_path))[0],
                    "carbon_source": friendly_name,
                    "growth_rate": float("nan"),
                    "carbon_uptake_flux": float("nan"),
                    "note": "opt_failed"
                })

    # Generate output files
    df = pd.DataFrame(records)
    csv_path = os.path.join(output_dir, f"{output_name}.csv")
    plot_path = os.path.join(output_dir, f"{output_name}_barplot.png")

    # Check for existing files
    if os.path.exists(csv_path) and not force:
        raise FileExistsError(f"CSV output exists: {csv_path}. Use --force to overwrite.")

    # Save CSV results
    df.to_csv(csv_path, index=False)
    if verbose:
        print(f"\nResults saved to: {csv_path}")

    # Generate comparative plot
    plot_success = False
    try:
        pivot_data = df.pivot_table(
            index="carbon_source", 
            columns="model", 
            values="growth_rate", 
            aggfunc="first"
        )
        # Maintain source order
        pivot_data = pivot_data.reindex(list(sources_mapping.keys()))
        
        axis = pivot_data.plot(kind="bar", figsize=(12, 7))
        axis.set_ylabel("Growth rate (objective value)")
        axis.set_xlabel("Carbon source")
        axis.set_title("Predicted growth rates by carbon source and model")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()
        plot_success = True
        if verbose:
            print(f"Plot saved to: {plot_path}")
            
    except Exception as e:
        print(f"[WARN] Could not create visualization: {e}", file=sys.stderr)
        plot_path = None

    return {
        "csv": csv_path,
        "plot": plot_path if plot_success else None,
        "dataframe": df
    }


def main_carbon_sources(model_paths: list, configuration: dict) -> dict:
    """
    Pipeline-friendly interface for carbon source comparison.
    
    Args:
        model_paths (list): List of paths to SBML model files
        configuration (dict): Configuration parameters including:
            - sources: Carbon sources to test
            - uptake: Uptake rate value
            - output_dir: Output directory path
            - output_name: Base name for output files
            - verbose: Enable detailed output
            - force: Overwrite existing files
            
    Returns:
        dict: Results dictionary from run_carbon_sources()
    """
    # Extract configuration with defaults
    carbon_config = configuration.get("comparison", {}).get("carbon_sources", {})
    
    sources_spec = carbon_config.get("sources", DEFAULT_SOURCES)
    uptake_rate = carbon_config.get("uptake", 5.0)
    output_directory = carbon_config.get("output_dir", "outputs/comparison/carbon_sources")
    output_basename = carbon_config.get("output_name", "carbon_sources")
    verbose_output = carbon_config.get("verbose", True)
    force_overwrite = carbon_config.get("force", False)

    # Build sources mapping
    sources_mapping = {}
    for source in sources_spec:
        if isinstance(source, dict):
            # Handle {name: token} format
            sources_mapping.update(source)
        elif source in FRIENDLY_SOURCES:
            # Use predefined token
            sources_mapping[source] = FRIENDLY_SOURCES[source]
        else:
            # Use source name as token
            sources_mapping[source] = source

    return run_carbon_sources(
        models=model_paths,
        sources_mapping=sources_mapping,
        uptake=uptake_rate,
        output_dir=output_directory,
        output_name=output_basename,
        verbose=verbose_output,
        force=force_overwrite
    )


# -------------------------------
# Command Line Interface
# -------------------------------
def main():
    """
    Standalone command-line interface for carbon source comparison.
    """
    # Build help epilog
    epilog_text = ("Available carbon sources (default and others):\n" +
                  list_available_sources() +
                  "\n\nCustom sources: use --add-source name=token (repeatable). "
                  "CLI sources override YAML and default configurations.")

    parser = argparse.ArgumentParser(
        description="Compare model growth across carbon sources using FBA",
        epilog=epilog_text,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("--config", "-c", type=str, default="config.yaml",
                       help="YAML configuration file path")

    parser.add_argument("--models", "-m", nargs="+", 
                       help="Lista de modelos a comparar (sobre-escribe el YAML)")
    parser.add_argument("--sources", "-s", nargs="+",
                       help="Carbon sources to test (names or name=token pairs)")
    parser.add_argument("--add-source", action="append",
                       help="Add custom carbon source: name=token")
    parser.add_argument("--list-sources", action="store_true",
                       help="List available carbon sources and exit")
    parser.add_argument("--uptake", "-u", type=float,
                       help="Carbon uptake rate (default: 5.0)")
    parser.add_argument("--outdir", type=str,
                       help="Output directory for results")
    parser.add_argument("--output-name", type=str,
                       help="Base name for output files")
    parser.add_argument("--force", action="store_true",
                       help="Overwrite existing output files")
    parser.add_argument("--verbose", action="store_true",
                       help="Enable detailed progress output")

    args = parser.parse_args()

    # Handle list sources option
    if args.list_sources:
        print(list_available_sources())
        sys.exit(0)

    # Validate required dependencies
    if read_sbml_model is None:
        print("ERROR: cobrapy is required but not installed.", file=sys.stderr)
        print("Install it with: pip install cobra", file=sys.stderr)
        sys.exit(1)

    # Load configuration
    try:
        config_data = load_config(args.config)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    # Get model list - FIXED LOGIC
    models = args.models or \
             config_data.get("comparison", {}).get("models") or \
             config_data.get("comparison", {}).get("carbon_sources", {}).get("models") or \
             config_data.get("preparation", {}).get("models") or []

    if not models:
        print("ERROR: No models provided via --model or config file.", file=sys.stderr)
        print("Please specify models using --model or in the configuration file.", file=sys.stderr)
        sys.exit(2)

    # Build carbon sources mapping
    carbon_config = config_data.get("comparison", {}).get("carbon_sources", {})
    
    # Start with predefined sources
    sources_mapping = dict(FRIENDLY_SOURCES)
    
    # Add custom sources from CLI
    try:
        custom_sources = parse_custom_sources(args.add_source)
        sources_mapping.update(custom_sources)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    # Determine final source selection
    selected_sources = []
    if args.sources:
        # CLI sources take precedence
        for source in args.sources:
            if "=" in source:
                name, token = source.split("=", 1)
                selected_sources.append((name.strip(), token.strip()))
            elif source in sources_mapping:
                selected_sources.append((source, sources_mapping[source]))
            else:
                selected_sources.append((source, source))
    elif carbon_config.get("sources"):
        # Use YAML configuration
        yaml_sources = carbon_config["sources"]
        for source in yaml_sources:
            if isinstance(source, dict):
                for name, token in source.items():
                    selected_sources.append((name, token))
            elif source in sources_mapping:
                selected_sources.append((source, sources_mapping[source]))
            else:
                selected_sources.append((source, source))
    else:
        # Use defaults
        selected_sources = [(name, FRIENDLY_SOURCES[name]) for name in DEFAULT_SOURCES]

    final_mapping = {name: token for name, token in selected_sources}

    # Determine parameters (CLI overrides YAML)
    uptake_rate = args.uptake or carbon_config.get("uptake", 5.0)
    output_directory = args.outdir or carbon_config.get("output_dir", "outputs/comparison/carbon_sources")
    output_basename = args.output_name or carbon_config.get("output_name", "carbon_sources")
    verbose_output = args.verbose if args.verbose is not None else carbon_config.get("verbose", True)
    force_overwrite = args.force or carbon_config.get("force", True)

    # Display configuration
    if verbose_output:
        print("Carbon source comparison configuration:")
        print(" Models:", models)
        print(" Selected sources (friendly -> token):")
        for name, token in final_mapping.items():
            print(f"   {name:20} -> {token}")
        print(" Uptake rate:", uptake_rate)
        print(" Output directory:", output_directory)
        print(" Output base name:", output_basename)
        print(" Force overwrite:", force_overwrite)

    # Execute comparison
    try:
        results = run_carbon_sources(
            models=models,
            sources_mapping=final_mapping,
            uptake=uptake_rate,
            output_dir=output_directory,
            output_name=output_basename,
            verbose=verbose_output,
            force=force_overwrite
        )
        print(f"\nComparison completed successfully.")
        print(f"Results: {results['csv']}")
        if results['plot']:
            print(f"Plot: {results['plot']}")
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()