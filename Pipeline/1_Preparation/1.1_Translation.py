#!/usr/bin/env python3
"""
Model translation between namespace formats.
Converts metabolic models between different database namespaces (BiGG, SEED, KEGG).

Supported input formats: SBML (.xml, .sbml), MATLAB (.mat), JSON (.json)
Supported namespaces: BiGG, SEED, KEGG, and others supported by mergem
"""

import os
import argparse
import sys
import yaml
import cobra.io as io
from cobra.io import write_sbml_model
from mergem import translate


def load_model(model_path):
    """
    Load metabolic model with automatic format detection.
    
    Args:
        model_path: Path to model file
        
    Returns:
        Loaded cobra.Model object
        
    Raises:
        FileNotFoundError: If model file does not exist
        RuntimeError: If file format is not supported
    """
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")

    ext = os.path.splitext(model_path)[1].lower()
    
    # SBML formats
    if ext in (".xml", ".sbml"):
        return io.read_sbml_model(model_path)
    
    # MATLAB format
    if ext == ".mat":
        matlab_loader = getattr(io, "load_matlab_model", None)
        if callable(matlab_loader):
            return matlab_loader(model_path)
        else:
            raise RuntimeError(
                "MATLAB format not supported. Convert to SBML or update COBRApy."
            )
    
    # JSON format
    if ext == ".json":
        try:
            return io.load_json_model(model_path)
        except Exception:
            raise RuntimeError("Could not read JSON format. Convert to SBML.")

    # Fallback: try SBML
    try:
        return io.read_sbml_model(model_path)
    except Exception as e:
        raise RuntimeError(f"Unrecognized file format for {model_path}: {e}")


def build_output_path(model_path, outdir, outname=None):
    """
    Generate output file path for translated model.
    
    Args:
        model_path: Path to original model file
        outdir: Output directory
        outname: Custom output filename (optional)
        
    Returns:
        Full path to output file
    """
    os.makedirs(outdir, exist_ok=True)
    
    if outname:
        return os.path.join(outdir, outname)
    
    # Generate automatic filename
    original_name = os.path.splitext(os.path.basename(model_path))[0]
    output_filename = f"{original_name}_translated.xml"
    return os.path.join(outdir, output_filename)


def translate_model(model_path, target_namespace, outdir, outname=None, 
                   force_overwrite=False, verbose=True):
    """
    Translate model to target namespace and save result.
    
    Args:
        model_path: Path to input model file
        target_namespace: Target namespace (e.g., 'bigg', 'seed', 'kegg')
        outdir: Output directory
        outname: Custom output filename (optional)
        force_overwrite: Overwrite existing files
        verbose: Print progress messages
        
    Returns:
        Path to created output file
        
    Raises:
        FileExistsError: If output exists and force_overwrite is False
        RuntimeError: If translation fails
    """
    if verbose:
        print(f"[INFO] Loading model: {model_path}")

    # Load model
    model = load_model(model_path)

    if verbose:
        print(f"[INFO] Translating to {target_namespace} namespace")

    # Perform translation
    try:
        translated_model = translate(model, trans_to_db=target_namespace)
    except TypeError:
        # Fallback for different mergem versions
        translated_model = translate(model, to=target_namespace)

    if translated_model is None:
        raise RuntimeError("Translation failed: mergem returned None. "
                          "Check input model or mergem compatibility.")

    # Generate output path
    output_path = build_output_path(model_path, outdir, outname)

    if os.path.exists(output_path) and not force_overwrite:
        raise FileExistsError(f"Output already exists: {output_path}. Use --force to overwrite.")

    # Save translated model
    write_sbml_model(translated_model, output_path)

    if not os.path.exists(output_path):
        raise RuntimeError("Failed to save translated model: output file not created.")

    if verbose:
        print(f"[OK] Model translated and saved to: {output_path}")

    return output_path


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Translate GEM model to different namespace (BiGG/SEED/KEGG)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--config", type=str, default="config.yaml",
                       help="Configuration YAML file")
    parser.add_argument("--model", "--models", nargs="+", dest="models",
                       help="Model file path(s) - supports multiple models (overrides config)")
    parser.add_argument("--to", type=str,
                       help="Target namespace (e.g., bigg, seed, kegg)")
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

    # Get translation section from config
    prep_config = config.get("preparation", {})
    trans_config = prep_config.get("translation") or {}

    # Priority: CLI > config
    model_paths = args.models if args.models else prep_config.get("models") or prep_config.get("model_path")
    
    # Convert single model to list for uniform processing
    if isinstance(model_paths, str):
        model_paths = [model_paths]
    
    target_namespace = args.to if args.to else trans_config.get("translation_db", "bigg")
    outdir = args.outdir if args.outdir else trans_config.get("output_dir", "outputs/translated_models")
    outname = args.outname if args.outname else trans_config.get("output_name")
    verbose = args.verbose or trans_config.get("verbose", False)
    force_overwrite = args.force or trans_config.get("force", False)

    if not model_paths:
        parser.error("No models specified. Use --model or define in config.")

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
            
            output_path = translate_model(
                model_path=model_path,
                target_namespace=target_namespace,
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
        print(f"\n[INFO] Translation completed for {len(output_paths)}/{len(model_paths)} model(s)")
    
    return 0 if len(output_paths) == len(model_paths) else 1


if __name__ == "__main__":
    exit(main())