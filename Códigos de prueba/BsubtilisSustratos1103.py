# multiproceso_cobrapy.py
# -*- coding: utf-8 -*-
import os
import cobra
import pandas as pd
from typing import List, Tuple, Dict, Any
from concurrent.futures import ProcessPoolExecutor, as_completed


def load_model_flexible(path: str) -> cobra.Model:
    """
    Intenta cargar un modelo con diferentes readers según extensión.
    Soporta: .xml/.sbml (SBML) y .json (Cobra JSON).
    Si necesitas .mat u otro formato, añade el loader correspondiente.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in [".xml", ".sbml"]:
        return cobra.io.read_sbml_model(path)
    elif ext in [".json"]:
        # cobra tiene load_json_model en algunas versiones; si falla, probar read_json_model
        try:
            return cobra.io.load_json_model(path)   # preferible si existe
        except Exception:
            return cobra.io.load_json_model(path)   # fallback (si tu versión lo implementa así)
    else:
        raise ValueError(f"Extensión '{ext}' no soportada por este loader. "
                         "Añade soporte para .mat u otros formatos si lo necesitas.")

# -------------------------
# Fijar medio y condiciones
# -------------------------
def set_bounds_if_present(m: cobra.Model, rxn_id: str, lb: float, ub: float):
    try:
        rxn = m.reactions.get_by_id(rxn_id)
        rxn.bounds = (lb, ub)
    except KeyError:
        # No detener ejecución, pero avisar para debug
        # (puedes comentar el print si no quieres mensajes)
        print(f"[WARN] reacción '{rxn_id}' no encontrada en el modelo '{m.id}'.")

def preparar_medio_base(m: cobra.Model):
    # Cerrar uptake en todos los exchanges (por seguridad)
    for ex in m.exchanges:
        ex.lower_bound = 0.0

    # Abrir las compuertas del medio base según tu tabla
    set_bounds_if_present(m, "EX_h2o_e",  -1000.0, 1000.0)  # H2O
    set_bounds_if_present(m, "EX_o2_e",      -18.0,    0.0)  # O2
    set_bounds_if_present(m, "EX_pi_e",       -5.0,  1000.0)  # Orthophosphate
    set_bounds_if_present(m, "EX_co2_e",   -1000.0, 1000.0)  # CO2
    set_bounds_if_present(m, "EX_nh4_e",      -5.0,  1000.0)  # NH3
    set_bounds_if_present(m, "EX_so4_e",      -5.0,  1000.0)  # Sulfate
    set_bounds_if_present(m, "EX_ca2_e",   -1000.0, 1000.0)  # Calcium
    set_bounds_if_present(m, "EX_h_e",     -1000.0, 1000.0)  # H+
    set_bounds_if_present(m, "EX_k_e",     -1000.0, 1000.0)  # Potassium
    set_bounds_if_present(m, "EX_mg2_e",   -1000.0, 1000.0)  # Magnesium
    set_bounds_if_present(m, "EX_na1_e",   -1000.0, 1000.0)  # Sodium
    set_bounds_if_present(m, "EX_fe3_e",   -1000.0, 1000.0)  # Fe3+

    # Apagar glucosa base para que la condición la ponga explícitamente
    set_bounds_if_present(m, "EX_glc__D_e", 0.0, 1000.0)

def aplicar_fuentes_carbono(m: cobra.Model, pares_rxn_uptake: List[Tuple[str, float]]):
    for rxn_id, uptake in pares_rxn_uptake:
        set_bounds_if_present(m, rxn_id, -float(uptake), 1000.0)

# -------------------------
# Condiciones (tal como definiste)
# -------------------------
CONDICIONES = {
    "Glucosa":                    [("EX_glc__D_e", 7.63)],
    "Gluconato":                  [("EX_glcn__D_e",   5.13)],
    "Glicerol":                   [("EX_glyc_e",   6.22)],
    "Malato":                     [("EX_mal__L_e",26.51)],
    "Malato; Glucosa":            [("EX_mal__L_e",14.6), ("EX_glc__D_e", 5.95)],
    "Piruvato":                   [("EX_pyr_e",    8.26)],
    "Succinato; L-Glutamato":     [("EX_succ_e",   3.35), ("EX_glu__L_e", 2.21)],
    "Fructosa":                   [("EX_fru_e",    5.72)],
}

# -------------------------
# Función por modelo
# -------------------------
def run_for_model(path: str, condiciones: Dict[str, List[Tuple[str, float]]]) -> pd.DataFrame:
    """
    Carga un modelo desde 'path', corre todas las condiciones y
    devuelve un DataFrame con columnas: model_path, model_id, condition, status, biomass.
    """
    model = load_model_flexible(path)
    # Opcional: asignar nombre legible
    model.id = os.path.basename(path)
    resultados = []
    for cond_name, pares in condiciones.items():
        mcopy = model.copy()
        preparar_medio_base(mcopy)
        aplicar_fuentes_carbono(mcopy, pares)
        sol = mcopy.optimize()
        resultados.append({
            "model_path": path,
            "model_id": model.id,
            "condition": cond_name,
            "status": sol.status,
            "biomass": sol.objective_value
        })
    return pd.DataFrame(resultados)

# -------------------------
# Runner maestro
# -------------------------
def run_multiple_models(model_paths: List[str],
                        condiciones: Dict[str, List[Tuple[str, float]]] = CONDICIONES,
                        parallel: bool = False,
                        max_workers: int = None) -> pd.DataFrame:
    """
    Corre run_for_model para cada ruta en model_paths.
    Si parallel=True intentará usar multiprocessing (ProcessPoolExecutor).
    Devuelve un DataFrame combinado.
    """
    dfs = []
    if not parallel:
        # Modo secuencial
        for p in model_paths:
            print(f"[INFO] Procesando (secuencial) {p} ...")
            df = run_for_model(p, condiciones)
            dfs.append(df)
    else:
        # Modo paralelo (cada proceso leerá su modelo desde disco)
        print(f"[INFO] Procesando en paralelo {len(model_paths)} modelos ...")
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            futures = {ex.submit(run_for_model, p, condiciones): p for p in model_paths}
            for fut in as_completed(futures):
                p = futures[fut]
                try:
                    df = fut.result()
                    dfs.append(df)
                    print(f"[INFO] Terminado {p}")
                except Exception as e:
                    print(f"[ERROR] fallo en {p} -> {e}")

    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
    else:
        combined = pd.DataFrame(columns=["model_path", "model_id", "condition", "status", "biomass"])
    return combined

if __name__ == "__main__":
    # Lista rutas aquí:
    model_paths = [
        "B.-subtilis-FBA\Models\iYO844.xml",
        "B.-subtilis-FBA\Models\iBsu1103_bigg.xml",
        "B.-subtilis-FBA\Models\iBsu1103v2_bigg.xml",
        "B.-subtilis-FBA\Models\iBsu1147.xml",
        "B.-subtilis-FBA\Models\ec_iYO844_bigg.xml",
        "B.-subtilis-FBA\Models\iBB1018.xml",
        "B.-subtilis-FBA\Models\ecBSU1.xml",
        "B.-subtilis-FBA\Models\iBsu1147R.xml",
        "B.-subtilis-FBA\Models\pan_model_bigg.xml",
       
    ]

    resultados_df = run_multiple_models(model_paths, parallel=False)
    print(resultados_df)

    # Guardar resultado combinado
    resultados_df.to_csv("resultados_biomasa_por_modelo.csv", index=False)
