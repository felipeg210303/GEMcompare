# -*- coding: utf-8 -*-
import cobra
import pandas as pd

# ================================
# 1) Cargar tu modelo
# ================================
# Cambia la ruta/archivo a tu modelo real (.xml/.json/.mat)
model = cobra.io.read_sbml_model("B.-subtilis-FBA\Models\ecBSU1_Converteddd.xml")

# ================================
# 2) Utilidades
# ================================
def set_bounds_if_present(m, rxn_id, lb, ub):
    """Asigna bounds si la reacción existe; si no, lo reporta."""
    try:
        rxn = m.reactions.get_by_id(rxn_id)
        rxn.bounds = (lb, ub)
    except KeyError:
        print(f"[ADVERTENCIA] No existe la reacción '{rxn_id}' en el modelo.")

def preparar_medio_base(m):
    """
    Cierra uptake en TODOS los exchanges y configura el medio base
    exactamente como en tu tabla, pero con D-Glucosa en 0 para
    'reemplazarla' por la(s) fuente(s) que probaremos en cada condición.
    """
    # 2.1 Cerrar uptake por defecto
    for ex in m.exchanges:
        ex.lower_bound = 0.0  # sin entrada a menos que la abramos explícitamente

    # 2.2 Abrir lo que definiste en el medio
    set_bounds_if_present(m, "EX_h2o_e",  -1000.0, 1000.0)  # H2O
    set_bounds_if_present(m, "EX_o2_e",      -10.0,    0.0)  # O2
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

    # Importante: apagar glucosa del medio base para "reemplazarla"
    set_bounds_if_present(m, "EX_glc__D_e", 0.0, 1000.0)     # D-Glucose

def aplicar_fuentes_carbono(m, pares_rxn_uptake):
    """
    Activa la(s) fuente(s) de carbono de una condición dada.
    'pares_rxn_uptake' es lista de (rxn_id_BiGG, uptake_positivo),
    y aquí se aplica como lower_bound = -uptake (convención FBA).
    """
    for rxn_id, uptake_pos in pares_rxn_uptake:
        set_bounds_if_present(m, rxn_id, -float(uptake_pos), 1000.0)

# ================================
# 3) Tabla de condiciones EXACTAS que pediste
#    (IDs BiGG entre paréntesis)
# ================================
condiciones = {
    "Glucosa":                    [("EX_glc__D_e", 7.63)],             # D-Glucose
    "Gluconato":                  [("EX_glcn_e",   5.13)],             # Gluconate
    "Glicerol":                   [("EX_glyc_e",   6.22)],             # Glycerol
    "Malato":                     [("EX_mal__L_e",26.51)],             # L-Malate
    "Malato; Glucosa":            [("EX_mal__L_e",14.6),
                                   ("EX_glc__D_e", 5.95)],
    "Piruvato":                   [("EX_pyr_e",    8.26)],             # Pyruvate
    "Succinato; L-Glutamato":     [("EX_succ_e",   3.35),              # Succinate
                                   ("EX_glu__L_e", 2.21)],             # L-Glutamate
    "Fructosa":                   [("EX_fru_e",    5.72)],             # D-Fructose
}

# ================================
# 4) Ejecutar FBA por condición
# ================================
resultados = []
for nombre, pares in condiciones.items():
    m = model.copy()
    preparar_medio_base(m)            # medio como el adjunto, sin glucosa base
    aplicar_fuentes_carbono(m, pares) # activar la(s) fuente(s) de la condición
    sol = m.optimize()
    resultados.append({
        "condicion": nombre,
        "status": sol.status,
        "biomasa": sol.objective_value
    })

df = pd.DataFrame(resultados).sort_values("biomasa", ascending=False)
print(df.to_string(index=False))

# (Opcional) inspeccionar flujos de una condición concreta:
# m_dbg = model.copy()
# preparar_medio_base(m_dbg)
# aplicar_fuentes_carbono(m_dbg, condiciones["Malato; Glucosa"])
# sol_dbg = m_dbg.optimize()
# print(sol_dbg.fluxes.sort_values(ascending=False).head(20))
