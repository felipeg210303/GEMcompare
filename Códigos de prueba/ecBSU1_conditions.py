import cobra
# Cargar los mimport cobra

# ================================
# Cargar tu modelo
# ================================
# Ojo: cambia la ruta a tu modelo real (.xml/.json/.mat)
model = cobra.io.read_sbml_model("tu_modelo.xml")

# ================================
# Configurar el medio
# ================================

# 1. Primero apagamos todas las fuentes de carbono que mencionaste
carbon_sources = [
    "EX_glc__D_e",   # Glucosa
    "EX_glcn_e",     # Gluconato
    "EX_glyc_e",     # Glicerol
    "EX_mal__L_e",   # Malato
    "EX_pyr_e",      # Piruvato
    "EX_succ_e",     # Succinato
    "EX_glu__L_e",   # L-Glutamato
    "EX_fru_e"       # Fructosa
]

for rxn in carbon_sources:
    if rxn in model.reactions:
        model.reactions.get_by_id(rxn).bounds = (0, 1000)

# 2. Ahora añadimos la glucosa como única fuente de carbono
# El flujo que mencionaste fue -8.7 → lo pongo como límite inferior
if "EX_glc__D_e" in model.reactions:
    model.reactions.get_by_id("EX_glc__D_e").bounds = (-8.7, 1000)

# 3. Ajustamos el resto del medio según lo que diste
# Oxígeno
model.reactions.get_by_id("EX_o2_e").bounds = (-18, 0)
# Fosfato
model.reactions.get_by_id("EX_pi_e").bounds = (-5, 1000)
# Amonio
model.reactions.get_by_id("EX_nh4_e").bounds = (-5, 1000)
# Sulfato
model.reactions.get_by_id("EX_so4_e").bounds = (-5, 1000)

# Los demás (H2O, CO2, H+, etc.) ya suelen estar abiertos en BiGG,
# pero si quieres ser estricto puedes setearlos:
model.reactions.get_by_id("EX_h2o_e").bounds = (-1000, 1000)
model.reactions.get_by_id("EX_co2_e").bounds = (-1000, 1000)
model.reactions.get_by_id("EX_ca2_e").bounds = (-1000, 1000)
model.reactions.get_by_id("EX_h_e").bounds = (-1000, 1000)
model.reactions.get_by_id("EX_k_e").bounds = (-1000, 1000)
model.reactions.get_by_id("EX_mg2_e").bounds = (-1000, 1000)
model.reactions.get_by_id("EX_na1_e").bounds = (-1000, 1000)
model.reactions.get_by_id("EX_fe3_e").bounds = (-1000, 1000)

# ================================
# Optimizar (FBA)
# ================================
solution = model.optimize()

print("Valor de la función objetivo (biomasa):", solution.objective_value)
print("Principales flujos:")
print(solution.fluxes.sort_values(ascending=False).head(15))

model1 = cobra.io.read_sbml_model("ruta/al/modelo.xml")