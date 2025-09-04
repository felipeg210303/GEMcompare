import cobra

# Cargar el modelo metabólico (
model = cobra.io.read_sbml_model("B.-subtilis-FBA\\iYO844\\iYO844_bigg.xml")

# Flujos de condiciones experimentales
model.reactions.get_by_id("ex_cellb_e").lower_bound = -7.71  # Consumo de glucosa
model.reactions.get_by_id("o2tu").lower_bound = -18  # Consumo de O2
    
# Optimizar el modelo usando FBA

solution = model.optimize()
print(f"Objective value: {solution.objective_value}")
print(f"Acetate production: {solution.fluxes['ACts']}")

