from cobra.io import read_sbml_model
from cobra.flux_analysis import single_gene_deletion
import pandas as pd

if __name__ == '__main__':
    # Cargar el modelo
    cobra_model = read_sbml_model("C:\\Users\\felip\\Desktop\\python\\B.-subtilis-FBA\\Models\\iYO844.xml")
    
    # Configurar el medio de cultivo específico
    print("Configurando medio de cultivo...")
    
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
        if rxn in cobra_model.reactions:
            cobra_model.reactions.get_by_id(rxn).bounds = (0, 1000)

    # 1. Glucosa
    if "EX_glc__D_e" in cobra_model.reactions:
        cobra_model.reactions.get_by_id("EX_glc__D_e").bounds = (-8.7, 1000)
    
    # 2. Oxígeno
    cobra_model.reactions.get_by_id("EX_o2_e").bounds = (-18, 0)
    
    # 3. Nutrientes esenciales
    cobra_model.reactions.get_by_id("EX_pi_e").bounds = (-5, 1000)    # Fosfato
    cobra_model.reactions.get_by_id("EX_nh4_e").bounds = (-5, 1000)   # Amonio
    cobra_model.reactions.get_by_id("EX_so4_e").bounds = (-5, 1000)   # Sulfato
    
    # 4. Otros componentes (abiertos)
    cobra_model.reactions.get_by_id("EX_h2o_e").bounds = (-1000, 1000)
    cobra_model.reactions.get_by_id("EX_co2_e").bounds = (-1000, 1000)
    cobra_model.reactions.get_by_id("EX_ca2_e").bounds = (-1000, 1000)
    cobra_model.reactions.get_by_id("EX_h_e").bounds = (-1000, 1000)
    cobra_model.reactions.get_by_id("EX_k_e").bounds = (-1000, 1000)
    cobra_model.reactions.get_by_id("EX_mg2_e").bounds = (-1000, 1000)
    cobra_model.reactions.get_by_id("EX_na1_e").bounds = (-1000, 1000)
    cobra_model.reactions.get_by_id("EX_fe3_e").bounds = (-1000, 1000)
    
    # Verificar el estado del modelo con el nuevo medio
    solution = cobra_model.optimize()
    print('Valor objetivo con el medio especificado:', solution.objective_value)
    
    # Realizar deleción génica
    print("Realizando deleción génica...")
    deletion_results = single_gene_deletion(cobra_model)
    
    # Filtrar resultados: crecimiento entre 0.1 y 0.6
    filtered_results = deletion_results[
        (deletion_results['growth'] >= solution.objective_value*0.2) & 
        (deletion_results['growth'] <= solution.objective_value*0.8)
    ].copy()
    
    # Calcular porcentaje de crecimiento
    filtered_results['growth_percentage'] = (filtered_results['growth'] / solution.objective_value * 100).round(2)
    
    # Mostrar resultados
    print(f"\nGenes con crecimiento entre 0.2X y 0.8X después de deleción: {len(filtered_results)}")
    print(filtered_results[['growth', 'growth_percentage', 'status']])
    
    # Exportar a CSV
    if len(filtered_results) > 0:
        output_path = "C:\\Users\\felip\\Desktop\\python\\B.-subtilis-FBA\\Results\\knockout_results_filtered2080.csv"
        filtered_results.to_csv(output_path, index=True, encoding='utf-8-sig')
        print(f"\nResultados exportados a: {output_path}")
        
    else:
        print("No se encontraron genes que cumplan los criterios.")
    
    # Información adicional
    print(f"\nInformación del medio de cultivo aplicado:")
    print(f"- Glucosa: {-cobra_model.reactions.EX_glc__D_e.lower_bound} mmol/gDW/h")
    print(f"- Oxígeno: {-cobra_model.reactions.EX_o2_e.lower_bound} mmol/gDW/h")
    print(f"- Fosfato: {-cobra_model.reactions.EX_pi_e.lower_bound} mmol/gDW/h")
    print(f"- Amonio: {-cobra_model.reactions.EX_nh4_e.lower_bound} mmol/gDW/h")
    print(f"- Sulfato: {-cobra_model.reactions.EX_so4_e.lower_bound} mmol/gDW/h")