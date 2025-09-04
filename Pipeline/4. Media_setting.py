import cobra
import pandas as pd

model_path=r"B.-subtilis-FBA\Models\ecBSU1.xml" #Path to your model
model = cobra.io.read_sbml_model(model_path) #Load the model

# Close all exchange reactions
for rxn in model.exchanges:
    rxn.bounds = (0, 0)

# Load culture medium data from Excel
# Every sheet corresponds to a different culture medium
excel_file = r"B.-subtilis-FBA\Pipeline\Culture_media.xlsx"
sheet_name = "Minimal_iYO844"   # Change this to the desired culture medium

df = pd.read_excel(excel_file, sheet_name=sheet_name)

# --- Apply the limits to the model ---
for _, row in df.iterrows():
    rxn_id = row["ID"]
    lower = float(str(row["LOWER_BOUND"]).replace(",", ".")) 
    upper = float(str(row["UPPER_BOUND"]).replace(",", "."))
    
    if rxn_id in model.reactions:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = lower
        rxn.upper_bound = upper
    else:
        print(f"⚠️ Reacción {rxn_id} no encontrada en el modelo")

# --- Test with FBA ---
solution = model.optimize()
print("Crecimiento predicho:", solution.objective_value)
