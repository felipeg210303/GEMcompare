# You can also translate to other formats such as BiGG, KEGG, or MetaCyc. 

from mergem import translate
import cobra
from cobra.io import write_sbml_model
import os

model_path = "C:\\Users\\felip\\Desktop\\python\\B.-subtilis-FBA\\Pan-genomic model\\modelo_salida.xml" #Path to your model file
output_path = "C:\\Users\\felip\\Desktop\\python\\B.-subtilis-FBA\\Translated Models (to BiGG)\\mergem_pan_model_bigg.xml" #Path where you want to save the translated model file

model = cobra.io.read_sbml_model(model_path) #Charge the model

translated_model = translate(model, trans_to_db="bigg") #Translate the model to Bigg format using the translate function

write_sbml_model(translated_model, output_path) #Save the translated model to a file
#translated_model is a cobra model object so write_sbml_model can be used to save it to a file

if os.path.exists(output_path): #Verifying correct creation of the file
    print(f" Model translated and saved to the following path:\n{output_path}")
else:
    print(" Error ocurred while saving the translated model")
