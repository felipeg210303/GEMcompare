#You must have memote installed to run this script. 
#You can install it using pip install memote in the terminal.

from memote.suite.api import test_model
import cobra

model_path = "C:/Users/felip/Desktop/python/TESIS/iBsu1103.xml" #Path to your model file

model = cobra.io.read_sbml_model(model_path) #Charge the model

report = test_model(model) #Run the validation process with memotes function test_model

print("Results summary:") #Print the results summary
print(report.scores)
