#You must have memote installed to run this script. 
#You can install it using pip install memote in the terminal.

from memote.suite.api import test_model
import cobra

model_path = r"B.-subtilis-FBA\Models\iYO844.xml" #Path to your model file

model = cobra.io.read_sbml_model(model_path) #Read the model

print("Loading... (this might take some time)") 
report = test_model(model) #Run the validation process with memotes function test_model

print("Results summary:") #Print the results summary
print(report.scores)

# You can also run and save a memote report in your terminal with: memote report snapshot --filename "YOUR_FILENAME" [PATH TO YOUR MODEL]