#You can also validate models online at: https://sbml.bioquant.uni-heidelberg.de/validator_servlet/
#I recommend to validate online rather than locally, as it is faster and more visual.

import libsbml

model_path = "C:\\Users\\felip\\Desktop\\python\\B.-subtilis-FBA\\Initial models\\iYO844\\iYO844.xml" #Path to your model file

reader = libsbml.SBMLReader() #Create a new reader object
document = reader.readSBML(model_path) # Read the SBML file

if document.getNumErrors() > 0: #Verifying errors
    print("Errors found in the SBML document:")
    document.printErrors()
else:
    print("Correctly read SBML file.")

num_errors = document.checkConsistency() # Check the consistency of the model 

if num_errors > 0: #Show the errors found in the model
    print(f"{num_errors} errors found in the model:")
    for i in range(num_errors):
        error = document.getError(i)
        print(f"[{i+1}] Level {error.getSeverity()} - {error.getMessage()}")
else:
    print("Model is consistent")
