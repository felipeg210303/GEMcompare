# B. subtilis FBA

Here we are going to develop a repository of genome-scale metabolic models (GEMs) available for *Bacillus subtilis* which are reviewed and well documented. We are also going to include some use cases and examples for each model and compare the results obtained between models with the computational tool of Flux Balance Analysis (FBA) performed in python using the library "cobrapy" based in COnstraint Based Reconstruction and Analysis (COBRA) and with experimental data.

Objectives of this repository: 

	*To make genome-scale metabolic models (GEMs) of Bacillus subtilis easily accesible and comprehensible in one single place even for people with low programming skills and allow them to make hypothesis with the use of this models and Flux Balance Analysis (FBA).

	*To elaborate a in silico experimental protocol for the use of genome-scale metabolico models (GEMs) of *Bacillus subtilis*.

	*To compare some of the available GEMs of *Bacillus subtilis* using the proposed protocol as pipeline.

You can find the tutorials and documentation for the use of this repository in the following Jupyter Notebooks: <>

The associated scientfic paper can be found in the following link: <>

Important definitions and concepts:
	*FBA (Flux Balance Analysis): Predicts metabolic fluxes by optimizing an objective function (e.g., growth) under steady-state and constraint-based assumptions.
	*COBRA (COsntraint Based Reconstruction and Analysis): A framework for the analysis of metabolic networks. It is based in the use of linear programming to solve optimization problems.
	*GEM (Genome-scale metabolic model): A model of a metabolic network that includes all the reactions and metabolites of a genome.
	*GPR (Gene-Protein-Reaction): Logical links connecting genes to the reactions they enable via enzymes.
	*Objective Function: A mathematical goal (usually biomass production) that FBA tries to maximize or minimize.
	*Constraints: Biological and physical limits applied to reactions, like nutrient availability or reaction directionality.
	*Stoichiometric Matrix (S-matrix): A matrix showing how metabolites participate in reactions. It’s the foundation for FBA calculations.
	*Exchange Reactions: Represent metabolite import/export between the cell and environment. Define the "media" of the model.
	*Biomass reaction: A pseudo-reaction that represents the synthesis of all cellular components needed for growth.
	*Metabolic Flux: The rate at which metabolites flow through biochemical reactions in a network. It reflects how active a pathway is and is typically measured in mmol/gDW/h (millimoles per gram dry weight per hour).
	*FVA (Flux Variability Analysis): Determines the minimum and maximum possible fluxes for each reaction while still achieving the objective.
	*pFBA (parsimonious FBA): Finds the simplest set of active reactions that still meet the optimal objective, mimicking real biological efficiency.
	*Gene knock out: Simulates deleting genes to study their effect on metabolism or growth

In case you need any help you can contact me at: <gr.gabriel@javeriana.edu.co>