MetaC

These files implement the simulation study of the manuscript entitled "Accounting for Time Dependency in Meta-Analyses of Concordance Probability Estimates" by Matthias Schmid, Tim Friede, Nadja Klein and Leonie Weinhold.
Comparison of meta-analysis and -regression methods for C-indices of studies with varying follow-up times. 

This folder contains R-Code to 
	(i) simulate C-indices of multiple studies , 
	(ii) to run different meta-analysis and -regression models to aggregate the C-indices of multiple studies
	(iii) to evaluate the meta-analysis and -regression models.
	
	
S0_functions.R: contains all functions needed for above described steps (i)-(ii)
S1_Simulation_generate_meta_data.R: contains function to generate meta-analysis data (study specific C-indices)
S2_Simulation_run_meta_analysis.R: contains function to  run all meta-analysis and -regression models
S3_Results.R: generates all plots and tables for result presentation
