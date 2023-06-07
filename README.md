REFERENCE: _"Intermittent migration can induce pulses of speciation in a two-island system"_, by
Debora Princepe, Simone Czarnobai, Rodrigo A. Caetano,Flavia M. D. Marquitti, Marcus A.M. de Aguiar, and Sabrina B. L. Araujo

We provide one code, written in FORTRAN: "Program_intermittent.f90"

	This program has one input file, "input_sea.txt", and one output, "Temporal.txt":
	
	a) The input_sea.txt provides in the first line five parameters (the values in the parentheses refer to the values in the file):
	a.1) Number of iterations (2000) 
	a.2) Initial migration probability (0.4)
	a.3) Number of repetitions for the set of parameters (1)
	a.4) Mutation rate (0.001)
	a.5) Time interval that species identification is made (10 iterations)
	The remaining lines provide the step iteration and the value of migration. The migration variation period (written over the subsequent lines) is based on the sea level data assuming a seabed of h=-50m.

Other parameters, such as minimal genetic similarity, genome size, and population size,  can be adjusted in the main program.
	
	b) The output "Temporal.txt" has 7 columns labeled as: 
	b.1) pop - the number of total individuals inhabiting the islands (parameter)
	b.2) rep - the number of the simulation repetition (parameter)
	b.3) mig - the migration rate (parameter)
	b.4) tem - the time iteration; saved as multiple of (a.5)
	b.5) esp - the ID of the species (initially, all individuals are conal, and their species ID is 1)
	b.6) abund1 - abundance of esp on island 1 at time iteration tem
	b.7) abund2 - abundance of esp on island 2 at time iteration tem

	For additional information, please, contact Sabrina Araujo (araujosbl@fisica.ufpr.br). 
