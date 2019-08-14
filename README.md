# chemotaxis_simulation

## Growth expansion model.

This python code can be used to numerically solve the Growth Expansion Model to investigate the growth expansion dynamics of a population of chemotactic cells. Code can also used to study competition and selection within the 
expanding documentions.
The Growth-Expansion model and its biological context is introduced our manuscript:
- J.Cremer, T.Honda, Y.Tang, J. Wong-Ng, M.Vergassola, T.Hwa. "Chemotaxis as a navigation strategy to thrive in nutrient-replete environments

Code can also be used to study competition and selection dynamics, as in our manuscript:
- W.Liu, J.Cremer, D.Li, T.Hwa, C.Liu. "An evolutionary stable strategy to colonize  spatially extended habitats"

August 2019, Jonas Cremer and all coauthors.

## Required modules 

Partial differential equations are numerically solved using the python module FiPy. 

- https://www.ctcms.nist.gov/fipy/
- https://github.com/usnistgov/fipy

Instructions how to install FiPy are provided on these websites.

Python 2.7 was used to run and test the code. Additional modules required include: NumPy, jason, os. Major parts of this code are provided in the file "gp_chemotaxis.py".

## Setup the parameters
Parameters can be provided as a JSON parameter file. An example parameter file, "examplesim.par", is provided. OD is used as unit for cellular densities, SI units are used for all other basic quantities. To model the competition of different strains, parameters describing cellular properties are provided two times, e.g. "lambda_1" and "lambda_2" describe the growth rate of strains 1 and 2. To vary the details of the dynamics, different dynamics can be set. Detailed equations for the different modes are defined in the major code file "gp_chemotaxis.py".

## Run simulations
To run simulations the file "chemotaxis_simulation.py" can be modified. This file shows an example how to run parameter scans, independently varying two parameters. Simulations of different parameter combinations can be run inparallel jobs and results can subsequently be merged with the merging function as provided in the script. 
