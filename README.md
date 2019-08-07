# chemotaxis_simulation

## Growth expansion model.

This python code can be used to numerically solve the growth expansion model to investigate the growth expansion dynamics of a population of chemotactic cells. The model is introduced in Cremer et al. Variations including the competition dynamics of different individuals as introduced in Liu et al. can also be analyzed with this script.

Partial differential equations are numerically solved using the python module FiPy. 
- https://www.ctcms.nist.gov/fipy/
- https://github.com/usnistgov/fipy

Instructions how to install FiPy are provided on these websites.

## Setup the parameters

Parameters are provided in a jason parameter file. An example parameter file "parameter.par" is provided. 

## Run simulations

To run simulations use chemotaxis_simulation.py. This file shows an example how to run parameter scans, independently varying two parameters. Results can subsequently be merged with the merging function as provided in the script. 

