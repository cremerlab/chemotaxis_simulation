###################################################################################
#Chemotaxis simulation
###################################################################################
#Code to simulate chemotaxis driven expansion and population growth,
#following the Growth-Expansion Model.
#Code can also used to study competition and selection within the 
#expanding documentions.
#The Growth-Expansion model and its biological context is introduced
#our manuscript:
#J.Cremer, T.Honda, Y.Tang, J. Wong-Ng, M.Vergassola, T.Hwa 
#"Chemotaxis as a navigation strategy to thrive in nutrient-replete environments"
#
#Code can also be used to study competition and selection dynamics, as
#in our manuscript:
#W.Liu, J.Cremer, D.Li, T.Hwa, C.Liu
#"An evolutionary stable strategy to colonize  spatially extended habitats"
#
#August 2019, Jonas Cremer and all coauthors.
#
###################################################################################
#Simulations using Python 2.7 and the partial differential equation solver FiPy
#Code available via GitHub at. See provided README file for additional information.
###################################################################################

from fipy import *
from datetime import datetime
import numpy as np
import json
import time
import os, sys

###################################################################################
#Import major part of the simulation code
###################################################################################
import gp_chemotaxis

###################################################################################
#set name of simulation
###################################################################################
name_sim='examplesim' #make sure parameter file with same name is present within (e.g. 'examplesim.par')

###################################################################################
#load parameter file
###################################################################################
parameters=gp_chemotaxis.load_parameters_simulation(name_sim+'.par')

###################################################################################
#run parameter scan, varying two parameters
###################################################################################          
#set which variables should be changed (variables have to exist in parameter file)
#1st parameter varied
xpar="xhi_1"
x_values=np.array([2.0,4.0])*np.power(10.,-9.)
#2nd parameter varied
ypar="D_nutrients"
y_values=np.array([10])*np.power(10.,-12.) 

#run parameter scan
gp_chemotaxis.run_sweep(name_sim,parameters,x_values,xpar,y_values,ypar,scanparameter=True,nametag="1")
#note: if scanparameter=True, data is not overwritten
#note: to simulate selection dynamics over several iterations, selection of previous iteration can be read in as initial information.

###################################################################################
###merge different realizations into one file
###run this to combine simulations into one major output file if parameter scans are run independently (for example within a HPCC environment)
#####################################################################################
gp_chemotaxis.merge_simulation(name_sim,diagonalonly=-1,clusterrun=False,changetodifference=False)
