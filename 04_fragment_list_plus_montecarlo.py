
#from pepBabel.PDBfiles import *
from MolecularSystem.molecule import *
from MolecularSystem.Fragment import get_fragment

from pepCore.Geometry import *
from pepCore.MonteCarlo import *

from pepBabel.XYZFiles import *
from pprint import pprint
import random as random

#system     =  Molecule(pdb = '/home/fernando/programs/pepdice3/examples/1zdd.pdb')
system     =  Molecule(pdb = '/home/fernando/programs/pepdice3/examples/1zdd_no_side.pdb')







"""
construindo um biblioteca de fragmentos com base na estrutura de 1zdd

FNMQCQRRFYEALHDPNLNEEQRNAKIKSIRDDCX
XXXXX------------------------------
-XXXXX-----------------------------
--XXXXX----------------------------
---XXXXX---------------------------
----XXXXX--------------------------
.
.
.
"""

fragments = []
for i in range (0, len(system.residues) -5):
	fragments.append([])

for i in range (0, len(system.residues) -5):
	fragment = get_fragment(molecule = system, size = 5, position = i )
	fragments[i].append(fragment)

	
	
	#variant = {}
	#for index in fragment.keys():
	#	
	#	phi   = fragment[index]["PHI"]
	#	if  phi:
	#		phi   =+ random.random()
	#	
	#	psi   = fragment[index]["PSI"]
	#	if psi:
	#		psi   = random.random()
	#	
	#	
	#	omega = fragment[index]["OMEGA"]
	#	name  = fragment[index]["NAME"]
	#	variant[index] ={
	#	                 'NAME'  : name  ,
	#	                 'PHI'   : phi   , 
	#	                 'PSI'   : psi   ,
	#	                 'OMEGA' : omega ,
	#	                 }
	#
	#fragments[i].append(variant)

system.fragments = fragments







for residue  in system.residues:
	residue.set_phi(180)
	residue.set_psi(180)
	residue.set_omega(180)

	#print(system.residues[5].psi ,psi+i, system.energy())










monte_carlo_cycle (molecule            = system                   ,
				   random              = None                   ,
				   
				   temperature          = 1                      ,
				   #final_temperature   = False                  ,
				   #gamma               = False                  ,
				   
				   
				   Kb                  = 1                      , # 0.0083144621               ,
				   angle_range         = 30                     ,
				   fragment_rate       = 1.0                    , #between 0  and 1
				   fragment_sidechain  = False                  ,
				   PhiPsi_rate         = 0.0                    ,
				   attempt_per_residue = 5                      ,
				   					
				   
				   #simulated_annealing = None                   , # exp, linear
				   cycle_size          = 1000                     ,
				   #number_of_cycles    = 10                     ,
				   
				   log_frequence       = 10                     ,
				   trajectory          = 'MonteCarlo_trajectory',
				   pn                  = 1                      ,
				   log                 = True                   ,)
