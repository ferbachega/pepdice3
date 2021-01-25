
#from pepBabel.PDBfiles import *
from MolecularSystem.molecule import *
from MolecularSystem.Fragment import get_fragment

from pepCore.Geometry import *
from pepCore.MonteCarlo import *

from pepBabel.XYZFiles import *
from pprint import pprint
import random as random


import json

#system     =  Molecule(pdb = '/home/fernando/programs/pepdice3/examples/1zdd.pdb')
system     =  Molecule(pdb = 'examples/1zdd_extended_structure_CB.pdb')







"""
construindo um biblioteca de fragmentos com base na estrutura de 1zdd

FNMQCQRRFYEALHDPNLNEEQRNAKIKSIRDDCX
XXXXX------------------------------
YYYYY------------------------------
-XXXXX-----------------------------
--XXXXX----------------------------
---XXXXX---------------------------
----XXXXX--------------------------


fragmentos = [
			 
			 [] # FRAGMENTOS DA POSICAO 0 - 4
			 
			 [dict, dict, dict] ,
		 
			 []                 ,
			 []                 ,
			 []                 ,
			 []                 ,

             ]

.
.
.
"""

# -----------------   importing fragments  ----------------------
fragments = json.load(open('examples/fragments_1zdd_example.json'))
pprint(fragments)
system.fragments = fragments
# ---------------------------------------------------------------
pprint(system.fragments)


# ----------------   importing Contacts  ------------------------
system.import_contact_map_from_file (filein  = 'examples/1zdd_pairProteinNative_example' )

# checking contacts 
for contact in system.contacts:
	print (contact.residue1.name, contact.residue1.id ,contact.residue2.name, contact.residue2.id )
# ---------------------------------------------------------------



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
				   trajectory          = 'MonteCarlo_trajectory.xyz',
				   pn                  = 1                      ,
				   log                 = True                   ,)
#'''
