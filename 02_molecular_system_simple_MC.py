
#from pepBabel.PDBfiles import *
from MolecularSystem.molecule import *
from pepCore.Geometry import *
from pepBabel.XYZFiles import *

from pepCore.MonteCarlo import monte_carlo_cycle

system     =  Molecule(pdb = '/home/fernando/programs/pepdice3/examples/poliAla.pdb')
trajout    = '/home/fernando/programs/pepdice3/examples/poliAla_MC_helix_restraint'


for residue in system.residues:
	residue.phi_restraint_angle  = -60.0
	residue.phi_restraint_weight = 100.0
	
	residue.psi_restraint_angle  = -40.0
	residue.psi_restraint_weight = 100.0
	
energy = system.energy()
print ('energy', energy)



#'''
monte_carlo_cycle (molecule            = system                  ,
				   random              = None                    ,
				   
				   temperature          = 1                      ,
				   #final_temperature   = False                  ,
				   #gamma               = False                  ,
				   
				   
				   Kb                  = 10                      , # 0.0083144621               ,
				   angle_range         = 30                     ,
				   fragment_rate       = 0.0                    , #between 0  and 1
				   fragment_sidechain  = False                  ,
				   PhiPsi_rate         = 1.0                    ,
				   attempt_per_residue = 5                      ,
				   					
				   
				   #simulated_annealing = None                   , # exp, linear
				   cycle_size          = 500                     ,
				   #number_of_cycles    = 10                     ,
				   
				   log_frequence       = 10                     ,
				   trajectory          = trajout                ,
				   pn                  = 1                      ,
				   log                 = True                   ,)
#'''
