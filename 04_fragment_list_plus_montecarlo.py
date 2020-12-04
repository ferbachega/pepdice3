
#from pepBabel.PDBfiles import *
from MolecularSystem.Fragment import *
from MolecularSystem.molecule import *
from pepCore.Geometry import *
from pepBabel.XYZFiles import *

from pepCore.MonteCarlo import monte_carlo_cycle

#system =  Molecule(pdb = '/home/fernando/Downloads/predicted_structure.pdb')
system     =  Molecule(pdb = '/home/fernando/programs/pepdice3/examples/1zdd.pdb')



#pdbs =  ['/home/fernando/programs/pepdice3/examples/1zdd.pdb']
#system  = build_fragment_library_from_pdbs (
#                                            molecule             = system ,
#                                            frag_size            = 7      ,
#                                            number_of_fragments  = 100    ,
#                                            pdblist              = pdbs   ,
#                                            )




for residue in system.residues:
    residue.set_phi  (angle = 180)
    residue.set_psi  (angle = 180)
    residue.set_omega(angle = 180)
    
save_XYZ_to_file(system , 'xyz_out_1zdd.xyz')


pprint(system.fragments)





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
				   cycle_size          = 10000                     ,
				   #number_of_cycles    = 10                     ,
				   
				   log_frequence       = 10                     ,
				   trajectory          = 'MonteCarlo_trajectory',
				   pn                  = 1                      ,
				   log                 = True                   ,)
