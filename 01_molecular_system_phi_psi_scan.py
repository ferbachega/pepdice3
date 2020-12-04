
#from pepBabel.PDBfiles import *
from MolecularSystem.molecule import *
from pepCore.Geometry import *
from pepBabel.XYZFiles import *




system     =  Molecule(pdb = '/home/fernando/programs/pepdice3/examples/poliAla.pdb')
#trajout    =  '/home/fernando/programs/pepdice3/examples/poliAla_traj_from_restraint'


psi = system.residues[5].get_psi()
system.residues[5].psi_restraint_angle  = psi
system.residues[5].psi_restraint_weight = 10.0

for i in range (1, 360, 5):
	system.residues[5].set_psi(psi+i)
	print(system.residues[5].psi ,psi+i, system.energy())
	#save_XYZ_to_file(system , trajout+'.xyz')


phi = system.residues[5].get_phi()
system.residues[5].phi_restraint_angle  = phi
system.residues[5].phi_restraint_weight = 10.0

for i in range (1, 360, 5):
	system.residues[5].set_phi(phi+i)
	print(system.residues[5].phi ,phi+i, system.energy())
	#save_XYZ_to_file(system , trajout+'.xyz')
	
	
	
	
	
	
	
	
	
	
	
	
