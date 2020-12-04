
#from pepBabel.PDBfiles import *
from MolecularSystem.molecule import *
from pepCore.Geometry import *
from pepBabel.XYZFiles import *




system     =  Molecule(pdb = '/home/fernando/programs/pepdice3/examples/poliAla.pdb')
trajout    =  '/home/fernando/programs/pepdice3/examples/poliAla_rotation_traj'


psi = system.residues[5].get_psi()


for i in range (1, 360, 5):
	system.residues[5].set_psi(psi+i)
	save_XYZ_to_file(system , trajout+'.xyz')


phi = system.residues[5].get_phi()
for i in range (1, 360, 5):
	system.residues[5].set_phi(phi+i)
	save_XYZ_to_file(system , trajout+'.xyz')

