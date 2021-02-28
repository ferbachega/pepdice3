
#from pepBabel.PDBfiles import *
from MolecularSystem.molecule import *
from pepCore.Geometry import *
from pepBabel.XYZFiles import *



system     =  Molecule(pdb = '/home/fernando/programs/pepdice3/examples/1zdd.pdb')


system.energy()
print (system.energy_components)

#print (system.energy())

#phi_psi_list = system.get_phi_psi_list()
#
#for residue in system.residues:
#	print ("residue:{:>}   phi: {}    psi: {}".format(residue.id, phi_psi_list[residue.id][0] , phi_psi_list[residue.id][1]) )
