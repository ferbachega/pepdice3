
#from pepBabel.PDBfiles import *
from MolecularSystem.molecule import *
from pepCore.Geometry import *
from pepBabel.XYZFiles import *
from pprint import pprint


system     =  Molecule(pdb = '/home/fernando/programs/pepdice3/examples/1zdd.pdb')





def get_fragment (molecule = None, size = 5, position = 0):
	""" Function doc 
	"""

	fragment  = []
	for residue in molecule.residues[position:position+size]:
		500,00
		name = residue.name
		phi  = residue.get_phi()
		psi  = residue.get_psi()
		omega= residue.get_omega()
		inex = residue.id

		residue_data =  {
		                 'index' : inex  ,
		                 'name'  : name  ,
		                 'phi'   : phi   , 
		                 'psi'   : psi   ,
		                 'omega' : omega ,
		                 } 
		
		fragment.append(residue_data)

	return fragment






# taking only one fragment 
fragment = get_fragment(molecule = system, size = 5, position = 0 )
pprint(fragment)






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
	fragment = get_fragment(molecule = system, size = 5, position = i )
	fragments.append(fragment)

pprint(fragments[1])



