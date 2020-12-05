
#from pepBabel.PDBfiles import *
from MolecularSystem.molecule import Molecule
from MolecularSystem.Fragment import get_fragment
from pepCore.Geometry import *
from pepBabel.XYZFiles import *
from pprint import pprint


system     =  Molecule(pdb = '/home/fernando/programs/pepdice3/examples/1zdd.pdb')


#taking just one fragment
fragment = get_fragment(molecule = system, size = 5, position = 0 )






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

pprint(fragments)


