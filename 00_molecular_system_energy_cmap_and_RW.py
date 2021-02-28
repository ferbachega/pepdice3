
#from pepBabel.PDBfiles import *
from MolecularSystem.molecule import *
from MolecularSystem.Fragment import get_fragment

from pepCore.Geometry import *
from pepCore.MonteCarlo import *

from pepBabel.XYZFiles import *
from pprint import pprint
import random as random


import json

system     =  Molecule(pdb = 'examples/1zdd_CB.pdb')

'''
# -----------------   importing fragments  ----------------------
fragments = json.load(open('examples/fragments_1zdd_example.json'))
system.fragments = fragments
# ---------------------------------------------------------------
'''

# ----------------   importing Contacts  ------------------------
system.import_contact_map_from_file (filein  = 'examples/1zdd_pairProteinNative' )
# ---------------------------------------------------------------

system.energy()
print ('n_atoms:', len(system.atoms),'n_residues:' ,len(system.residues) , system.energy_components)
