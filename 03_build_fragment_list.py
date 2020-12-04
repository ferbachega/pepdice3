
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



def get_fragments_ (_):
	""" Function doc """
	




fragment = get_fragment(molecule = system, )
pprint(fragment)















#phipsilist = system.get_phi_psi_list()
#e_rama     = system.get_phi_psi_energy()




#print (phipsilist[1:-1])


# testing the function that centers on the x, y and z coordinates  
#print (system.residues[0].atoms[0].pos)
#print (system.residues[0].atoms[2].pos)
#center_atom(system, 
#            system.residues[0].atoms[0].pos[0], 
#            system.residues[0].atoms[0].pos[1], 
#            system.residues[0].atoms[0].pos[2])
#print (system.residues[0].atoms[0].pos)
#print (system.residues[0].atoms[2].pos)


#print (system.residues[0].atoms[2].pos)

'''
for residue in system.residues:
	
	phi = computePhiPsi (molecule = system, 
							   resi = residue.id, 
							   bond = 'PHI')
	
	psi = computePhiPsi (molecule = system, 
							   resi =  residue.id, 
							   bond = 'PSI')
							   
	print (residue.name, residue.id, phi, psi)
#'''

'''
for residue in system.residues:
	phi = residue.get_phi()
	psi = residue.get_psi()
	ome = residue.get_omega()
	print (residue.name, residue.id, phi, psi, ome)
	
	#print(residue.atoms[0].resn,  residue.atoms[0].pos, residue.atom_dic['N'].pos,  residue.atoms[0].pos[0] - residue.atom_dic['N'].pos[0] )
#'''

save_XYZ_to_file(system , 'xyz_out_1zdd.xyz')


'''
for residue in system.residues:
	phi = residue.get_phi()
	psi = residue.get_psi()
	ome = residue.get_omega()
	try:
		residue.set_phi  (round(residue.get_phi()  ))
		residue.set_psi  (round(residue.get_psi()  ))
		residue.set_omega(round(residue.get_omega()))
	except:
		pass
	#set_phi_psi_dihedral (molecule = system, resi= residue.id , bond='PHI', angle = -60 )
	#set_phi_psi_dihedral (molecule = system, resi= residue.id , bond='PSI', angle = -40 )
#'''	
save_XYZ_to_file(system , 'xyz_out_1zdd.xyz')


'''
for residue in system.residues:
	phi = residue.get_phi()
	psi = residue.get_psi()
	ome = residue.get_omega()
	print (residue.name, residue.id, phi, psi, ome)
	
	#print(residue.atoms[0].resn,  residue.atoms[0].pos, residue.atom_dic['N'].pos,  residue.atoms[0].pos[0] - residue.atom_dic['N'].pos[0] )
#'''


	
'''
for residue in system.residues:
	phi = residue.get_phi()
	psi = residue.get_psi()
	ome = residue.get_omega()
	residue.set_phi  (angle = -60)
	residue.set_psi  (angle = -40)
	residue.set_omega(angle = -40)

	#set_phi_psi_dihedral (molecule = system, resi= residue.id , bond='PHI', angle = -60 )
	#set_phi_psi_dihedral (molecule = system, resi= residue.id , bond='PSI', angle = -40 )
#'''	
	
'''
for residue in system.residues:
	phi = residue.get_phi()
	psi = residue.get_psi()
	ome = residue.get_omega()
	
	
	
	
	#print (residue.name, residue.id, phi, psi, ome)
	#print(residue.atoms[0].resn,  residue.atoms[0].pos, residue.atom_dic['N'].pos,  residue.atoms[0].pos[0] - residue.atom_dic['N'].pos[0] )
#'''

#
#psi = system.residues[17].get_psi()
#
#system.residues[17].psi_restraint_angle  = psi
#system.residues[17].psi_restraint_weight = 10.0
#
#for i in range (1, 360, 5):
#	system.residues[17].set_psi(psi+i)
#	print(system.residues[17].psi ,psi+i, system.energy())
#	save_XYZ_to_file(system , 'xyz_out_1zdd.xyz')
#
#phi = system.residues[17].get_phi()
#
#system.residues[17].phi_restraint_angle  = phi
#system.residues[17].phi_restraint_weight = 10.0
#
#for i in range (1, 360, 5):
#	system.residues[17].set_phi(phi+i)
#	print(system.residues[17].phi ,phi+i, system.energy())
#	save_XYZ_to_file(system , 'xyz_out_1zdd.xyz')
	
	
	
	
	
	
	
	
	
	
	
	
