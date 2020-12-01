







filein =  open('ffG43A1-S3.atp', 'r')
ATOM_TYPES     = {}
NONBOND_PARAMS = {}
PAIRTYPES      = {}
BONDS          = {}
ANGLES         = {}
#nonbond_params


for line in filein:
	line2 = line.split(';')
	line3 = line2[0].split()
	
	atom_name = line3[0]
	atom_mass = float(line3[1])
	#print (atom_name, atom_mass)
	ATOM_TYPES[line3[0]] = {}
	ATOM_TYPES[line3[0]]['mass'] = atom_mass

#print ATOM_TYPES



filein =  open('ffG43A1-S3.02.itp', 'r')

atomtypes_tag      = False
nonbond_params_tag = False
pairtypes_tag      = False
bonds_tag          = False 
angle_tag          = False
for line in filein:
	line2 =  line.split()
	
	if 'atomtypes' in line2:
		atomtypes_tag = True

	if atomtypes_tag:
		if len(line2) >=6:
			#print (line2)
			if 'mass' in line2:
				pass
			else:
				#print (line2)
				ATOM_TYPES[line2[0]]['C6']     = float(line2[4])
				ATOM_TYPES[line2[0]]['C12']    = float(line2[5])
				ATOM_TYPES[line2[0]]['ptype']  = line2[3]
				ATOM_TYPES[line2[0]]['charge'] = float(line2[2])
		else:
			pass
	
	
	
	
	
	if 'nonbond_params' in line2:
		atomtypes_tag      = False
		nonbond_params_tag =True
	
	
	if nonbond_params_tag:
		if len(line2) >=5:
			if 'func' in line2:
				pass
			else:
				print line2
				NONBOND_PARAMS[(line2[0],line2[1])]  = {'C6'  : float(line2[3]), 'C12' : float(line2[4]), 'func':int(line2[2]) }
		else:
			pass	
	
	
		
	if 'pairtypes' in line2:
		atomtypes_tag      = False
		nonbond_params_tag = False
		pairtypes_tag      = True
		
	if pairtypes_tag:
		if len(line2) >=5:
			if 'func' in line2:
				pass
			else:
				print line2
				PAIRTYPES[(line2[0],line2[1])]  = {'cs6'  : float(line2[3]), 'cs12' : float(line2[4]),  'func':int(line2[2]) }
		else:
			pass	





	if 'bond-stretching' in line2:
		atomtypes_tag      = False
		nonbond_params_tag = False
		pairtypes_tag      = False
		bonds_tag          = True
		
	
	if bonds_tag:
		if len(line2) >= 4:
			if '#define' in line2:
				BONDS[line2[1]] = {'b0': line2[2], 'cb': line2[3]} 
			else:
				
				pass
				#print line2
				#PAIRTYPES[(line2[0],line2[1])]  = {'cs6'  : float(line2[3]), 'cs12' : float(line2[4]),  'func':int(line2[2]) }
		else:
			pass


	if 'bond-angle' in line2:
		atomtypes_tag      = False
		nonbond_params_tag = False
		pairtypes_tag      = False
		bonds_tag          = False
		angle_tag          = True
	
	if angle_tag:
		if len(line2) >= 3:
			if '#define' in line2:
				ANGLES[line2[1]] = {'b0': line2[2], 'cb': line2[3]} 
			else:
				
				pass
				#print line2
				#PAIRTYPES[(line2[0],line2[1])]  = {'cs6'  : float(line2[3]), 'cs12' : float(line2[4]),  'func':int(line2[2]) }
		else:
			pass
			
			
			
			
print ATOM_TYPES['OM']['C6']
#filein = filein.readlines()
print PAIRTYPES[('OM','OM')]
print NONBOND_PARAMS[('OM','OM')]
print BONDS
