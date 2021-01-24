




AA_atoms = {
                "GLY": { # 1
			'N'  : None, 
			'CA' : None, 
			'C'  : None, 
			'O'  : None,     
		       },
		       
		"ALA": { # 2
			'N'  : None,
			'CA' : None,
			'C'  : None,
			'O'  : None,
			'CB' : None,
			},
		
		"MET": {
			'N'  : None,
			'CA' : None,
			'C'  : None,
			'O'  : None,
			'CB' : None,
			'CG' : None,
			'SD' : None,
			'CE' : None,
			},
		
		"THR": {
			'N'   : None,
			'CA'  : None,
			'C'   : None,
			'O'   : None,
			'CB'  : None,
			'OG1' : None,
			'CG2' : None,
			},
		
		"LYS": {
			'N'   : None,  
			'CA'  : None,  
			'C'   : None,  
			'O'   : None,  
			'CB'  : None,  
			'CG'  : None,  
			'CD'  : None,  
			'CE'  : None,  
			'NZ'  : None,  
			},
		
		"LEU": {
			'N'   : None, 
			'CA'  : None, 
			'C'   : None, 
			'O'   : None, 
			'CB'  : None, 
			'CG'  : None, 
			'CD1' : None, 
			'CD2' : None, 
			},
		
		"ILE": {
			'N'   : None,
			'CA'  : None,
			'C'   : None,
			'O'   : None,
			'CB'  : None,
			'CG1' : None,
			'CG2' : None,
			'CD1' : None,
			},
			
		"GLU": {
			'N'   : None, 
			'CA'  : None, 
			'C'   : None, 
			'O'   : None, 
			'CB'  : None, 
			'CG'  : None, 
			'CD'  : None, 
			'OE1' : None, 
			'OE2' : None, 
			},	
		
		"VAL": {
			'N'   : None,   
			'CA'  : None,   
			'C'   : None,   
			'O'   : None,   
			'CB'  : None,   
			'CG1' : None,   
			'CG2' : None,
			},
			
		
		"ASP": {
			'N'   : None, 
			'CA'  : None, 
			'C'   : None, 
			'O'   : None, 
			'CB'  : None, 
			'CG'  : None, 
			'OD1' : None, 
			'OD2' : None,
			},			
		
		"PHE": {
			'N'   : None, 
			'CA'  : None, 
			'C'   : None, 
			'O'   : None, 
			'CB'  : None, 
			'CG'  : None, 
			'CD1' : None, 
			'CE1' : None, 
			'CD2' : None, 
			'CE2' : None, 
			'CZ'  : None, 
			},			
					
		
		"GLN": {
			'N'   : None, 
			'CA'  : None, 
			'C'   : None, 
			'O'   : None, 
			'CB'  : None, 
			'CG'  : None, 
			'CD'  : None, 
			'OE1' : None, 
			'NE2' : None, 
			},			
							
		"TYR": {
			'N'   : None, 
			'CA'  : None, 
			'C'   : None, 
			'O'   : None, 
			'CB'  : None, 
			'CG'  : None, 
			'CD1' : None, 
			'CE1' : None, 
			'CD2' : None, 
			'CE2' : None, 
			'CZ'  : None, 
			'OH'  : None, 
			},			
		
		"ASN": {
			'N'   : None, 
			'CA'  : None, 
			'C'   : None, 
			'O'   : None, 
			'CB'  : None, 
			'CG'  : None, 
			'OD1' : None, 
			'ND2' : None, 
			},			
		
		"TRP": {
			'N'   : None, 
			'CA'  : None, 
			'C'   : None, 
			'O'   : None, 
			'CB'  : None, 
			'CG'  : None, 
			'CD1' : None, 
			'CD2' : None, 
			'NE1' : None, 
			'CE2' : None, 
			'CE3' : None, 
			'CZ2' : None, 
			'CZ3' : None, 
			'CH2' : None,
			},				
		
		
		
		
		"ARG": {
			'N'   : None,
			'CA'  : None,
			'C'   : None,
			'O'   : None,
			'CB'  : None,
			'CG'  : None,
			'CD'  : None,
			'NE'  : None,
			'CZ'  : None,
			'NH1' : None,
			'NH2' : None,
			},
		    
		"SER": {
			'N'   : None, 
			'CA'  : None, 
			'C'   : None, 
			'O'   : None, 
			'CB'  : None, 
			'OG'  : None,
			},
		
		"PRO": {
			'N'   : None, 
			'CA'  : None, 
			'C'   : None, 
			'O'   : None, 
			'CB'  : None, 
			'CG'  : None, 
			'CD'  : None, 
			},
		
		"HIS": {
			'N'   : None,  
			'CA'  : None,  
			'C'   : None,  
			'O'   : None,  
			'CB'  : None,  
			'CG'  : None,  
			'CD2' : None,  
			'ND1' : None,  
			'CE1' : None,  
			'NE2' : None, 
			},

		"CYS": {
			'N'  : None, 
			'CA' : None, 
			'C'  : None, 
			'O'  : None, 
			'CB' : None, 
			'SG' : None, 
			},
		}





def PDBParser(pdbin):
    """ Function doc 

ATOM      1  N   THR A   1      -1.820  24.919  -5.344  1.00  0.00           N  
ATOM      2  CA  THR A   1      -1.256  24.379  -4.074  1.00  0.00           C  
ATOM     61  CB  ILE A   4      -7.386  -0.466   0.343  1.00  0.00           C  
ATOM  16     HO2 GLC  1       2.188   -0.704  0.939   1.00  300.00          H 0.0000   
ATOM     27  NH1 ARG A   5      68.029  23.029  29.719  1.00 38.75           N1+
    """
    pdb_file_lines  = open(pdbin, 'r')
    #pdb_file_lines  = rawframe.split('\n')   
    pdb_file_lines = pdb_file_lines.readlines()
    atoms           = []
    index           = 0
    
    for line in pdb_file_lines:
    
        if line[:4] == 'ATOM' or line[:6] == 'HETATM':
            
            at_name    = line[12:16].strip()
            at_pos     = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            at_resi    = int(line[22:27])
            at_resn    = line[17:20].strip()
            at_ch      = line[21]             
            at_symbol  = line[76:78]
                
            try:
                at_occup   = float(line[54:60])   #occupancy
            except:
                at_occup   = 0.00

            try:
                at_bfactor = float(line[60:66])
            except:
                at_bfactor = 0.00

            
            at_charge  = 0.0
            
            cov_rad  = None#at.get_cov_rad (at_name)
            gridpos  = None #[int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
            #ocupan   = float(line[54:60])
            #bfactor  = float(line[60:66])
            
                            #0      1        2        3       4        5        6       7       8       9       10          11        12      
            atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
            index += 1

    return atoms



def write_PDB_file (molecule , pdbout):
	""" Function doc 
	1 	"ATOM " or "HETATM" 	        				6 	{:6s} 	01-06 	[0:6]
	2 	atom serial number 								5 	{:5d} 	07-11 	[6:11]
	3 	atom name 										4 	{:^4s} 	13-16 	[12:16]
	4 	alternate location indicator 					1 	{:1s} 	17 	[16:17]
	5 	residue name 									3 	{:3s} 	18-20 	[17:20]
	6 	chain identifier 								1 	{:1s} 	22 	[21:22]
	7 	residue sequence number 						4 	{:4d} 	23-26 	[22:26]
	8 	code for insertion of residues 					1 	{:1s} 	27 	[26:27]
	9 	orthogonal coordinates for X (in Angstroms) 	8 	{:8.3f} 	31-38 	[30:38]
	10 	orthogonal coordinates for Y (in Angstroms) 	8 	{:8.3f} 	39-46 	[38:46]
	11 	orthogonal coordinates for Z (in Angstroms) 	8 	{:8.3f} 	47-54 	[46:54]
	12 	occupancy 										6 	{:6.2f} 	55-60 	[54:60]
	13 	temperature factor 								6 	{:6.2f} 	61-66 	[60:66]
	14 	element symbol 									2 	{:>2s} 	77-78 	[76:78]
	15 	charge on the atom 								2 	{:2s} 	79-80 	[78:80]
	
	
	
	
	
	"""
	
	text = []
	
	for atom in molecule.atoms:
		line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format("ATOM"     ,
																															   atom.id+1  , 
																															   atom.name  , 
																															   ""         ,
																															   atom.resn  ,
																															   "A"        ,
																															   atom.resi  ,
																															   ""         ,
																															   atom.pos[0],
																															   atom.pos[1],
																															   atom.pos[2],
																															   0.00       ,
																															   0.00       ,
																															   atom.element,
																															   ""
																															   		
																															  )
		#print (line)
	
		text.append(line)
	
	
	fileout = open(pdbout , 'w')
	fileout.writelines(text)
	fileout.close()
	
	
	
	





#atoms = PDBParser('/home/fernando/Downloads/predicted_structure.pdb')
#print (atoms)








#cpdef get_list_of_atoms_from_rawframe(rawframe, gridsize = 3, at =  None):
#    """ Function doc 
#
#ATOM      1  N   THR A   1      -1.820  24.919  -5.344  1.00  0.00           N  
#ATOM      2  CA  THR A   1      -1.256  24.379  -4.074  1.00  0.00           C  
#ATOM     61  CB  ILE A   4      -7.386  -0.466   0.343  1.00  0.00           C  
#ATOM  16     HO2 GLC  1       2.188   -0.704  0.939   1.00  300.00          H 0.0000   
#ATOM     27  NH1 ARG A   5      68.029  23.029  29.719  1.00 38.75           N1+
#
#iCode = ""
#( serial       ,
#  name         ,
#  altLoc       ,
#  resName      ,
#  chainID      ,
#  resSeq       ,
#  u00          ,
#  u11          ,
#  u22          ,
#  u01          ,
#  u02          ,
#  u12          ,
#  segID        ,
#  atomicNumber ,
#  formalCharge ) = self._ParseFixedFormatLine ( line                    ,
#                                                (  6, 11, int  , None ) ,
#                                                ( 12, 16, None , ""   ) ,
#                                                ( 16, 17, None , ""   ) ,
#                                                ( 17, 20, None , ""   ) ,
#                                                ( 21, 22, None , ""   ) ,
#                                                ( 22, 27, int  , None ) ,
#                                                ( 28, 35, float, None ) ,
#                                                ( 35, 42, float, None ) ,
#                                                ( 42, 49, float, None ) ,
#                                                ( 49, 56, float, None ) ,
#                                                ( 56, 63, float, None ) ,
#                                                ( 63, 70, float, None ) ,
#                                                ( 72, 76, None , ""   ) ,
#                                                ( 76, 78, PeriodicTable.AtomicNumberFromSymbol, -1 ) ,
#                                                ( 78, 80, None , ""   ) )
#
#    """
#
#    pdb_file_lines  = rawframe.split('\n')   
#    atoms           = []
#    cdef int index           = 0
#    for line in pdb_file_lines:
#        if line[:4] == 'ATOM' or line[:6] == 'HETATM':
#            
#            at_name    = line[12:16].strip()
#            at_pos     = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
#            at_resi    = int(line[22:27])
#            at_resn    = line[17:20].strip()
#            at_ch      = line[21]             
#            at_symbol  = line[76:78]
#            at_occup   = float(line[54:60])   #occupancy
#            at_bfactor = float(line[60:66])
#            at_charge  = 0.0
#            
#            cov_rad  = at.get_cov_rad (at_name)
#            gridpos  = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
#            #ocupan   = float(line[54:60])
#            #bfactor  = float(line[60:66])
#            
#                            #0      1        2        3       4        5        6       7       8       9       10          11        12      
#            atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
#            index += 1
#
#    return atoms
