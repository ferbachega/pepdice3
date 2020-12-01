
import numpy as np
import math

from pprint import pprint
from pepCore.Vectors import *



def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation
    about
    the given axis by theta radians.
    """
    axis  = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2)
    b, c, d = -axis * math.sin(theta / 2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def center_atom(molecule, x, y, z):
    '''the function that centers on the x, y and z coordinates  '''
    for atom_i in molecule.atoms:
        atom_i.pos[0] = atom_i.pos[0] - x
        atom_i.pos[1] = atom_i.pos[1] - y
        atom_i.pos[2] = atom_i.pos[2] - z
    return molecule


def move_atom(molecule=None, index=None, delta_x=0, delta_y=0, delta_z=0):
    """ Function that makes the straight displacement of an atom """
    atom = molecule.atoms[index]
    atom.pos[0] = atom.pos[0] + delta_x
    atom.pos[1] = atom.pos[1] + delta_y
    atom.pos[2] = atom.pos[2] + delta_z


def Compute_Rab(atom_a, atom_b):
    """ Function doc
    
    Calculates the vector R that is the subtraction: vector a - vector b    
    
    atom1, atom2  are atom objects

    """
    x_a = atom_a.pos[0]
    x_b = atom_b.pos[0]

    y_a = atom_a.pos[1]
    y_b = atom_b.pos[1]

    z_a = atom_a.pos[2]
    z_b = atom_b.pos[2]

    x = x_a - x_b
    y = y_a - y_b
    z = z_a - z_b
    
    return [x, y, z]


def distance_ab (atom1, atom2):
    """ Calculates the Euclidean distance between two atoms
    
    atom1, atom2  are atom objects
    
    """
    da = atom1.pos[0] - atom2.pos[0]
    db = atom1.pos[1] - atom2.pos[1]
    dc = atom1.pos[2] - atom2.pos[2]

    distante  = da**2 + db**2 + dc**2
    return distante**0.5



"""

        M A I N     C H A I N      R O T A T I O N S

"""

def rotate_Calpha_dihedral(molecule, axis, theta, window):
    """ Function doc """
    #print molecule, axis, theta, window
    for residue_i in molecule.residues:
        for atom_i in residue_i.atoms:

            if atom_i.id in window:
                #print atom_i.id
                atom_i.pos = np.dot(rotation_matrix(axis,
                                                    theta),
                                                    atom_i.pos
                                                    )


def rotate_backbone(molecule=None, resi=1, bond='PSI', theta=0, steps=1):
    """ Function doc """

    index =  resi
    if resi in molecule.fixed_residues:
        return 0


    #print molecule, resi, bond

    sideChain = []
    #print resi
    residue = molecule.residues[resi]
    for atom_i in residue.atoms:
        H = False
        if atom_i.name == 'CA':
            CA = atom_i

        elif atom_i.name == 'N':
            N = atom_i
                                                            #|---amber---|
        elif atom_i.name in ['H']:#, 'HT1', 'HT2', 'HT3', 'HN','H1','H2','H3']:
            H = atom_i

        elif atom_i.name == 'C':
            C = atom_i

        elif atom_i.name == 'O':
            O = atom_i
        
        elif atom_i.name == 'OXT':
            O = atom_i


        else:
            sideChain.append(atom_i.id)

    if bond == 'PSI':
        subcoord = CA.pos
        x = subcoord[0]
        y = subcoord[1]
        z = subcoord[2]

        molecule = center_atom(molecule=molecule, x=x, y=y, z=z)
        axis = C.pos
        window = range(0, CA.id + 1)
        window = list(window)
        for i in sideChain:
            window.append(i)
        
        if H:
            window.append(H.id)
        
        #try:
        #    for i in sideChain:
        #        window.append(i)
        #except:
        #    pass

        for i in range(0, steps):
            rotate_Calpha_dihedral(molecule, axis, theta, window)

    if bond == 'PHI':
        subcoord = CA.pos  # N - nitrogenio
        x = subcoord[0]
        y = subcoord[1]
        z = subcoord[2]
        molecule = center_atom(molecule=molecule, x=x, y=y, z=z)

        axis = N.pos     # CA - C alphal
        window = range(0, N.id + 1)
        try:
            window.append(H.id)
        except:
            pass

        for i in range(0, steps):
            rotate_Calpha_dihedral(molecule, axis, theta, window=window)

    if bond == "OMEGA":
        residue = molecule.residues[resi+1]
        for atom_i in residue.atoms:
            if atom_i.name == 'CA':
                CA2 = atom_i
            elif atom_i.name == 'N':
                N2 = atom_i
            elif atom_i.name in ['H', 'HT1', 'HT2', 'HT3', 'HN']:
                H2 = atom_i
            elif atom_i.name == 'C':
                C2 = atom_i
            elif atom_i.name == 'O':
                O2 = atom_i

            else:
                sideChain.append(atom_i.id)

        subcoord = C.pos  # N - nitrogenio
        x = subcoord[0]
        y = subcoord[1]
        z = subcoord[2]
        molecule = center_atom(molecule=molecule, x=x, y=y, z=z)

        axis = N2.pos
        window = range(0, N2.id + 1)
        #try:
        #    window.append(H2.id)
        #except:
        #    pass

        for i in range(0, steps):
            #print axis, theta, window
            rotate_Calpha_dihedral(molecule, axis, theta, window=window)


def set_phi_psi_dihedral (molecule=None, resi=1, bond='PSI', angle = 0.0 ):
    """ Function doc """

    #print resi, molecule.residues[resi].name
    #print bond
    #print angle


    initial_angle =  computePhiPsi(molecule = molecule,
                                   resi     = resi,
                                   bond      = bond)
    #if initial_angle == 0.0:
    #    initial_angle = 180.0

    if initial_angle == None:
        #print bond, 'False'
        return False

    if bond == 'PSI':
        delta_angle = (angle -  initial_angle) * -1
        theta       = math.radians(delta_angle)

    elif bond == 'OMEGA':
        delta_angle = (angle -  initial_angle) * -1
        theta       = math.radians(delta_angle)


    else:
        delta_angle = (angle -  initial_angle)
        theta       = math.radians(delta_angle)

    rotate_backbone(molecule = molecule   ,
                    resi     = resi      ,
                    bond     = bond       ,
                    theta    = theta,
                    steps    = 1)

    final_angle = computePhiPsi     (molecule  = molecule,
                                      resi     = resi,
                                      bond     = bond)

    #print bond, resi, molecule.residues[resi].name, initial_angle, final_angle, angle, delta_angle
    return final_angle






def computePhiPsi (molecule = None,  # obsolete
                   resi     =    1, 
                   bond     = 'PSI'):
    
    """ 
    obsolete
    
    """
    C1  = None
    N2  = None
    CA2 = None
    C2  = None
    N3  = None

    #----------------------------------------------
    # obtaining the C positions residue n 
    residue = molecule.residues[resi]
    CA2     = residue.CA
    N2      = residue.N
    C2      = residue.C
    #----------------------------------------------

    #CA2     = residue.atom_dic['CA']
    #N2      = residue.atom_dic['N']
    #C2      = residue.atom_dic['C']




    #----------------------------------------------
    # obtaining the C positions residue n - 1
    if resi == 0:
        phi = False
        C1  = None
        CA1 = None
        pass

    else:
        residue = molecule.residues[resi - 1]
        C1   = residue.C
        CA1  = residue.CA
        
        #C1   = residue.atom_dic['C']
        #CA1  = residue.atom_dic['CA']
        
        phi = True
    #----------------------------------------------


    #----------------------------------------------
    # obtaining the C positions residue n + 1
    try:
        residue = molecule.residues[resi + 1]
        N3   = residue.N
        CA3  = residue.CA
        
        #N3   = residue.atom_dic['N']
        #CA3  = residue.atom_dic['CA']
        
        psi = True
        ome = True
    except:
        psi = False
        ome = False
    #----------------------------------------------


    #print phi, psi, bond
    if phi:
        if bond == 'PHI':
            angle = dihedral(C1.pos  ,
                             N2.pos  ,
                             CA2.pos ,
                             C2.pos  )
            return math.degrees(angle)
    
    if psi:
        if bond == 'PSI':
            angle = dihedral(
                             N2.pos  ,
                             CA2.pos ,
                             C2.pos  ,
                             N3.pos)
            
            return math.degrees(angle)

    if ome:
        if bond == 'OMEGA':
            angle = dihedral(
                             CA2.pos, 
                             C2.pos ,
                             N3.pos ,
                             CA3.pos)
            
            
            
            if math.degrees(angle) <= 1 and math.degrees(angle) >= -1:
                # se o diedro é igual a 180 e funcao dihedral retorna 0
                return 180
            else:
                return math.degrees(angle)#angle*(180/math.pi)#57.29577951308232













