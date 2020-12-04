#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Atom.py
#  
#  Copyright 2016 farminf <farminf@farminf-3>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

from pprint import pprint
from pepBabel.PDBfiles import *
from pepCore.Geometry import *
from MolecularSystem.Energy import Energy


class Atom(object):
    """class to store info about an atom"""

    def __init__(self, 
                 molSystem = None     ,
                 id        = 0        ,
                 a_number  = None     ,
                 name      = 'Xx'     ,
                 element   = "UNK"    ,   
                 resn      = "UNK"    ,
                 resi      = 1        ,
                 sigma     = 0.00     ,
                 epsilon   = 0.00     ,   
                 charge    = 0.00     ,
                 pos       = [0, 0, 0]):

        # 119.8 K and  0.3405 nm)

        # PDB data
        
        self.id   = id    # index
        self.name = name  # atom name  "OH" or "CA"

        self.resn             = resn          #  residue name
        self.resi             = resi          #  residues index
        self.element          = element
        self.resobj           = None     # residue object
        
        self.molSystem = molSystem
        
        
        # FF paramters
        #self.charge         = charge  # carga
        #self.sigma          = sigma
        #self.epsilon        = epsilon
        #self.ff_atom_number = None
        #self.ff_atom_name   = None
        #self.atomic_number  = None
        #self.mass           = None

        # coordinates and energy
        
        self.pos           = pos
        self.actual_energy = 0.00
        self.sigma_ab      = 0.00


class Residue(object):
    """ Class doc """

    def __init__(self, molSystem = None, id = 0,  name = 'UNK', resi = 0):
        """ Class initialiser """
        self.id    = id
        self.atoms = []
        self.bonds = []
        self.name  = name
        self.resi  = resi
        self.molSystem = molSystem
        
        self.N  = None
        self.CA = None
        self.C  = None
        self.O  = None 
        self.H  = None
        self.CB = None
         
       
        self.e_rama  = 0
        self.phi = None
        self.psi = None
        
        
        # PHI
        self.phi_restraint_angle  = 0.0
        self.phi_restraint_weight = 0.0
        
        # PSI
        self.psi_restraint_angle  = 0.0
        self.psi_restraint_weight = 0.0

        # OMEGA
        self.omega_restraint_angle  = 0.0
        self.omega_restraint_weight = 0.0




        # CHI's
        
        self.CHI1 = None
        self.CHI2 = None
        self.CHI3 = None
        self.CHI4 = None
        self.CHI5 = None
        
        #self.chi1_restraint_angle  = 0.0
        #self.chi2_restraint_angle  = 0.0
        #self.chi3_restraint_angle  = 0.0
        #self.chi4_restraint_angle  = 0.0
        #self.chi5_restraint_angle  = 0.0
        #self.chi1_restraint_weight = 0.0
        #self.chi2_restraint_weight = 0.0
        #self.chi3_restraint_weight = 0.0
        #self.chi4_restraint_weight = 0.0
        #self.chi5_restraint_weight = 0.0

        #self.ss_restraint   = [None, None, None]
        #self.w_ss_restraint = [0.0 , 0.0,  0.0 ]

    def add_new_atom (self, atom = None):
        """ Function doc """
        self.atoms.append(atom)
        
        if atom.name == 'N':
            self.N =  atom
        
        if  atom.name == 'CA':
            #print (atom.name, atom.pos, atom.resi )
            self.CA =  atom
            #print (self.CA.name, self.CA.pos, self.CA.resi )
        if atom.name == 'C':
            self.C =  atom
        
        if  atom.name == 'O':
            self.O =  atom

        if  atom.name == 'CB':
            self.CB =  atom

        if  atom.name == 'H':
            self.H =  atom


    def build_atoms_dic (self):
        """ Function doc """
        self.atom_dic = AA_atoms[self.name]
        
        for atom in self.atoms:
            self.atom_dic[atom.name] = atom 
     
    def get_phi (self, radian = False):
        """ Function doc 
      

            Amino Acid 1              Amino Acid 2              Amino Acid 3           
        |-----------------------|-------------------------|------------------------


        
                  H       O               H        O               H        O
                  |       |               |        |               |        |
         N - - - CA - - - C - - - N - - - CA - - - C - - - N - - - CA - - - C - - - N
                  |                       |                        |
                  R                       R                        R


                         (C1)----(N2)---(CA2)----(C2)     
                                      ^        ^       ^ 
                                     PHI      psi     ome
        



        C1  = self.molSystem.residues[self.id-1].C  <---- residue before
        N2  = self.N 
        CA2 = self.CA
        C2  = self.C 

        """

        if self.id == 0:
            return None

        else:
            C1  = self.molSystem.residues[self.id-1].C
            N2  = self.N 
            CA2 = self.CA
            C2  = self.C 
            
            
            
            angle = dihedral(C1.pos  ,
                             N2.pos  ,
                             CA2.pos ,
                             C2.pos  )
            if radian:
                return angle
                self.phi = math.degrees(angle)
            else:
                self.phi = math.degrees(angle)
                return math.degrees(angle)

    def get_psi (self,  radian = False):
        """
      

            Amino Acid 1              Amino Acid 2              Amino Acid 3           
        |-----------------------|-------------------------|------------------------


        
                  H       O               H        O               H        O
                  |       |               |        |               |        |
         N - - - CA - - - C - - - N - - - CA - - - C - - - N - - - CA - - - C - - - N
                  |                       |                        |
                  R                       R                        R


                         (C1)----(N2)---(CA2)----(C2)----(N3)----(CA3)    
                                      ^        ^       ^ 
                                     phi      PSI     ome
        

        """
        if self.id == len(self.molSystem.residues) -1:
            return None

        else:
            N2  = self.N 
            CA2 = self.CA
            C2  = self.C 
            N3  = self.molSystem.residues[self.id+1].N
            
            #print (N2, 
            #       CA2,
            #       C2, 
            #       N3, self.id, self.resi, self.CA, self.atoms) 
            
            angle = dihedral(
                         N2.pos  ,
                         CA2.pos ,
                         C2.pos  ,
                         N3.pos)
            if radian:
                return angle
                self.psi = math.degrees(angle)
            else:
                self.psi = math.degrees(angle)
                return math.degrees(angle)      
        
    def get_omega (self,  radian = False):
        """ Function doc """
        if self.id == len(self.molSystem.residues) -1:
            return None
        
        else:
            CA2 = self.CA
            C2  = self.C 
            N3  = self.molSystem.residues[self.id+1].N
            CA3 = self.molSystem.residues[self.id+1].CA

            angle = dihedral(
                             CA2.pos, 
                             C2.pos ,
                             N3.pos ,
                             CA3.pos)
            if radian:
                return angle
            else:           
                if math.degrees(angle) <= 1 and math.degrees(angle) >= -1:
                    # se o diedro é igual a 180 e funcao dihedral retorna 0
                    return 180
                else:
                    return math.degrees(angle)#angle*(180/math.pi)#57.29577951308232

    def set_phi (self, angle = 0.0 ):
        """ Function doc """
        
        initial_angle =  self.get_phi()
        if initial_angle == None:
            return False
        
        else:
            delta_angle = (angle -  initial_angle)
            theta       = math.radians(delta_angle)

        rotate_backbone(molecule = self.molSystem ,
                        resi     = self.id      ,
                        bond     = "PHI"     ,
                        theta    = theta,
                        steps    = 1)

        final_angle =  self.get_phi()
        
        return final_angle
            
    def set_psi (self, angle = 0.0 ):
        initial_angle =  self.get_psi()
        if initial_angle == None:
            return False
        else:
            delta_angle = (angle -  initial_angle) * -1
            theta       = math.radians(delta_angle)

        rotate_backbone(molecule = self.molSystem ,
                        resi     = self.id      ,
                        bond     = "PSI"     ,
                        theta    = theta,
                        steps    = 1)

        final_angle =  self.get_psi()
        return final_angle
   
    def set_omega (self, angle = 0.0):
        """ Function doc """
        initial_angle =  self.get_omega()
        if initial_angle == None:
            return False
        else:
            delta_angle = (angle -  initial_angle) * -1
            theta       = math.radians(delta_angle)
        rotate_backbone(molecule = self.molSystem ,
                        resi     = self.id      ,
                        bond     = "OMEGA"     ,
                        theta    = theta,
                        steps    = 1)
        final_angle =  self.get_omega()
        return final_angle


class Molecule(#Atom       ,
               #Residue    ,
               #Coordinates,
               Energy     ,
               #ModelAB    ,
               ):
    """ 
    class to store info about a molecule

    energy models:
    amber   - amber ff03ua energy components
    Calpha  - C alpha model  - like AB model
    Contact - Contact map only
    LPFSF   - LABIO protein folding scoring function

    """

    def __init__(self, 
                 id   = 0,
                 name = 'protein',
                 pdb  = None,
                 ):
                     
        self.id       = id
        self.name     = name
        self.atoms    = [] # a list of atom objects
        self.residues = []
        self.chains   = []

        self.load_PDB_to_system(filename = pdb)
        

        #self.top      = None
        #self.psf      = None
        #self.param    = None
        #self.ff_type  = None

        #self.torsions       = None
        #self.FIX_atoms_CHI  = None


        # Parameters and Restraints
        # -------------------------------------------- #
        self.fixed_residues                     = []
        self.fragments                          = []
        self.cmap                               = None 

        self.hamonical_potential_restraint_list = []


        # MC important atributes
        self.actual_energy   = None
        self.previous_energy = None 


        self.energy_model = 'LABIO'
        
        self.phi_psi_list = []
        self.e_rama_total = None



        self.aa_dic = { 
          'A' : 'ALA',
          'R' : 'ARG',
          'N' : 'ASN',
          'D' : 'ASP',
          'C' : 'CYS',
          'E' : 'GLU',
          'Q' : 'GLN',
          'G' : 'GLY',
          'H' : 'HIS',
          'I' : 'ILE',
          'L' : 'LEU',
          'K' : 'LYS',
          'M' : 'MET',
          'F' : 'PHE',
          'P' : 'PRO',
          'S' : 'SER',
          'T' : 'THR',
          'W' : 'TRP',
          'Y' : 'TYR',
          'V' : 'VAL'
          }


        self.AA_dic = { 
          'ALA': 'A',
          'ARG': 'R',
          'ASN': 'N',
          'ASP': 'D',
          'CYS': 'C',
          'GLU': 'E',
          'GLN': 'Q',
          'GLY': 'G',
          'HIS': 'H',
          'HIE': 'H',
          'HID': 'H',
          'ILE': 'I',
          'LEU': 'L',
          'LYS': 'K',
          'MET': 'M',
          'PHE': 'F',
          'PRO': 'P',
          'SER': 'S',
          'THR': 'T',
          'TRP': 'W',
          'TYR': 'Y',
          'VAL': 'V'
         }

    
    def get_phi_psi_list (self):
        """ Function doc """
        phi_psi_list = []
        
        for residue in self.residues:
            phi = residue.get_phi()
            psi = residue.get_psi()
            
            phi_psi_list.append((phi, psi))
        
        self.phi_psi_list = phi_psi_list
        return phi_psi_list
            

    def load_PDB_to_system(self, filename = None):
        
        atoms_raw_list  = PDBParser(filename)
        self.build_the_atom_object_list(atoms_raw_list)
        self.build_the_residue_object_list()

    
    def build_the_atom_object_list (self, atoms_raw_list):
        """ Function doc 
        
        
        
           - - - -  building the atoms list  - - - - 
        
        
        """
        
    
        n = 0
        for atom_raw in atoms_raw_list:
            #[424, 'OG1', None, [-41.551, 14.604, 22.51], 55, 'THR', 'A', ' O', [], None, 1.0, 0.0, 0.0]
            index    =  atom_raw[1]
            name     =  atom_raw[1]
            pos      =  atom_raw[3]             
            resi     =  atom_raw[4]             
            resn     =  atom_raw[5]
            chain    =  atom_raw[6]
            element  =  atom_raw[7]
            
            atom  =  Atom( molSystem = self     ,
                           id        = n         ,
                           a_number  = index     ,
                           name      = name      ,
                           element   = element   ,   
                           resn      = resn      ,
                           resi      = int(resi) ,
                           sigma     = 0.00      ,
                           epsilon   = 0.00      ,   
                           charge    = 0.00      ,
                           pos       = pos       )
            
            
            self.atoms.append(atom)
            n += 1
            #print (atom)  

    def build_the_residue_object_list (self):
        """ Function doc """
        #
        #    - - - -  building the residues list  - - - - 
        #
        residue       = None
        new_residue   = True
        index_counter = 0
        id_counter    = 0 
        for atom  in self.atoms:
            if new_residue:
                # verifica se um novo resíduo deve ser criado
                resn          =  atom.resn 
                resi          =  atom.resi
                index_counter =  resi
                residue =  Residue ( molSystem = self , id = id_counter,  name = resn, resi = resi)

                #residue =  Residue ( molSystem = self , id = resi,  name = resn)
                residue.add_new_atom(atom)
                
                new_residue = False
                id_counter    += 1
            
            else:

                resn          =  atom.resn 
                resi          =  atom.resi
                
                if index_counter == resi:
                    residue.add_new_atom(atom)
                    


                else:
                    # aciona o resíduo já devidamente preenchido na lista de resíduos
                    self.residues.append(residue)
                    index_counter =  resi
                    residue =  Residue ( molSystem = self , id = id_counter,  name = resn , resi = resi)
                    residue.add_new_atom(atom)
                    id_counter    += 1

        self.residues.append(residue)
        

        for residue in self.residues:
            residue.build_atoms_dic()
            
            for atom in residue.atoms:
                
                if atom.name in ['N','CA', 'C', 'O']:
                    
                    if atom.name == 'N':
                        print (residue.id, residue.resi, atom.resi, atom.resn,  atom.name, atom.pos, residue.name, residue.N.name, residue.N.pos, atom.pos[0] - residue.N.pos[0]    , atom.pos[0] - residue.atom_dic["N"].pos[0])
                    
                    elif atom.name == 'CA':
                        print (residue.id, residue.resi, atom.resi, atom.resn,  atom.name, atom.pos, residue.name, residue.CA.name, residue.CA.pos, atom.pos[0] - residue.CA.pos[0], atom.pos[0] - residue.atom_dic["CA"].pos[0])
                    
                    elif atom.name == 'C':
                        print (residue.id, residue.resi,atom.resi, atom.resn,  atom.name, atom.pos, residue.name, residue.C.name, residue.C.pos,   atom.pos[0] - residue.C.pos[0], atom.pos[0] - residue.atom_dic["C"].pos[0])
                        
                    elif atom.name == 'O':
                        print (residue.id, residue.resi,atom.resi, atom.resn,  atom.name, atom.pos, residue.name, residue.O.name, residue.O.pos,   atom.pos[0] - residue.O.pos[0] , atom.pos[0] - residue.atom_dic["O"].pos[0])
                        
                    else:
                        print (residue.id, residue.resi,atom.resi, atom.resn,  atom.name, atom.pos)
                
 
    def import_coordinates_to_system(self, coordinates):
        """ Function doc """
        # for key in self.atoms.keys():

        counter = 0
        for  atom in self.atoms:
            atom.pos[0] = coordinates[counter][0]
            atom.pos[1] = coordinates[counter][1]
            atom.pos[2] = coordinates[counter][2]
            counter += 1
            
            
        
        #for residue_i in self.residues:
        #    for atom_i in residue_i.atoms:
        #        index_i = atom_i.id
        #        atom_i.pos[0] = coordinates[index_i - 1][0]
        #        atom_i.pos[1] = coordinates[index_i - 1][1]
        #        atom_i.pos[2] = coordinates[index_i - 1][2]

    def get_coordinates_from_system(self):  # molecule = None):
        """ Function doc """
        # pprint(ff.biotypes)
        
        initial_coordinates = []
        for atom in self.atoms:
            initial_coordinates.append([atom.pos[0], atom.pos[1], atom.pos[2]])

        return initial_coordinates
        
        #initial_coordinates = []
        #for residue_i in self.residues:
        #    for atom_i in residue_i.atoms:
        #        initial_coordinates.append(
        #            [atom_i.pos[0], atom_i.pos[1], atom_i.pos[2]])
		#
        #return initial_coordinates
        

































