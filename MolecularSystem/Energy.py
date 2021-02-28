import math 
from pepBabel.PDBFiles import *
from pepCore.Geometry import distance_ab
#from parameters.RAMA.rama_left_ALL import RAMA
import os

class Energy:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
        #self.energy_components = {'e_cmap_short' : 0,
        #                          'e_cmap_medium': 0,
        #                          'e_cmap_long'  : 0,
        #                          'e_RW'         : 0,
        #                          }

    def read_calRW_logfile (self, logfile):
        """ Function doc """
        logfile = open(logfile, 'r')
        energy = None
        for line in logfile:
            line2 =  line.split()
            if line2[0] == "RW":
                energy = line2[-2]
            else:
                pass
        
        return energy
                

    def get_calRW_energy (self, theard_number = 0):
        """ Function doc """
        current_directory =  os.getcwd()
        

        os.chdir(self.RW_PATH)
        tmp_PDB = "TMP/" + str(theard_number) + "tmp.pdb"
        tmp_log = "TMP/" + str(theard_number) + "tmp.log"
        
        
        
        write_PDB_file(self, tmp_PDB)

        cmd  = './calRW ' + tmp_PDB + " > " + tmp_log
        os.system(cmd)
        energy = self.read_calRW_logfile(tmp_log)
        
        os.chdir(current_directory)
        self.energy_components['e_RW' ] = float(energy)
        return float(energy)


    def get_phi_psi_energy (self):
        """ Function doc 
        
        This function uses data from http://dunbrack.fccc.edu/ndrd/ndrd_PLOSCB.pdf
        
        """
        phi_psi_list = self.get_phi_psi_list()
        
        e_rama_total = 0
        c = 1
        for phi_psi in phi_psi_list[1:-1]:
            phi = phi_psi[0]
            psi = phi_psi[1]
            
            phi_round = (round (phi/10))*10
            psi_round = (round (psi/10))*10
            
            if phi_round == 180:
                phi_round = -180
            
            if psi_round == 180:
                psi_round = -180
            
            energy = RAMA[self.residues[c].name][(phi_round,psi_round)]['logP']
            
            self.residues[c].e_rama = energy
            
            e_rama_total += energy
            
            #print (c,self.residues[c].name, phi, psi, phi_round, psi_round, energy)
            c += 1
        
        self.e_rama_total = e_rama_total
        return e_rama_total
        
    
    def get_secondary_structure_restraint_energy (self, log = False):
        """ Function doc """
        
        
        '''
        if log:
            text = '%-14s ' %('RESIDUE')
            text+= '%-7s  ' %('PHI')  
            text+= '%-14s  ' %('Restrainted-PHI')  
            text+= '%-14s  ' %('DELTA PHI')  
            text+= '%-14s  ' %('weight')  
            text+= '%-14s  ' %('phi_energy')  

            
            text+= '%-7s  ' %('PSI')  
            text+= '%-14s  ' %('Restrainted-PSI')  
            text+= '%-14s  ' %('DELTA PSI') 
            text+= '%-14s  ' %('weight')  
            text+= '%-14s  ' %('psi_energy')  
            
            text+= '%-14s  ' %('TOTAL ENERGY')  

            print (text)
        '''
        
        total_energy = 0
        energy_list  = []        
        
        for residue in self.residues:
            
            phi_energy = 0
            psi_energy = 0
            Kd = 1.0
            
            #------------------------------------------------------------------------------------
            #                P H I
            if  residue.phi == None:
                pass
            else:
                if residue.phi_restraint_angle != None:
                    delta_phi  =  residue.phi - ( residue.phi_restraint_angle)
                    K_phi      =  Kd *  residue.phi_restraint_weight
                    #phi_energy =  K_phi*(math.radians(delta_phi))**2

                    phi_energy =  K_phi*(1 - math.cos(math.radians(delta_phi)))#**2
            #------------------------------------------------------------------------------------
            
            #------------------------------------------------------------------------------------
            #                P S I
            if  residue.psi == None:
                pass
            else:
                if  residue.psi_restraint_angle != None:
                    delta_psi     =  residue.psi - ( residue.psi_restraint_angle) 
                    K_psi         =  Kd * residue.psi_restraint_weight
                    
                    #psi_energy    =  K_psi*(math.radians(delta_psi))**2

                    psi_energy    =  K_psi*(1 - math.cos(math.radians(delta_psi)))#**2
            #------------------------------------------------------------------------------------
            
            
            energy     =  phi_energy + psi_energy 
            total_energy += energy
        
        if log:
            print ('TOTAL ENERGY = ', total_energy)
            print (energy_list)
        
        #print (total_energy)
        return total_energy        
    
    
    def get_distance_restraint_energy (self, log = False):
        """ Function doc """
        pass
    
    
    def get_contact_map_energy (self, cmap =  None):
        """ Function doc """
		
        energy_short  = 0
        energy_medium = 0
        energy_long   = 0  
        for contact in self.contacts:
        
            residue1 =contact.residue1
            residue2 =contact.residue2
            
            
            # ATOM 1
            if residue1.name == 'GLY':
                atom1 =residue1.CA 
            else:
                atom1 =residue1.CB 

            # ATOM 2    
            if residue2.name == 'GLY':
                atom2 =residue2.CA 
            else:
                atom2 =residue2.CB 

           
            # get distance
            distance = distance_ab (atom1, atom2)
            
            
            
            if distance <=  contact.cutoff:
                if contact.type == 'C':
                    energy_short += -10.0*contact.weight
                
                if contact.type == 'M':
                    energy_medium += -10.0*contact.weight            
                
                if contact.type == 'L':
                    energy_long += -10.0*contact.weight
            
            else:
                pass
        
            #print (distance, energy_short, energy_medium, energy_long, residue1.name , residue2.name)
             
        self.energy_components['e_cmap_short' ] = energy_short
        self.energy_components['e_cmap_medium'] = energy_medium
        self.energy_components['e_cmap_long'  ] = energy_long
        self.energy_components['e_cmap_total' ] = energy_short +  energy_medium +  energy_long
        return  energy_short, energy_medium, energy_long  
    
 
    def energy (self, rw = True, cmap = True):
        """ Function doc """

        self.get_phi_psi_list()
        #e_rama = 0
        #e_rama         = self.get_phi_psi_energy ()
        #e_ss_restraint = self.get_secondary_structure_restraint_energy()
        
        self.current_energy = 0
        
        if cmap:
            energy_short, energy_medium, energy_long = self.get_contact_map_energy()
            self.current_energy += energy_short +  energy_medium +  energy_long
        
        
        if rw:
            e_RW  = self.get_calRW_energy()
            self.current_energy += e_RW

        #for energy_term in self.energy_components:
        #    print(energy_term, self.energy_components[energy_term])
			
        #print (e_rama + e_ss_restraint + e_RW)
        #self.current_energy = e_rama + e_ss_restraint + e_RW + e_contact_map
        #return e_rama + e_ss_restraint + e_RW + e_contact_map
        #self.current_energy = e_contact_map
        return self.current_energy
        
        
    
