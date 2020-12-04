import math 
#from parameters.RAMA.rama_left_ALL import RAMA

class Energy:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
    
        pass



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
    
    
    def energy (self):
        """ Function doc """
        
        self.get_phi_psi_list()
        e_rama = 0
        #e_rama         = self.get_phi_psi_energy ()
        e_ss_restraint = self.get_secondary_structure_restraint_energy()
    
        
        
        return e_rama + e_ss_restraint
        
        
    
