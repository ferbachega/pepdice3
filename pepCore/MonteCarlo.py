
#import numpy as np
#import math
from pepCore.Geometry                 import *
from pepBabel.XYZFiles import *
#from pprint import pprint
#from pepCore.Vectors import *



#-------------------------------------------
import math
#-------------------------------------------

def insert_fragment (molecule = None, fragment = None, sidechain = False):
    """ Function doc """
    if fragment is None:
        return
    else:
    #for key in fragment.keys():
        #PSI = fragment[key]['PSI']
        #PHI = fragment[key]['PHI']
        print (fragment.keys(), fragment)
        
        for residue_index in fragment.keys():
            
            index = int(residue_index)
            phi   = fragment[residue_index]["PHI"]
            psi   = fragment[residue_index]["PSI"]
            omega = fragment[residue_index]["OMEGA"]
            #print ()
            #print ("\n\n", index, phi,psi,omega,"\n\n")
            
            #print ('\n\n\fragment')
            #print (fragment)
            #print ( molecule.residues[index],  index, phi, psi, omega)
            #print ('\n\n\ ')

            molecule.residues[index].set_phi  (angle = phi  )
            molecule.residues[index].set_psi  (angle = psi  )
            molecule.residues[index].set_omega(angle = omega)
        
        
        #for bond in ['PSI','PHI']:
        #    #print  key, bond,  fragment[key][bond]
        #    molecule.set_phi
        #    
        #    
        #    
        #    set_phi_psi_dihedral (molecule = molecule,
        #                              resi = key     ,
        #                              bond = bond    ,
        #                              angle = fragment[key][bond])

        
        #if sidechain:
        #    #try:
        #    for bond in ['CHI1','CHI2','CHI3','CHI4','CHI5']:
        #        #try:
        #        if bond in fragment[key]:
        #            try:
        #                set_chi_dihedral (molecule  = molecule,
        #                                      resi  = key,
        #                                      bond  = bond,
        #                                      angle = fragment[key][bond])
        #            except:
        #                pass
        
        
                #except:
                #    print 'fail', bond
            #except KeyError as error:
            #    logger.debug(error.message)
            #    logger.debug("Target: " + "".join([
            #        three_to_one(molecule.residues[i].name)
            #            if molecule.residues[i].name != 'HIE' else 'H'
            #            for i in fragment
            #    ]))
            #    logger.debug("Fragment: " + "".join([
            #        three_to_one(a['NAME']) for a in fragment.values()
            #    ]))

def fragment_selection (molecule           = None  ,
                        previous_fragment  = None  ,
                        fragment_rate      = None  , 
                        random             = None  ):
    
    #print random
    
    fragment_acceptance = random.uniform(0, 1)
    FRAGMENTS = False
    PHIPSI    = False
    
    #----------------------------------------------#
    #                  FRAGMENTS                   #
    #----------------------------------------------#
    # (1) se o numero for menor ou igual a chance, entao um novo fragmento eh atribuido a estrutura
    if fragment_acceptance <= fragment_rate:
             
        # (2) sorteia uma posicao do alinhamento
        resi = random.randint(0, len(molecule.fragments)-1)

        # (3) enquando nao houver fragmentos para esta posicao (ou seja, lista vazia) sortear novas posicoes
        while molecule.fragments[resi] == []:
            resi = random.randint(0, len(molecule.fragments)-1) #(-1) -> nao pegar a ultima posicao
            
        # (4) sorteia um fragmento para a posicao selecionada 
        fragment_index = random.randint(0, len(molecule.fragments[resi])-1)
        fragment       = molecule.fragments[resi][fragment_index]
        
        if fragment != previous_fragment:
            return fragment, fragment_index
        
        else:
            return False
    else:
        return False



def insert_fragment_from_dic (molecule = None, fragment = {}, sidechain = False):
    """ Function doc """
    for i in fragment:
        phi_final_angle = set_phi_psi_dihedral( molecule=molecule, resi=i, bond='PHI',angle   = fragment[i]['PHI']  )
        psi_final_angle = set_phi_psi_dihedral( molecule=molecule, resi=i, bond='PSI',angle   = fragment[i]['PSI']  )
        ome_final_angle = set_phi_psi_dihedral( molecule=molecule, resi=i, bond='OMEGA',angle = fragment[i]['OMEGA'])
        print (i, phi_final_angle, psi_final_angle, ome_final_angle)

    for i in fragment:
        if sidechain:
            for chi in molecule.torsions[molecule.residues[i].name]:
                
                try:
                    angle = fragment[i][chi]
                    set_chi_dihedral (molecule=molecule, resi=i, bond=chi, angle = angle, log = False)
                    #fragment[i][chi] = computeCHI (molecule= molecule, resi=i, bond=chi)
                    #print i, residue.name, chi,  fragment[i][chi] 
                except:
                    pass
                     
                    #fragment[i][chi] = None
                    #print i, residue.name, chi,  fragment[i][chi] 
            





def metropolis_acceptance_test (energy = None        ,
                       previous_energy = None        ,
                           temperature = 273.15      ,
                                    Kb = 1           ,
                                random = None): #0.0019872041):
    """ Function doc """
    if energy < previous_energy:
        return True
    
    else:
        dE = (energy - previous_energy)
        
        Px = math.exp(-1 * dE / (Kb * temperature))
        X  = random.uniform(0, 1)
        #print('metropolis_acceptance_test', X, Px)
        return X <= Px



def rotate_backbone_attempt (molecule  = None,
                                 resi  = None,
                                 bond  = None,
                                theta  = None,
                      previous_energy  = None,
                          temperature  = None,
                                   pn  = None,
                                   Kb  = None,
                                random = None):

    rotate_backbone(molecule = molecule,
                        resi = resi    ,
                        bond = bond    ,
                       theta = theta   )

    energy = molecule.energy()
    #print ('energy rotate_backbone_attempt:',energy )
    if energy is not False:
        if metropolis_acceptance_test (energy  = energy          ,
                              previous_energy  = previous_energy ,
                                  temperature  = temperature     ,
                                           Kb  = Kb              , 
                                       random  =  random):
            
            #print ('metropolis_acceptance_test', energy)
            return energy

        else:
            
            return False
    else:
        #print ('metropolis_acceptance_test', energy)
        return False



def monte_carlo_cycle (molecule            = None                   ,
					   random              = None                   ,
					   
					   temperature          = 1000                   ,
					   #final_temperature   = False                  ,
					   #gamma               = False                  ,
					   
					   
					   Kb                  = 1                      , # 0.0083144621               ,
					   angle_range         = 60                     ,
					   fragment_rate       = 1.0                    , #between 0  and 1
					   fragment_sidechain  = False                  ,
					   PhiPsi_rate         = 0.0                    ,
					   attempt_per_residue = 5                      ,
					   					
					   
					   #simulated_annealing = None                   , # exp, linear
					   cycle_size          = 10000                   ,
					   #number_of_cycles    = 10                     ,
					   
					   log_frequence       = 10                     ,
					   trajectory          = 'MonteCarlo_trajectory',
					   pn                  = 1                      ,
					   log                 = True                   ,):


    if random == None:
        import random as random

    sum_energy = 0 
    attempted  = 0
    previous_energy = molecule.energy()
    previous_coordinates = molecule.get_coordinates_from_system()
    #print 'cycle:', cycle
    waste          = trajectory+'.waste.xyz'  
    #--------------------------------#
    #       attempted_accepted       #
    #--------------------------------#
    attempted_fragment = 0.0         #
    accepted_fragment  = 0.0         #
			             #
    attempted_phi = 0.0              #
    accepted_phi  = 0.0              #
			             #
    attempted_psi = 0.0              #
    accepted_psi  = 0.0              #
    #--------------------------------#
    frame      = 0
    fragment = None 
    previous_fragment    = None






    # LOGFILES - FRAGMENTS
    fragments_failed        = trajectory+'.fragments_failed.xyz'
    fragments_failed_text   = ''
    fragments_failed        = open(fragments_failed, 'a')
    
    fragments_selected      = trajectory+'.fragments_selected.xyz'
    fragments_selected_text = ''
    fragments_selected      = open(fragments_selected, 'a')
    
    fragments_rejected      = trajectory+'.fragments_rejected.xyz'
    fragments_rejected_text = ''
    fragments_rejected      = open(fragments_rejected, 'a')
    
    fragments_accepted      = trajectory+'.fragments_accepted.xyz'
    fragments_accepted_text = ''
    fragments_accepted      = open(fragments_accepted, 'a')










    for i in range(0, cycle_size):
        #nSteps +=1 
        
        # parametros controla a taxa de tentativas com que novos fragmentos sao testatos
        fragment_acceptance = random.uniform(0, 1)
        FRAGMENTS = False
        PHIPSI    = False


        if fragment_rate != 0.0:
            fragment =  False
            fragment_acceptance = random.uniform(0, 1)
            #----------------------------------------------#
            #                  FRAGMENTS                   #
            #----------------------------------------------#
            # (1) se o numero for menor ou igual a chance, entao um novo fragmento eh atribuido a estrutura
            if fragment_acceptance <= fragment_rate:
                     
                # (2) sorteia uma posicao do alinhamento
                resi = random.randint(0, len(molecule.fragments)-1)
                
                # (3) enquando nao houver fragmentos para esta posicao (ou seja, lista vazia) sortear novas posicoes
                while molecule.fragments[resi] == []:
                    resi = random.randint(0, len(molecule.fragments)-1) #(-1) -> nao pegar a ultima posicao
                    
                # (4) sorteia um fragmento para a posicao selecionada 
                fragment_index = random.randint(0, len(molecule.fragments[resi])-1)
                #print ()
                fragment       = molecule.fragments[resi][fragment_index]


                #---------------------------------------------------
                if fragment != previous_fragment:
                    pass
                    #return fragment, fragment_index
                else:
                    fragment =  False
                #---------------------------------------------------
                
                
                
                
                '''----------------------------------------------------'''
                if fragment:
                    fragtext  = ''
                    for k in range(len(molecule.residues)):
                        if k in fragment:
                            fragtext    += 'X'                    
                        else:
                            fragtext    += '-'
                    
                    #print '%4d %3d %s'%(resi, fragment_index,  fragtext) #print resi, fragment_index#, fragment 
                    
                    fragments_selected_text = '%4d %3d %s\n'%(resi, fragment_index,  fragtext)
                    fragments_selected.write(fragments_selected_text)
                
                else:
                    fragments_failed_text = '%4d %3d \n'%(resi, fragment_index)
                    fragments_failed.write(fragments_failed_text)
                    
                    #print resi,fragment_index,  'failed'
                '''----------------------------------------------------'''

        if fragment: 
            
            attempted_fragment += 1.0
            #------- associa o fragmento selecionado com a estrutura -------
            previous_fragment = fragment
            
            insert_fragment (molecule   = molecule,
                             fragment   = fragment,
                             sidechain  = fragment_sidechain)
            
            #---------------------------------------------------------------
            energy = molecule.energy()            
            
            if energy  or energy != None:

                if metropolis_acceptance_test (energy          = energy         ,
                                               previous_energy = previous_energy,
                                               temperature     = temperature    , 
                                               random          = random         ,
                                               Kb              = Kb             ):
                    frame += 1
                    save_XYZ_to_file (molecule, trajectory)
                    previous_energy      = energy
                    previous_coordinates = molecule.get_coordinates_from_system()
                    
                    accepted_fragment += 1
                    
                    FRAGMENTS = True

                    '''----------------------------------------------------'''
                    fragtext = ''
                    for k in range(len(molecule.residues)):
                        if k in fragment:
                            fragtext    += 'X'                    
                        else:
                            fragtext    += '-'
                    
                    print ('%4d %3d %s'%(resi, fragment_index,  fragtext)) #print resi, fragment_index#, fragment 
                    fragments_accepted_text = '%4d %3d %s\n'%(resi, fragment_index,  fragtext)
                    fragments_accepted.write(fragments_accepted_text)
                    '''----------------------------------------------------'''


                
                
                
                else:
                    save_XYZ_to_file (molecule, waste)                            # salva as coordenadas no descarte
                    molecule.import_coordinates_to_system (previous_coordinates)  # Restaura as coordenadas entriores 
                    #print 'fragment: ',fragment_index, energy, len(fragment), 'failed'
                    
                    
                    '''----------------------------------------------------'''
                    fragtext = ''

                    for k in range(len(molecule.residues)):
                        if k in fragment:
                            fragtext    += 'X'                    
                        else:
                            fragtext    += '-'

                    print ('%4d %3d %s  - REJECTED'%(resi, fragment_index,  fragtext)) #print resi, fragment_index#, fragment
                    fragments_rejected_text = '%4d %3d %s\n'%(resi, fragment_index,  fragtext)
                    fragments_rejected.write(fragments_rejected_text)
                    '''----------------------------------------------------'''

            else:
                save_XYZ_to_file (molecule, waste)                                # salva as coordenadas no descarte
                molecule.import_coordinates_to_system (previous_coordinates)      # Restaura as coordenadas entriores 

                '''----------------------------------------------------'''
                fragtext = ''
                for k in range(len(molecule.residues)):
                    if k in fragment:
                        fragtext    += 'X'                    
                    else:
                        fragtext    += '-'
                print ('%4d %3d %s  - REJECTED'%(resi, fragment_index,  fragtext)) #print resi, fragment_index#, fragment
                fragments_rejected_text = '%4d %3d %s\n'%(resi, fragment_index,  fragtext)
                fragments_rejected.write(fragments_rejected_text)
                '''----------------------------------------------------'''






        #----------------------------------------------#
        #               PHI/PSI sampling               #
        #----------------------------------------------#

        PhiPsi_acceptance = random.uniform(0, 1)

        if PhiPsi_acceptance <= PhiPsi_rate:
            resi = random.randint(0, len(molecule.residues)-1)

            if resi in molecule.fixed_residues:
                pass
            
            else:
                
                for attempt in range(0, attempt_per_residue ):
                    for bond in ['PSI','PHI']:

                        if bond == 'PSI':
                            attempted_phi += 1
                        if bond == 'PHI':
                            attempted_psi += 1

                        #attempted_psi += 1
                        theta  = random.uniform(-1 * angle_range, angle_range)
                        theta = math.radians(theta)

                        #------------------------------------------#
                        #         rotate_backbone_attempt          #
                        #------------------------------------------#
                        energy = rotate_backbone_attempt (molecule = molecule      ,
                                                            resi = resi            ,
                                                            bond = bond            ,
                                                           theta = theta           ,
                                                  previous_energy = previous_energy,
                                                     temperature = temperature     ,
                                                              pn = pn              ,
                                                          random = random          ,
                                                              Kb = Kb              )
                        #energy = False
                        #print(previous_energy,energy )
                        if energy is not False:
                            frame += 1
                            save_XYZ_to_file (molecule, trajectory+".xyz")
                            previous_energy      = energy
                            previous_coordinates = molecule.get_coordinates_from_system()

                            if bond == 'PSI':
                                accepted_psi += 1
                            if bond == 'PHI':
                                accepted_phi += 1
                            PHIPSI = True
                            #text.append("pn: {:<3d} step: {:5d} energy: {:<20.7f}rotate_backbone theta: {:<6.3f}\n".format(pn, i, energy , theta*57.324))
                            print ("step: {:5d}, resi: {:5d} energy: {:<20.7f}rotate_backbone theta: {:<6.3f}".format( i, resi, energy , theta*57.324))
                        else:
                            #salva as coordenadas que foram descartadas
                            #save_XYZ_to_file (molecule, waste)
                            
                            molecule.import_coordinates_to_system (previous_coordinates)


    if energy:
        if energy <= previous_energy:
            #frame +=1
            
            if FRAGMENTS:
                frag   = "fragment: {:<4d}  ".format(fragment_index)
            else:
                frag   = 'fragment: -     ';format('-')
            if PHIPSI:
                phipsi = "rotate_backbone: {:<6.3f}".format(theta*57.324)
            else:
                phipsi = 'rotate_backbone: None'
            text.append("pn: {:<3d} frame: {:5d} cycle: {:5d} step: {:5d} energy: {:<20.10f}".format(pn, frame, cycle,  nSteps, energy) + frag + phipsi+'\n')
            
            sum_energy        += energy
            attempted         += 1
    #---------------------------------------------#
    #                LOGFILEWRITE                 #
    #---------------------------------------------#
    logfile_counter += logfile_counter

    if logfile_counter >= log_frequence:
        logfile.writelines(text)
        logfile.close()
        text    = []
        logfile = open(logfilename, 'a')
        logfile_counter = 1



























