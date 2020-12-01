
#import numpy as np
#import math

#from pprint import pprint
#from pepCore.Vectors import *



def monte_carlo_cycle (molecule            = None                   ,
					   random              = None                   ,
					   
					   initial_temperature = 1000                   ,
					   #final_temperature   = False                  ,
					   #gamma               = False                  ,
					   
					   
					   Kb                  = 1                      , # 0.0083144621               ,
					   angle_range         = 60                      ,
					   fragment_rate       = 0.0                    , #between 0  and 1
					   fragment_sidechain  = True                   ,
					   PhiPsi_rate         = 1.0                    ,
					   attempt_per_residue = 5                      ,
					   					
					   
					   #simulated_annealing = None                   , # exp, linear
					   cycle_size          = 10000                  ,
					   #number_of_cycles    = 10                     ,
					   
					   log_frequence       = 10                     ,
					   trajectory          = 'MonteCarlo_trajectory',
					   pn                  = 1                      ,
					   log                 = True                   ,):
    
    
    if random == None:
        import random as random
	
	sum_energy = 0 
	attempted  = 0
	
	#print 'cycle:', cycle
	
	
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
	
	
	
	for i in range(0, cycle_size):
		#nSteps +=1 

		# parametros controla a taxa de tentativas com que novos fragmentos sao testatos
		fragment_acceptance = random.uniform(0, 1)
		FRAGMENTS = False
		PHIPSI    = False
		
		
		#----------------------------------------------#
		#                 FRAGMENTS                    #
		#----------------------------------------------#
		
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
			
			#try:
			energy = molecule.energy(pn = pn)            
			#except:
			#    energy = False
			#print energy
			
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
					#text.append("pn: {:<3d} step: {:5d} energy: {:<20.7f}fragment: {:<3d}\n".format(pn, i, energy , fragment_index))
					#print "pn: {:<3d} step: {:5d} energy: {:<20.7f}fragment: {:<3d}".format(pn, i, energy , fragment_index)
				
				
				
				
					'''----------------------------------------------------'''
					fragtext = ''
					for k in range(len(molecule.residues)):
						if k in fragment:
							fragtext    += 'X'                    
						else:
							fragtext    += '-'
					
					#print '%4d %3d %s'%(resi, fragment_index,  fragtext) #print resi, fragment_index#, fragment 
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

					#print '%4d %3d %s  - REJECTED'%(resi, fragment_index,  fragtext) #print resi, fragment_index#, fragment
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
				#print '%4d %3d %s  - REJECTED'%(resi, fragment_index,  fragtext) #print resi, fragment_index#, fragment
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
						
						if energy:
							frame += 1
							save_XYZ_to_file (molecule, trajectory)
							previous_energy      = energy
							previous_coordinates = molecule.get_coordinates_from_system()

							if bond == 'PSI':
								accepted_psi += 1
							if bond == 'PHI':
								accepted_phi += 1
							PHIPSI = True
							#text.append("pn: {:<3d} step: {:5d} energy: {:<20.7f}rotate_backbone theta: {:<6.3f}\n".format(pn, i, energy , theta*57.324))
							#print "pn: {:<3d} step: {:5d} energy: {:<20.7f}rotate_backbone theta: {:<6.3f}".format(pn, i, energy , theta*57.324)
						else:
							save_XYZ_to_file (molecule, waste)
							molecule.import_coordinates_to_system (previous_coordinates)

		
		#if side_chain == True:
		#    monte_carlo_side_chain (molecule  = molecule        ,
		#                           rotamers   = rotamers        ,
		#                           random     = random          ,
		#                           initial_T  = temperature     ,
		#                           final_T    = temperature     ,
		#                           angle_range= angle_range     ,
		#                           nSteps     = side_chain_steps,
		#                           trajectory = trajectory      )
		
		
		#        got changes in energy?
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
		
		
   
