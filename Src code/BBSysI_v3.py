"""
## Updated March 2023 Mohamed.Mousa
## Mohamed.mousa@wright.edu


Dec 2021
# Trpazoidal wave form is now enabled, use the 'hold_Time' variable in BBGlobals.py to set the hold time.
# """

from neuron import h
from BBGlobals import *
# from matplotlib import pyplot as plt  # cause problems on neuroscience gateway
import numpy as np


## Global variables....................
SYSIFLAG 	= 0
NewFactor 	= 0  # modified to zero ,changed each time by PERCENT_CHANGE till it reach 1
STEPNUMBER 	= 0	 #Parameters for adjusting synaptic NetCon weights during runtime.
LASTTIME 	= STARTTIME + INTERVAL*(HOWMANY_UP-1) + INTERVAL*(HOWMANY_DOWN) + hold_Time

# print("Syn sys V3 loaded")



## CreateSysISyn
def CreateSysISyn(cell,Tau1,Tau2,Erev,Total_Gsyn):
	"""
	Process to create system I synapses
	Arg 1: Cell Object to apply action on.
	Arg 2: Value of tau1 for synapses
	Arg 3: Value of tau2 for synapses
	Arg 4: Value of Reversal Potential for synapses
	Arg 5: Total Conductance for Sys I for this cell class."""
	# print "synapses Created"
	cell.TotalSysISynNumber = 0
	cell.TotalArea 			= 0
	cell.SconductanceVec_SysI = h.Vector()			#

	# Compute total area
	for sec in cell.den:
		for seg in sec:     # for sec.allseg() for acessing the X=0 and X =1 seg which area is zero
			cell.TotalArea += seg.area()					# Add up membrane area of dendritic compartments


	# Compute total number of synapses to be created
	for sec in cell.den:
		cell.TotalSysISynNumber += sec.nseg			# Add one to the counter for every dendritic segment.

	## for debugging purposes :-------
	# print "Cell Total Synapses = %d" % cell.TotalSysISynNumber 										# At this point, we're basically printing the number of segmets for all dendrites.
	# print "Cell Total Area = %.2f" % cell.TotalArea

	# ssA = "%s.SysISyn[%d]" %(cell,cell.TotalSysISynNumber)  		# define an obj variable
	# exec(ssA)							.							# execute fuction bounces its argument to the top of the stack.
	# ssA = "%s.SysI_NC[%d]" %(cell,cell.TotalSysISynNumber)  		# define an obj variable for net con
	# exec(ssA)

	#Blocks for actual synaptic creation below.
	cell.TotalSysISynNumber = 0											# Reset counter of total synapses
	for sec in cell.den:												# Iterate over the entire list of cell.den
		cell.SectionArea = 0											# Reset variable
		cell.Gsynapse = 0												# Reset variable
		
		for seg in sec:
			cell.SectionArea += seg.area()				# Compute the area of the currently accessed section

		cell.NumSysISynapsesPerSection = sec.nseg		# Number of synapses to be placed is equal to the number of sections in that segment

		for seg in sec:
			cell_secName = sec.name().split(".")
			ssA = "cell.SysISyn[{0:d}] = h.Exp2Syn({1}({2}))".format(int(cell.TotalSysISynNumber),"cell."+cell_secName[1],seg.x)		# create one synapse at the middle of each section
			# print ssA
			exec(ssA)
			if(h.n3d(sec) > 0):												# Confirm that 3D data is present
				if (h.x3d(0,sec) >= 0):											# Check if synapse is on Rostral or Caudal side
					cell.Rostral_Syn_sysI.append(cell.TotalSysISynNumber)	# Add synapse to Rostral synapse list
				else:
					cell.Caudal_Syn_sysI.append(cell.TotalSysISynNumber)	# Add synapse to Caudal synapse list

				if (h.y3d(0,sec) >= 0):										# Check if synapse is on Dorsal or Ventral side
					cell.Dorsal_Syn_sysI.append(cell.TotalSysISynNumber)	# Add synapse to Dorsal synapse list
				else:
					cell.Ventral_Syn_sysI.append(cell.TotalSysISynNumber)	# Add synapse to Ventral synapse list

				if (h.z3d(0,sec) >= 0):										# Check if synapse is on Medial or Lateral side
					cell.Medial_Syn_sysI.append(cell.TotalSysISynNumber)	# Add synapse to Medial synapse list
				else :
					cell.Lateral_Syn_sysI.append(cell.TotalSysISynNumber)	# Add synapse to Lateral synapse list


			ssA = "cell.SysISyn[%d].tau1 = %f" %(cell.TotalSysISynNumber, Tau1)  # [ms], an alpha conductance (as tau1=tau2) that peaks after 0.2 ms
			exec(ssA)																# move to the top of the stack.
			ssA = "cell.SysISyn[%d].tau2 = %f" %(cell.TotalSysISynNumber, Tau2)  # [ms], an alpha conductance (as tau1=tau2) that peaks after 0.2 ms
			exec(ssA)
			ssA = "cell.SysISyn[%d].e = %f" %(cell.TotalSysISynNumber,Erev) 		# [mV], the reversal potential of the synapses
			exec(ssA)

			cell.Gsynapse = (Total_Gsyn*cell.SectionArea)/(cell.TotalArea*cell.NumSysISynapsesPerSection)		# Compute the conductance of synapses in that section
			cell.SconductanceVec_SysI.append(cell.Gsynapse)								# Store synapse conductance to a vector

			cell.TotalSysISynNumber +=1													# update counter of total synapses

	global SYSIFLAG
	SYSIFLAG = 1																# Turn flag ON
# End of Function
#############################################################################################################


# Arg 1: Cell Object.
# Maybe set up additional args to allow various parameters to be adjusted more freely?
def SysI_ON(cell):
	# print "synapses ON"
	if (SYSIFLAG == 1):			# Check if the Ia-Afferents synapses have been created or not

		# Create 4 sources of stimulation to each group of the Ia-Afferents Synapses
		cell.NS_A_sysI = h.NetStim(0.5)   						# creates the NetStim that will drive Group A
		cell.NS_A_sysI.interval = 1000.0/SysIVibrationFreq					# Periodic time [ms] --dividing by [Hz] yields [s]
		cell.NS_A_sysI.number = SysITimeofVibration*SysIVibrationFreq/1000.0		# time of stimulation

		cell.NS_B_sysI = h.NetStim(0.5)  		 					# creates the NetStim that will drive Group B
		cell.NS_B_sysI.interval = 1000.0/SysIVibrationFreq					# Periodic time
		cell.NS_B_sysI.number = SysITimeofVibration*SysIVibrationFreq/1000.0		# time of stimulation

		cell.NS_C_sysI = h.NetStim(0.5)   						# creates the NetStim that will drive Group C
		cell.NS_C_sysI.interval = 1000.0/SysIVibrationFreq					# Periodic time
		cell.NS_C_sysI.number = SysITimeofVibration*SysIVibrationFreq/1000.0		# time of stimulation

		cell.NS_D_sysI = h.NetStim(0.5)   						# creates the NetStim that will drive Group D
		cell.NS_D_sysI.interval = 1000.0/SysIVibrationFreq					# Periodic time
		cell.NS_D_sysI.number = SysITimeofVibration*SysIVibrationFreq/1000.0		# time of stimulation

		# Determine the pattern of activation
		if (SysISynActivation == 0):								# All-group activation pattern
			cell.NS_A_sysI.start = STARTTIME							# Gr A: start at time = 0
			cell.NS_B_sysI.start = STARTTIME							# Gr B: start at time = 0
			cell.NS_C_sysI.start = STARTTIME							# Gr C: start at time = 0
			cell.NS_D_sysI.start = STARTTIME							# Gr D: start at time = 0
		elif(SysISynActivation == 2):							# 2-group activation pattern
			cell.NS_A_sysI.start = STARTTIME							# Gr A: start at time = 0
			cell.NS_B_sysI.start = STARTTIME + 1000.0/SysIVibrationFreq/2.0	# Gr B: Phase shift to start at 50% of cycle
			cell.NS_C_sysI.start = STARTTIME							# Gr C: start at time = 0
			cell.NS_D_sysI.start = STARTTIME + 1000.0/SysIVibrationFreq/2.0	# Gr D: Phase shift to start at 50% of cycle
		elif(SysISynActivation == 4):							# 4-group activation pattern
			cell.NS_A_sysI.start = STARTTIME									# Gr A: start at time = 0
			cell.NS_B_sysI.start = STARTTIME + 1000.0/SysIVibrationFreq/4.0			# Gr B: Phase shift to start at 25% of cycle
			cell.NS_C_sysI.start = STARTTIME + 1000.0/SysIVibrationFreq/2.0			# Gr C: Phase shift to start at 50% of cycle
			cell.NS_D_sysI.start = STARTTIME + 1000.0/SysIVibrationFreq/4.0*3.0		# Gr D: Phase shift to start at 75% of cycle

		# Set synapse parameters
		iii = 0
		while(iii<cell.TotalSysISynNumber):					# This functions as a WHILE loop
			# Link synapses to net con
			if(iii<cell.TotalSysISynNumber):
				ssA  = "%s.SysI_NC[%d] = h.NetCon(%s.NS_A_sysI, %s.SysISyn[%d])" %("cell",iii,"cell","cell",iii)		# link 1st synapse to net con A
				exec(ssA)													# top of stack
				ssA  = "%s.SysI_NC[%d].weight[0] = %f" %("cell",iii,cell.SconductanceVec_SysI.x[iii])  			# read the synapse conductance from vector
				exec(ssA)																		# top of stack
				ssA  = "%s.SysI_NC[%d].delay = %f" %("cell",iii,SysISyn_Delay)
				exec(ssA)																		# top of stack
				iii+=1


			if(iii<cell.TotalSysISynNumber):
				ssA  = "%s.SysI_NC[%d] = h.NetCon(%s.NS_B_sysI, %s.SysISyn[%d])" %("cell",iii,"cell","cell",iii)	# link 2nd synapse to net con B
				exec(ssA)
				ssA  = "%s.SysI_NC[%d].weight[0] = %f" %("cell",iii, cell.SconductanceVec_SysI.x[iii])  		# read the synapse conductance from vector
				exec(ssA)
				ssA  = "%s.SysI_NC[%d].delay = %f" %("cell",iii,SysISyn_Delay)
				exec(ssA)
				iii+=1


			if(iii<cell.TotalSysISynNumber):
				ssA  = "%s.SysI_NC[%d] = h.NetCon(%s.NS_C_sysI, %s.SysISyn[%d])" %("cell",iii,"cell","cell",iii)	# link 3rd synapse to net con C
				exec(ssA)
				ssA  = "%s.SysI_NC[%d].weight[0] = %f" %("cell",iii, cell.SconductanceVec_SysI.x[iii])  		# read the synapse conductance from vector
				exec(ssA)
				ssA  = "%s.SysI_NC[%d].delay = %f" %("cell",iii,SysISyn_Delay)
				exec(ssA)
				iii+=1


			if(iii<cell.TotalSysISynNumber):
				ssA  = "%s.SysI_NC[%d] = h.NetCon(%s.NS_D_sysI, %s.SysISyn[%d])" %("cell",iii,"cell","cell",iii)	# link 4th synapse to net con D
				exec(ssA)
				ssA  = "%s.SysI_NC[%d].weight[0] = %f" %("cell",iii, cell.SconductanceVec_SysI.x[iii])  		# read the synapse conductance from vector
				exec(ssA)
				sssA  = "%s.SysI_NC[%d].delay = %f" %("cell",iii,SysISyn_Delay)
				exec(ssA)
				iii+=1
		# print "SYS Synapses ON"

# End of Procedure
#####################################################################################################


# This procedure store the initial weight of all synapses. This initial weight will be used to grade the synaptic activation
# If we need to reset synapses or modulate their weights, we do this here
# modified by M.H Feb 2018 to apply zero weights at the beginning
# Arg 1: Cell object
def StoreInitialSysIWeights(cell, reset = True):
	"""This function store the initial weights of the synaptic input, This initial weight will be used to grade the synaptic activation
	then initialize the weights to zero
	:Param cell : cell object"""
	# print("Intial weights Stored")
	global NewFactor
	NewFactor = 0  # modified to zero

	cell.InitialSysIWeights.resize(cell.TotalSysISynNumber)

	for i in range(int(cell.TotalSysISynNumber)):
		ssA =  "%s.InitialSysIWeights.x[%d] = %s.SysI_NC[%d].weight[0]" %("cell", i,"cell", i)  	# Save the initial synaptic weight of all synapses
		# print(ssA)
		exec(ssA)

	if(reset):
		## Reset Synaptic input weights to zero , to start the simulation with zero effet.
		for ii in range(int(cell.TotalSysISynNumber)):
			ssA = "%s.SysI_NC[%d].weight[0] = %f*%f"  %("cell", ii, cell.InitialSysIWeights.x[ii], NewFactor)  			# Update the synaptic weight of that synapse
			exec(ssA)

	# reset the CVODE integrator - a required step for variable time step simulations
	if (h.cvode_active()):
		h.cvode.re_init()
	# print("Initial Sys I weights stored.")
# End of Procedure
#####################################################################################################


# Update the synaptic weights during simulation time to grade the synaptic activation
# Arg 1: Cell Object
# May update this to allow better control of steps...?
def UpdateSysIWeights(cell):
	global  NewFactor , STEPNUMBER

	# Update step number counter-------->
	if ( (STEPNUMBER < HOWMANY_UP) and (h.t <= (STARTTIME + 0.5 + INTERVAL*(HOWMANY_UP - 1))) ): ## add 0.5 to solve for time roundoff ERROR
		STEPNUMBER += 1
		NewFactor += PERCENT_CHANGE_UP				# Newfactor is increasing during upward ramp
		
	elif( h.t > (STARTTIME + 0.5 + INTERVAL*(HOWMANY_UP - 1) + hold_Time )):
		STEPNUMBER += 1
		NewFactor -= PERCENT_CHANGE_DN				# Newfactor is increasing during upward ramp
		#comment above line for ramp to hold

	## print next synaptic weights Parameters
	# print("NewFactor %f , PERCENT_CHANGE = %f , STEPNUMBER %d ,at time %0.0f" %(NewFactor,PERCENT_CHANGE_UP,STEPNUMBER,h.t))


	# Update Synaptic Weights
	for ii in range(int(cell.TotalSysISynNumber)):
		ssA = "%s.SysI_NC[%d].weight[0] = %f*%f"  %("cell", ii, cell.InitialSysIWeights.x[ii], NewFactor)  			# Update the synaptic weight of that synapse
		exec(ssA)

	#print "New NC weight is %f" %(cell.SysI_NC[0].weight[0])
	# reset the CVODE integrator - a required step for variable time step simulations
	if (h.cvode_active()):
		h.cvode.re_init()

	# Determine when the next event will occur
	if (h.t+INTERVAL <= LASTTIME+0.5):
		h.cvode.event(h.t+INTERVAL,(UpdateSysIWeights,(cell)))
# End of Procedure
#################################################################################################




def fi_update(cell):
	"""This function is to create the event that will update the synaptic weights"""
	margin = 0 # int(500/SysIVibrationFreq)
	h.cvode.event(0+STARTTIME-margin,(UpdateSysIWeights,(cell)))
	print("Called fi_update")
	
##-------------End of Procedure----------------------------------------------


def restore_initialWeights(cell):
	for ii in range(int(cell.TotalSysISynNumber)):
		cell.SysI_NC[ii].weight[0] = cell.InitialSysIWeights.x[ii]
	print("Synaptic System I weights has been restored to their initial weights")
##-------------End of Procedure----------------------------------------------


def saveData(timeVec, CurrentVec, Cell_type, Cell_num):
	"""Saving the time-voltage-current file"""
	filename = "SynapticI_%d_%d.dat" %(Cell_type, Cell_num)
	np.savetxt(filename, list(zip(timeVec, CurrentVec)),
                fmt='%0.3f', header='Time Current')
				
	print(("#%s file has been saved" %(filename)))
#####-----------------End of Function------------------------------------------------



def AverageFilter(vector,window):
	"""Three stage Average Filter"""

	filtered = np.convolve(vector,np.ones(window)/window, 'same')
	filtered = np.convolve(filtered, np.ones(window)/window, 'same')
	filtered = np.convolve(filtered, np.ones(window)/window, 'same')

	return filtered
#####-----------------End of Function------------------------------------------------

def RMS(vec, window):
	"""Calculate RMS value for the vector"""
	RMSvec	= h.Vector()

	NofW	= int(len(vec)/window)

	for i in range(0,NofW-1,1):
		startind 		= i*window
		endind		= (i*window) + window
		square		= 0

		for j in range(startind,endind):
			square	+= vec[j]**2
		
		mean2 = square/float(window)
		RMSvec.append(h.sqrt(mean2))

	return RMSvec
#####-----------------End of Function------------------------------------------------

