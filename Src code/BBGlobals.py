# this file define the Global parameters for the MN pool simulation

from neuron import h

h.tstop 		= 35000						#Time for which to run each cell
h.v_init 		= -70						#Model initialization voltage
h.celsius 	= 36						# temperature at which simulations are performed

Model_Shape 	= 0

S_numcells 		= 39						#Number of cells of each type #S, FR, FF 
R_numcells 		= 15						# VL muscle cells			
F_numcells 		= 5						#
Tot_numcells 	= S_numcells+R_numcells+F_numcells

S_Counter 		= 0						#Track number of cells for each class.
R_Counter 		= 0
F_Counter 		= 0

# band distribution parameters in the LLVA channels on the S cells dendritic tree
S_CaProxDend 	= 0.42					#start point ( electronics distance )
S_CaDistDend 	= 0.9					#end point ( electronics distance )
S_GCaDend 	= 0.00014				#channels conductance
S_CaThetaLow 	= -44.55				#Activation curve range
S_CaThetaHigh 	= -45.35

R_CaProxDend 	= 0.4
R_CaDistDend 	= 0.9
R_GCaDend 	= 0.00014      # basic value : 0.00014, enhanced 0.000201
R_CaThetaLow 	= -43		#-40.68
R_CaThetaHigh 	= -44		#-46.54

F_CaProxDend 	= 0.4
F_CaDistDend 	= 0.9
F_GCaDend 	= 0.00014
F_CaThetaLow 	= -39.4
F_CaThetaHigh 	= -41.52

NeuroMod 		= 1.0	# Neuromodulatory State.  Default is 1.( change Llva conductance)

# firing properties of the S Cells - shifting the Na , K , etc channels activation curves to shift the firing threshold
S_LowThresh 	= -65				#Shift Threshold Parameters for S Cells
S_HighThresh 	= -61.5				# voltage thresh

R_LowThresh 	= -56
R_HighThresh 	= -60

F_LowThresh 	= -60
F_HighThresh 	= -59.7

S_LowVary 		= 1							#Values for - mV to generate cells to cover lower range of pool variability
R_LowVary 		= 2  					# need to vary the S, and FF like the FR, 2 up, and 2 down instead of 1 down, and 3 up
F_LowVary 		= 1

S_HighVary 		= 4							# + mV to cover higher range of pool variability
R_HighVary 		= 2
F_HighVary 		= 4

S_LowSoma 		= 230						# Numbers for membrane resistance values for cells.
S_HighSoma 		= 496.4
S_LowDen 			= S_LowSoma*30
S_HighDen 		= S_HighSoma*30

R_LowSoma 		= 77.875
R_HighSoma 		= 245.125
R_LowDen 			= R_LowSoma*50
R_HighDen 		= R_HighSoma*50

F_LowSoma 		= 22
F_HighSoma 		= 94
F_LowDen 			= F_LowSoma*250
F_HighDen 		= F_HighSoma*250

# Synapses system No #I parameters
# uniform input at 80% setting

S_SysISyn_Erev 	= 0					#reversal voltage ...
S_SysISyn_tau1 	= 0.2				# activation time constant
S_SysISyn_tau2 	= 0.2				# deactivation time constant
S_TotalSysIGsyn 	= 6.936  # 13.5nA  # must be float
									
R_SysISyn_Erev 	= 0					#...FR...
R_SysISyn_tau1 	= 0.2
R_SysISyn_tau2 	= 0.2				# deactivation
R_TotalSysIGsyn 	= 8.468  # 13.5nA
									
F_SysISyn_Erev 	= 0					#...and FF Cells.
F_SysISyn_tau1 	= 0.2
F_SysISyn_tau2 	= 0.2
F_TotalSysIGsyn 	= 8.92  # 13.5nA

HOWMANY_UP 		= 10000					#Maybe need to split these up into one group for each system...?
hold_Time			= 10000
HOWMANY_DOWN 		= 10000
STARTTIME 		= 100
INTERVAL 			= 1
PERCENT_CHANGE_UP 	= 1.0/(HOWMANY_UP)
PERCENT_CHANGE_DN 	= 1.0/(HOWMANY_DOWN) # made zero to ramp then hold


SysIVibrationFreq 	= 180	#180			# [Hz], activation frequency of synaptic input
#5000 Hz will produce a perfectly smooth input!
SysITimeofVibration = 1e9			# [ms], how long the synapses will be activated
SysISyn_Delay 		= 1	       		# [ms], the delay before activating the synapses
SysISynActivation 	= 4             	# Pattern of activation of the Ia-synapses, Choose one of the following:
								# 0: for All-group activation, all synapses are activated at the same time.
								# 2: for 2-group activation, synapses are divided into 2 groups activated with 50% phase shift
								# 4: for 4-group activation, synapses are divided into 4 groups activated with 25% phase shift

