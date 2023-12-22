"""BulletinBoard python version Measuring passive, and Rheobase
with dynamic synaptic input"""

## --- modules imports

# import matplotlib  # to work on NSG
# matplotlib.use('Agg')  # to work on NSG. must be called the first thing
## Rheobase Measurements..........................................
# from Rheobase_NSG import Rheobase, Rheobase_binary
from mpi4py import MPI # must be called to initialize the mpi correctly
import time
import numpy as np
from neuron import h
h.load_file("stdrun.hoc")
from BBGlobals import *
from BBSysI_v3  import *


h.nrnmpi_init()

# templates that define cell classes..............................
h.load_file("FRcellTemplate.hoc") 		# FR cell template
h.load_file("STemp2.hoc")				# S cell template
h.load_file("FFTemp.hoc")				# FF cell template
# h.load_file("BBAxon.hoc")				# Axon Module

# insert Hoc functions to modify the cell template:
h.load_file("BBModifyCell.hoc")

# create ParallelContext instance here so it exists whenever we need it
pc = h.ParallelContext()


## functions to print computation time
def print_time(seconds):
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	print("\n## Runtime = %d:%02d:%02d\n" % (h, m, s))
##End of Function
###################################################################

def getType_Thresh(layer = "base"):

	# print("Layer = %s"%(layer))

	low_thresh  = np.array([S_LowThresh, R_LowThresh, F_LowThresh])
	high_thresh = np.array([S_HighThresh, R_HighThresh, F_HighThresh])

	low_vary 	= np.array([S_LowVary, R_LowVary, F_LowVary])
	high_vary 	= np.array([S_HighVary, R_HighVary, F_HighVary])

	if(layer == "low"):
		layerCode = 10
		return (low_thresh-low_vary, high_thresh-low_vary,layerCode)
	elif(layer == "high"):
		layerCode = 100
		return (low_thresh+high_vary, high_thresh+high_vary,layerCode)
	elif(layer == 'base'):
		layerCode = 1
		return (low_thresh, high_thresh,layerCode)
	else:
		print("invalid option for Threshold layer")
	##End of Function
#########################################################################



def SingleRun(Cell_type, Cell_num, layer):

	# print("start with CellType %2d , CellNumber %2d" %(Cell_type,Cell_num))

	Low_thresh, High_thresh, layerCode = getType_Thresh(layer)

	returnvec 	= h.Vector()
	tempvec 		= h.Vector()

	# // This block sets up the cell based on the input arg.
	if(Cell_type ==1):
		tempcell = h.S_Cell()

		# set the cell properties
		h.SpreadRin(tempcell, S_LowSoma, S_HighSoma, S_LowDen, S_HighDen, S_numcells, Cell_num)
		h.SpreadThresh(tempcell, Low_thresh[0], High_thresh[0], S_numcells, Cell_num)
		h.makeLLVA(S_CaProxDend, S_CaDistDend, S_GCaDend, tempcell, S_CaThetaLow, S_CaThetaHigh, S_numcells, Cell_num, NeuroMod)

		
		## create synapses and turn them on
		CreateSysISyn(tempcell, S_SysISyn_tau1, S_SysISyn_tau2, S_SysISyn_Erev, S_TotalSysIGsyn)
		SysI_ON(tempcell)
		StoreInitialSysIWeights(tempcell)


	elif(Cell_type ==2):
		tempcell = h.FR_Cell()

		# set the cell properties
		h.SpreadRin(tempcell, R_LowSoma, R_HighSoma, R_LowDen, R_HighDen, R_numcells, Cell_num)
		h.SpreadThresh(tempcell, Low_thresh[1], High_thresh[1], R_numcells, Cell_num)
		h.makeLLVA(R_CaProxDend, R_CaDistDend, R_GCaDend, tempcell, R_CaThetaLow, R_CaThetaHigh, R_numcells, Cell_num, NeuroMod)

		## create synapses and turn them on
		CreateSysISyn(tempcell, R_SysISyn_tau1, R_SysISyn_tau2, R_SysISyn_Erev, R_TotalSysIGsyn)
		SysI_ON(tempcell)
		StoreInitialSysIWeights(tempcell)

	else:
		tempcell = h.FF_Cell()

		# set the cell properties
		h.SpreadRin(tempcell, F_LowSoma, F_HighSoma, F_LowDen, F_HighDen, F_numcells, Cell_num)
		h.SpreadThresh(tempcell, Low_thresh[2], High_thresh[2], F_numcells, Cell_num)
		h.makeLLVA(F_CaProxDend, F_CaDistDend, F_GCaDend, tempcell, F_CaThetaLow, F_CaThetaHigh, F_numcells, Cell_num, NeuroMod)

		## create synapses and turn them on
		CreateSysISyn(tempcell, F_SysISyn_tau1, F_SysISyn_tau2, F_SysISyn_Erev, F_TotalSysIGsyn)
		SysI_ON(tempcell)
		StoreInitialSysIWeights(tempcell)


	## the following line for synaptic input update
	tempfih2 = h.FInitializeHandler((fi_update,(tempcell)))


	# #This block records spike times to a vector for use.
	nc = h.NetCon(tempcell.soma(0.5)._ref_v,None,sec = tempcell.soma)
	nc.threshold 	= -35
	tempvec		= h.Vector()
	nc.record(tempvec)
	#---------------------------------------------------------------

	# tempic = h.TriangleIClamp(0.5)
	# tempic.delay = 500
	# tempic.durA = 5000
	# tempic.durB = 5000
	# tempic.AMPstart = 0
	# tempic.AMPmax = 20
	# tempic.AMPend = 20

	# tempic2 = h.IClamp(tempcell.soma(0.5)) # Used an IClamp as an easy example of what we can do with this.
	# tempic2.delay 	= 100
	# tempic2.dur 	= 100
	# tempic2.amp 	= 20


	h.run()

	returnvec.append(Cell_type * layerCode)
	returnvec.append(Cell_num)
	returnvec.append(tempvec)


	myid	= h.hoc_ac_
	return (myid, Cell_type * layerCode, Cell_num, returnvec)
##End of Function
#########################################################################

DataStore	=	[None]

# //Necessary to make this whole process work.
pc.runworker()

# For simulation timing
start_time	= time.time()
print("Simulation time is set to: %g ms" %h.tstop)

def run_pool():
	for run_id in range(S_numcells):
		pc.submit(SingleRun, 1, run_id,"base") 
		pc.submit(SingleRun, 1, run_id, "low") 
		pc.submit(SingleRun, 1, run_id, "high") 
	for run_id in range(R_numcells):
		pc.submit(SingleRun, 2, run_id,"base") 
		pc.submit(SingleRun, 2, run_id, "low") 
		pc.submit(SingleRun, 2, run_id, "high") 
	for run_id in range(F_numcells):
		pc.submit(SingleRun, 3, run_id,"base") 
		pc.submit(SingleRun, 3, run_id, "low") 
		pc.submit(SingleRun, 3, run_id, "high") 

	while pc.working():
		job_id,Cell_type,cell_num , returnvec 	= pc.pyret()
		print("Finished job %2d , CellType %2d , CellNumber %2d" %(job_id ,Cell_type,cell_num))
		DataStore.append(returnvec)
##-------------------------End of function--------------------

run_pool()  # run the bulletin board
print_time(time.time()-start_time)

pc.done()
# all simulations have been completed
# and the workers have been released
# but the boss still has things to do



##------------------SAVING FUNCTIONS-------------------------//

def save_results2(DataStore):
	"""This function is used to save the firing time of the cells
	in files named after the cells
	M.H Feb 2018"""
	print("\n||-----SAVING by Cell:Type_Name-----||\n")
	for i in range(len(DataStore)-1):
		tempV 		= np.array(DataStore[i+1])
		filename 	= "OutFile%d_%d.txt" %(tempV[0],tempV[1])
		np.savetxt(filename,tempV,fmt='%f')
		print("%s have been saved" % filename)

	print("\n||-----Data Saved-----||\n")
##End of Function
############################################################################

def saveTo_dat(DataStore):
	"""This function save to the .dat format (of NEURON)
	consider using tofile and fromfile of python"""
	print("\n||-----SAVING--TO .dat format---||\n")
	for i in range(len(DataStore)-1):
		filename 	= "OutVec%d.dat" %i
		fileObj		= h.File()
		fileObj.wopen(filename)
		DataStore[i+1].vwrite(fileObj) 
		fileObj.close()
		print("%s have been saved" % filename)
	print("\n||-----Data Saved-----||\n")
##End of Function
############################################################################

## calling Saving functions
save_results2(DataStore)
saveTo_dat(DataStore) # needed to generate force from the .dat files

print("Cells are coded, low Rheobase layer x10, High Rheobase layer x100")
