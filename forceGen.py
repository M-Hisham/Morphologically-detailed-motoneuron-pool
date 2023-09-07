"""
* Goal: this file generate force based on spikes timing for motor pool of 3 cell types (S, FR, FF)
* this script assign different force based on the cell type, as well as the recruitment order ( from cell name not spike time).
* can work for homogenous model or hetrogenous model, and number of cells do not affect the force limits of type
* Time: takes around 2 min for ~200 cell, with t.stop of 1000
* Time: takes around 20 min for ~200 cell, with t.stop of 20100


# Global Parameters for Fuglevand, et al, 1993 style force and EMG model
# ForceFac Developed from combination of our data in relation to Zengel 1985, and Krutki 2006
# Globally from the largest FF cell to the smallest S cell:
RT = 5.5 and TL = 110 taken from Burke, Levine, Zajac, Tsairis, and Engel, 1971

Developed initially by John Allen 2017, and rewritten in python and updated to
1) Differentiat force range based on the cell type
2) parallel vector computing ( on the time vector) instead of for loop to make it faster.
3) saves files based on celltype and cell number
4) print number of cell that did not fire,produce small files instead of full size vector filled with zeros

by Mohamed Hisham
Last updated: May 17, 2022.
Mohamed.Mousa@wright.edu.

"""
## --- modules imports
import numpy as np
from mpi4py import MPI # must be called to initialize the mpi correctly
from neuron import h
h.load_file("stdrun.hoc")
import time
import BBGlobals
import sys


## Global settings.>>>>
ForceCal       = 9.8 #Multiplicative Factor to go from mass units to actual force units. (g to milli Newton)
ForceFac       = ForceCal*np.array([1.2, 4.5, 30]) ## calibrated to young cat (gravitation Acce(9.8 m/s^2) * mass(g) = Force )
## Ranges are used from cat. scalling from human does not change it.>>>>
RP             = [10.5,12.22 ,4.33] # Range of twitch Forces [ S, FR, FF].  How many times the strongest twitch is the bigger.
RT             = [110/58 ,55/30, 47/20] # Range of Twitch Times [ S, FR, FF].  How many times weakest twitch is longer than the strongest. the weakest twitch duration is for strongest contraction.
## Calibrated to old cat.>>>>
TL             = [110, 55, 47] # Single longest twitch time [ S, FR, FF]. [ms]
Tot_numcells   = 3*np.array([39, 15, 5])   # total number of cells [ S, FR, FF]
h.tstop        = BBGlobals.h.tstop
CellnotFired   = 0

print("Hello to force Generator V2.5 /young Cat/")
print("tstop is",h.tstop)

# create ParallelContext instance here so it exists whenever we need it
pc = h.ParallelContext()

def print_time(seconds):
     m, s = divmod(seconds, 60)
     h, m = divmod(m, 60)
     print("\n## Runtime = %d:%02d:%02d\n" % (h, m, s))
##End of Function
###################################################################

def AssignTwitchForceP(RP, cell_num , CT ):
     # RP: Range of Twitch Forces(1: RP)[a.u.]
     # iFire: Index of Cell in Order of 1st to Fire(0: CELLS.count-1)
     # CT     : Cell type (0 or 1 or 2)
	# Equation for b, n = Number of Cells in Pool
	#b = log(RP)/Tot_numcells
	# Equation for P, i = ith Cell in Pool(i=1: n)
	# P = exp(b.*i)
     #return exp(b*($2 + 1)), then after simplification
     return RP[CT]**((cell_num)/Tot_numcells[CT])
##------end of function----------------------##          

def AssignTwitchForceT(RP, RT, TL, T_force, CT):
     """Return the twitch contraction time"""
     # $1 RP: Range of Twitch Forces(1: RP)[a.u.]
     # $2 RT: Range of Contraction Times(TL/RT: TL)[a.u.]
     # $3 TL: Longest Duration Contraction Time[ms]
     # $4 P: Twitch Force[a.u.]
     # CT  : Cell type (0 or 1 or 2)

     # Equation for c
     # c = log(RP)/log(RT)
     c = np.log(RP[CT])/np.log(RT[CT])
     # Equation for T, Given in [ms])
     # T = TL*((1/P) ^ (1/c))
     return TL[CT]*((1 / T_force) ** (1/c))  # Eq 14 in Fuglevand paper 
##------end of function----------------------##

def ForceModel2(isi, spiketime, T_Force, cont_T, time, CT):
     """ 
     This function calulate the gain between the force and stimulus rate.
     
     operate vectorwise, paster by O(n)
	# $1 isi: Inter-Spike Interval[ms]
	# $2 spiketime: Time of Current Spike[ms]
	# $3 T_Force: Twitch Force[a.u.]
	# $4 cont_T: Contraction Time[ms]
	# $5 time: Simulation Time[ms]
     # CT     : Cell type (0 or 1 or 2)
     """

     t = time - spiketime ## array

	# if T/isi <= 0.4 , gain = 1
     if (cont_T / isi <= 0.4):
          gij = 1
     else:
          gij0dot4 = (1-np.exp(-2*(0.4 ** 3)))/(0.4)
            # gij = (1-exp(-2*((T/isi) ^ 3)))/(T/isi)/gij0dot4
          gij = (1-np.exp(-2*((cont_T / isi) ** 3))) / (cont_T / isi)/gij0dot4

	# twitch = gij*((P*(time-spiketime)/T)*exp(1-((time-spiketime)/T)))
     return gij*((T_Force*t/cont_T)*np.exp(1-(t / cont_T))) *ForceFac[CT]
##------end of function----------------------##


def celltype(Cell_type_lf, cell_num):
     """Return cell type-1 = 0,1,2  and layer as 0, 1, 2
     cell layer of 0: low threshold
     cell layer of 1: normal threshold
     cell layer of 2: high threshold
     
     cell type of 0: S
     cell type of 1: FR
     cell type of 2: FF

     cell number fix, is the cell number from 0 to 3*[51, 9, 3]
     as the original cell_num  does only extend to 1/3 of these
     """
     layerMap = [1, 0, 2]


     cell_Layer     = int(np.log10(Cell_type_lf))
     cell_type      = int(Cell_type_lf/(10**cell_Layer) - 1  )
     cell_Layer     = int(layerMap[cell_Layer])

     # print("cell type = ", cell_type)
     # print("cell layer = ", cell_Layer)


     cell_num_fix   = cell_num + (cell_Layer*Tot_numcells[cell_type]/3)

     return cell_Layer, cell_type , int(cell_num_fix)
##------end of function----------------------##


def ForceGen(SpikeVec, Cell_type, cell_num):

     """
     ## Arg 1: Individual Cell firing time vector
     ## Arg 2: Cell number in pool based on order of first to fire.
     """
     Fvec_size      = int(h.tstop/h.dt)+1
     # print("Time Vector size = ", Fvec_size)
     # tempTwitch     = np.zeros(Fvec_size)
     Twitch         = np.zeros(Fvec_size)
     returnvec      = h.Vector()
     SimTime        = h.Vector()
     global CellnotFired

     ##fill simTime vector
     SimTime = np.array(SimTime.indgen(0, h.tstop, h.dt))

     # print("Cell type and num  = ", Cell_type,"   ", cell_num)

     _, MNtype,cell_numU = celltype(Cell_type, cell_num)

     # print("Cell num updated = ", cell_numU)

     ## calculate Twitch force, and twitch duration
     tempforce = AssignTwitchForceP(RP, cell_numU, MNtype)
     temptime  = AssignTwitchForceT(RP, RT, TL, tempforce, MNtype)


     isi           = (1/0.4)*temptime
     ## start with 2, as the first two elements of the vector are indentifiers
     for j in range(2,len(SpikeVec)):

          if(j>2):
               isi = SpikeVec[j]-SpikeVec[j-1]

          ## This Loop Iterates Over the Simulation Time Steps (by Index)
          ndx = int(SpikeVec[j]/h.dt)
          tempTwitch = np.zeros(Fvec_size) ## reset and recreate the vector as empty
          ### syntax: ForceModel(isi, spiketime, P, T, time)
          tempTwitch[ndx:] = ForceModel2(
              isi, SpikeVec[j], tempforce, temptime, SimTime[ndx:], MNtype)

          Twitch += tempTwitch  # append the temp vector to the main twitch.

     if(len(SpikeVec) <=2):
          print("cell  %2d:%2d did not fire" %(SpikeVec[0],SpikeVec[1]))
          Twitch         = (-1)*np.ones(2)
          CellnotFired   += 1


     ## return force calculations
     returnvec.append(SpikeVec[0])
     returnvec.append(SpikeVec[1])
     returnvec = np.append(returnvec,Twitch)

     myid = h.hoc_ac_
     return (myid, Cell_type, cell_num, returnvec)
##------end of function----------------------##


def BinarySpikesReader(fileIndex = -1, debug = False):


     DataStore = []

     if(fileIndex == -1 ):
          filelist  = range(np.sum(Tot_numcells))
     else:
         filelist = fileIndex

     for i in filelist:
          inVec     = h.Vector()
          inFile    = h.File()
          filename  = "OutVec%d.dat" %(i)
          inFile.ropen(filename)
          inVec.vread(inFile)

          DataStore.append(inVec)
          inFile.close()
          if(debug):
               print("read file no # %d"%(i))

     print("Cells read in DataStore: ", len(DataStore))

     return DataStore
##------end of function----------------------##




def save_results2(DataStore):
     """This function is used to save the firing time of the cells 
     in files named after the cells M.H Feb 2018"""

     print("\n||-----SAVING FORCE files by Cell:Type_Name-----||\n")
     for i in range(len(DataStore)):
          tempV 		= np.array(DataStore[i])
          filename 	= "Force_File%d_%d.txt" %(tempV[0],tempV[1])
          np.savetxt(filename,tempV,fmt='%f')
          print("%2d: %s has been saved" % (i,filename))
          
     print("\n||-----Data Saved-----||\n")
##End of Function
############################################################################

ForceStore = []


def ForceGentest(fileIndex):
     """Run on single cell or many, but not suitable for NSG"""

     DataStore = BinarySpikesReader(fileIndex)
    
     for i in range(0,len(DataStore)):
         print(i)
         job_id, Cell_type, cell_num, ReturnVec = ForceGen(
             DataStore[i], DataStore[i][0], DataStore[i][1])
         print("Finished job %2d , CellType %2d , CellNumber %2d" %(job_id, Cell_type, cell_num))
         sys.stdout.flush()
         ForceStore.append(ReturnVec)
     
     save_results2(ForceStore)
##------end of function----------------------##


# //Necessary to make this whole process work.
# code between pc.runworker() and PC.done()  is executed only by the master
pc.runworker()

    # For simulation timing
start_time = time.time()

def ForceGenMain(fileIndex = -1):
     """if given a list of files index, it will only open them, if no input is given it will open all files"""
     ## reading files
     print("Number of cells to read is [S, FR, F] = ",Tot_numcells)
     
     DataStore = BinarySpikesReader(fileIndex)

     global ForceStore
     ForceStore = []  # reset ForceStore to empty
     
     # Calls function within parallel context.  Arg goes after f(x) name
     for i in range(0, len(DataStore)):
          # print("send job", i)
          pc.submit(ForceGen, DataStore[i], DataStore[i][0], DataStore[i][1])
          
     while pc.working():
          job_id, Cell_type, cell_num, returnvec = pc.pyret()
          print("Finished job %3d , CellType %3d , CellNumber %3d" %(job_id ,Cell_type,cell_num))
          sys.stdout.flush()
          ForceStore.append(returnvec)

     ## no need to return as ForceStore is global variable
     
##------end of function----------------------##

## calling main routine, (work on NSG or single core PC)
ForceGenMain()

pc.done()
# all simulations have been completed
# and the workers have been released
# but the boss still has things to do


save_results2(ForceStore)
print("Number of cell that did not fire is %2d out of %2d" %(CellnotFired,len(ForceStore)))
print_time(time.time()-start_time)
