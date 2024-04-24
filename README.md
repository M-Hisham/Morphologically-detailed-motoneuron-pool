# Morphologically detailed motoneuron-pool
Morphologically detailed motoneuron pool for the vastus lateralis muscle. the model is composed of three cell types: fast fatigable (FF), slow (S), and fast fatigue resistant (FR). The model is based on the [Hodgkin-Huxley model](https://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model) with some modifications to the active channels and dendritic mechanisms.

The model is based on [Allen & Elbasiouny (2018)](https://iopscience.iop.org/article/10.1088/1741-2552/aa9db5)

<span style="display: block; margin-left: auto; margin-right: auto; width: 80%;">![Model structure](imgs/Model.png "Model structure")</span>  

<span style="font-size: 14px; display: block; margin-left: auto; margin-right: auto; width: 90%;">Model schematic showing three distinct cell types (S, FR, FF). model can be stimulated with somatic currents or synaptic pulses. the generated spikes times are saved in output files, where is can be used to generate the force of each motor unit.</span>  

### Cells classes:
| Cell Type | Cell Type Code |
| --- | --- |
| FF | 30 , 3 , 300|
| FR | 20, 2, 200 |
| S | 10, 1, 100 |




## Requirments:
1. Requires Python3.8 (Anaconda package is recommended)
2. Install the [NEURON simulator](https://neuron.yale.edu/neuron/what_is_neuron)
3. Install the packages: numpy, mpi4py, matplotlib, pickle

## Files/folders description
| File/folder | Description |
| --- | --- |
| Simulation #1-8 | Contain .dat and .txt files for the spike times of each MU. The .dat files are used to generate the force vectors |
| Src code | folder contains the source code for the model |
| Simulations data of panels C & D.xlsx | Excel files contains the summary data of all simulations presented in figures 3-8 panels C & D |
| `forceGen.py` | Contains the force generation functions |

<br>

### Simulations source code files:
| File | Function |
| --- | --- |
| `AutoStart.py` | main file to run the simulation |
| `BBGlobals.pyc` | contains the global variables |
| `BBSysI_v3.pyc` | contains synaptic input functions to create and set the synaptic pulses |
| `cellBuilder.pyc` | contain functions for setting up the motoneurons cell templates |
| other .pkl files | files contain information for the cell templates (classes) |
| .mod | files are for the active channels, and dendritic mechanisms|


## How to run the simulation:
### From the command line:
1. Run the simulation using the following command:
          
     ```shell
     python -i AutoStart.py
     ```

<br>

### To run on the [Neuroscience Gateway](https://www.nsgportal.org/):
**Note:** files are written to be compatible with parallel context.
1. Upload the files to the gateway (mod files must be in the parent folder)
2. Create a new task, select "NEURON on Expanse" as the tool.
3. Set the parameters as follow:
     - Run time : 10 Hours (estimate for Expanse)
     - Main file: AutoStart.py
     - Nodes & cores: 3 and 60 (for a 177 cell on Expanse)

<br>

### Force generations from spikes:
After the execution of the simulation, .txt, and  .dat files containing the spikes time of each motor unit will be generated. To generate the respective force of each motor unit:
1. Open the `forceGen.py`, and double check the number of the cells for each type. Variable (Tot_numcells) 
2. Open the terminal then run the `forceGen.py` script using the following command:
          
      ```shell
      python -i forceGen.py
      ```
3. force file for each motor unit will be generated in the same folder
<br>







## Biophysical properties of the model:
<span style="display: block; margin-left: auto; margin-right: auto; width: 70%;">![Scatter matrix for the model input resistance, Rheobase, AHP depth, and AHP half-decay](imgs/Model-properties-scatter-matrix-Rin.png "Summary Scatter matrix")</span>  

<span style="font-size: 16px; display: block; margin-left: auto; margin-right: auto; width: 80%;">Scatter matrix for the model input resistance, Rheobase, AHP depth, and AHP half-decay.</span>  


<span style="display: block; margin-left: auto; margin-right: auto; width: 70%;">![Summary of the model properties in parallel coordinates](imgs/Model-properties-Rin-PC.png "Summary Scatter matrix")</span> 

<span style="font-size: 16px; display: block; margin-left: auto; margin-right: auto; width: 80%;">Parallel coordinates visual for the model input resistance, Rheobase, AHP depth, and AHP half-decay.</span> 


<!-- ## Published article
* [JNP Article](https://journals.physiology.org/doi/full/10.1152/jn.00543.2020)

[![DOI](https://zenodo.org/badge/293670752.svg)](https://zenodo.org/badge/latestdoi/293670752) -->
