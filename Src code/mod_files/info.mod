COMMENT
This mechanism is intended to store some information about each compartment

DendDist [um]: is the dendritic path distance of the dendritic compartment from the soma (Note: not the radial distance)
               This distance is computed between the center of the compartment and the soma
		   This measurement should be equal to the value stored in DendPathDist(0.5)
               
SpaceConstant [um]: This is the space constant (lambda)

ElectroL - no units: is the electrotonic length of the section (ratio of the compartment length to the space constant)

ElectroDist - no units: is the electrotonic distance of the section from the soma

Terminal - Boolean (0 or 1): Zero means that the compartment is NOT a terminal dendrites
				     ONE means that the compartment is a terminal dendrites

x, y, z:    Allows local storage of xyz coordinates interpolated from the pt3d data.  These coordinates are used by hoc code to interpolate pt3d data
		along the section
ENDCOMMENT


NEURON {
	SUFFIX info
	RANGE DendDist, SpaceConstant, EculadianDist, ElectroDist, Terminal
	RANGE x, y, z
}

PARAMETER {
	DendDist = 0		
	SpaceConstant = 0	

	EculadianDist = 0
	ElectroDist = 0

	Terminal = 0

	x = 0 (1) : spatial coords
	y = 0 (1)
	z = 0 (1)
}