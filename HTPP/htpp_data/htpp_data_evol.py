# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
from odbAccess import *
from abaqusConstants import *
import numpy as np
import sys
from textRepr import prettyPrint

# ----------------------------------------------------------------- PEEQ OPTIONS
def htpp_data_evol_option(odb,option):
	step = '\n'.join(s for s in option if 'step=' in s).replace(",", "").replace("\n", "")[6:-1].strip().split(';')
	step = [int(s) for s in step]
	frame = '\n'.join(s for s in option if 'frame=' in s).replace(",", "").replace("\n", "")[7:-1].strip().split(';')
	frame = [int(i) for i in frame]
	elem_set = '\n'.join(s for s in option if 'elemSet=' in s).replace(",", "").replace("\n", "")[9:-1].split(';')

	return step,frame,elem_set

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------- PEEQ
# ------------------------------------------------------------------------------
def htpp_data_evol(odb,ioi,option,section,odb_name):
	step,frame,elem_set = htpp_data_evol_option(odb,option)

	soi = ioi.elementSets['SET-1']
	pos = WHOLE_ELEMENT

	dir = 'htpp_output/' + odb_name + '/'

	secs = 1
	if len(section) > 0:
		secs = len(section)

	sT = odb.steps['Step-1']
	fRa = sT.frames[52]
	odb_data = fRa.fieldOutputs['EVOL']
	evol = odb_data.getSubset(region=soi,position=pos)
	with open(dir + 'EVOL.dat','w') as file:
		for v in evol.values:
			file.write('%f\n' %(v.data))

