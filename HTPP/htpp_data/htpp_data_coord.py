# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
from re import I
from odbAccess import *
from abaqusConstants import *
import numpy as np
import sys
from textRepr import prettyPrint



# ----------------------------------------------------------------- PEEQ OPTIONS
def htpp_data_coord_option(odb,option):
	step = '\n'.join(s for s in option if 'step=' in s).replace(",", "").replace("\n", "")[6:-1].strip().split(';')
	step = [int(s) for s in step]
	frame = '\n'.join(s for s in option if 'frame=' in s).replace(",", "").replace("\n", "")[7:-1].strip().split(';')
	frame = [int(i) for i in frame]
	elem_set = '\n'.join(s for s in option if 'elemSet=' in s).replace(",", "").replace("\n", "")[9:-1].split(';')

	return step,frame,elem_set

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------- PEEQ
# ------------------------------------------------------------------------------
def htpp_data_coord(odb,ioi,option,section,odb_name):
	step,frame,elem_set = htpp_data_coord_option(odb,option)

	soi = ioi.elementSets['SET-1']
	pos = INTEGRATION_POINT

	dir = 'htpp_output/' + odb_name + '/'

	secs = 1
	if len(section) > 0:
		secs = len(section)
	for f in frame:
		sT = odb.steps['Step-1']
		fRa = sT.frames[f]
		odb_data = fRa.fieldOutputs['COORD']
		coords = odb_data.getSubset(region=soi,position=pos)
		with open(dir + 'COORD.dat','w') as file:
			file.write('COORD1, COORD2, COORD3\n')
			for v in coords.values:
				file.write('%f, %f, %f\n' %(v.data[0],v.data[1],v.data[2]))
    

	nset = odb.rootAssembly.nodeSets['SET-1']
	if len(section) > 0:
		secs = len(section)
	for f in frame:
		sT = odb.steps['Step-1']
		fRa = sT.frames[f]
		odb_data = fRa.fieldOutputs['U']
		disp = odb_data.getSubset(region=nset,position=NODAL)
		with open(dir + 'DISP.dat','w') as file:
			file.write('U1, U2, U3\n')
			for v in disp.values:
				file.write('%f, %f, %.15f\n' %(v.data[0],v.data[1],v.data[2]))

	secs = 1
	if len(section) > 0:
		secs = len(section)
	for f in frame:
		sT = odb.steps['Step-1']
		fRa = sT.frames[0]
		odb_data = fRa.fieldOutputs['COORD']
		coords = odb_data.getSubset(region=soi,position=pos)
		with open(dir + 'COORD_undeformed.dat','w') as file:
			file.write('COORD1, COORD2, COORD3\n')
			for v in coords.values:
				file.write('%f, %f, %f\n' %(v.data[0],v.data[1],v.data[2]))

	# --------------------------------------------------------------- WRITE DATA
