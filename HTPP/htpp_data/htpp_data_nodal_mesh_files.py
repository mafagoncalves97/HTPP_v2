# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
from re import I
from odbAccess import *
from abaqusConstants import *
import numpy as np
import sys
from textRepr import prettyPrint
import os


# ----------------------------------------------------------------- PEEQ OPTIONS
def htpp_data_nodal_mesh_files_option(odb,option):
	step = '\n'.join(s for s in option if 'step=' in s).replace(",", "").replace("\n", "")[6:-1].strip().split(';')
	step = [int(s) for s in step]
	frame = '\n'.join(s for s in option if 'frame=' in s).replace(",", "").replace("\n", "")[7:-1].strip().split(';')
	frame = [int(i) for i in frame]
	elem_set = '\n'.join(s for s in option if 'elemSet=' in s).replace(",", "").replace("\n", "")[9:-1].split(';')

	return step,frame,elem_set

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------- PEEQ
# ------------------------------------------------------------------------------
def htpp_data_nodal_mesh_files(odb,ioi,option,section,odb_name):
	step,frame,elem_set = htpp_data_nodal_mesh_files_option(odb,option)

	nset = odb.rootAssembly.nodeSets['SET-1']

	dir = 'htpp_output/' + odb_name + '/'
	path = os.getcwd()[:-4] + '/'
	secs = 1
	if len(section) > 0:
		secs = len(section)
	for f in frame:
		
		for i in range(0,f+1):
			
			try:
				os.remove(dir+'nodal_files_'+odb_name+'_%d.csv' %i)
			except:
				pass
			

			errors = odb.diagnosticData.numberOfAnalysisErrors
			if errors != 0:
				return

			fmt1 = '%d;%.6f;%.6f;%.6f'

			f = open(dir+'nodal_files_'+odb_name+'_%d.csv' %i,'a')

			last_frame = odb.steps['Step-1'].frames[i]  # get last frame of all

			Coords = last_frame.fieldOutputs['COORD'].getSubset(position=NODAL, region=nset)

			f.write('*Part, name=Part-1\n*Nodes\n')  # save variables to file
			l=1
			for v in Coords.values:
				coordX = v.data[0]
				coordY = v.data[1]
				coordZ = v.data[2]
				np.savetxt(f, np.column_stack((l,coordX, coordY, coordZ)), fmt=fmt1, delimiter=',')  # save variables to file
				l=l+1

			f.write('*End Part')  # save variables to file
			f.close()


		for i in range(0,1):
			
			try:
				os.remove(dir+odb_name+'_mesh.mesh')
			except:
				pass
			
			errors = odb.diagnosticData.numberOfAnalysisErrors
			if errors != 0:
				return

			fmt1 = '%d;%.6f;%.6f;%.6f'

			f = open(dir+odb_name+'_mesh.mesh','a')

			last_frame = odb.steps['Step-1'].frames[i]  # get last frame of all

			Coords = last_frame.fieldOutputs['COORD'].getSubset(position=NODAL, region=nset)

			f.write('*Part, name=Part-1\n*Nodes\n')  # save variables to file
			l=1
			for v in Coords.values:
				coordX = v.data[0]
				coordY = v.data[1]
				coordZ = v.data[2]
				np.savetxt(f, np.column_stack((l,coordX, coordY, coordZ)), fmt=fmt1, delimiter=',')  # save variables to file
				l=l+1

			f.write('*End Part')  # save variables to file
			f.close()

	# --------------------------------------------------------------- WRITE DATA
