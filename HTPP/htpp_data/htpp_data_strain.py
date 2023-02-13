# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
from odbAccess import *
from abaqusConstants import *
import numpy as np
import os
from textRepr import prettyPrint

# --------------------------------------------------------------- STRAIN OPTIONS
def htpp_data_strain_option(odb,option):
	step = '\n'.join(s for s in option if 'step=' in s).replace(",", "").replace("\n", "")[6:-1].strip().split(';')
	step = [int(s) for s in step]
	frame = '\n'.join(s for s in option if 'frame=' in s).replace(",", "").replace("\n", "")[7:-1].strip().split(';')
	frame = [int(i) for i in frame]
	elem_set = '\n'.join(s for s in option if 'elemSet=' in s).replace(",", "").replace("\n", "")[9:-1].split(';')

	return step,frame,elem_set

# ------------------------------------------------------------------- STRAIN AUX
def htpp_data_strain_aux(vals):
	E1 = np.empty(len(vals))
	E2 = np.empty(len(vals))
	for i in range(len(vals)):
		dim = len(vals[i].data)
		if dim <= 4:
			E1[i] = vals[i].maxInPlanePrincipal
			E2[i] = vals[i].minInPlanePrincipal
		elif dim > 4:
			strain = np.zeros((3,3))
			strain[0,0] = vals[i].data[0]
			strain[1,1] = vals[i].data[1]
			strain[2,2] = vals[i].data[2]
			strain[0,1] = vals[i].data[3]/2.0
			strain[1,0] = strain[0,1]
			strain[0,2] = vals[i].data[4]/2.0
			strain[2,0] = strain[0,2]
			strain[1,2] = vals[i].data[5]/2.0
			strain[2,1] = strain[1,2]
			w, v = np.linalg.eig(strain)
			E1[i] = max(w[0:2])
			E2[i] = min(w[0:2])

	return E1,E2

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRAIN
# ------------------------------------------------------------------------------
def htpp_data_strain(odb,ioi,option,section,odb_name, info, param):

	if info == 1: 
		step,frame,elem_set = htpp_data_strain_option(odb,option)

		step_names = odb.steps.keys()

		# --------------------------------------------------------------- WRITE DATA
		pos = CENTROID
		dir = 'htpp_output/' + odb_name + '/'
		for l in step:
			sT = odb.steps[step_names[l-1]]
			for set in elem_set:
				if set != '':
					soi = ioi.elementSets[set]
					for f in frame:
						fRa = sT.frames[f]
						odb_data = fRa.fieldOutputs['LE']
						if len(section) == 0:
							vals = odb_data.getSubset(region=soi, position=pos).values
						else:
							vals = odb_data.getSubset(region=soi,position=pos, sectionPoint=section[-1]).values
						E1,E2 = htpp_data_strain_aux(vals)
						strings = ['strain',set.lower(),'s'+str(l),'f'+str(f)]
						fname = '_'.join(strings)
						with open(dir + fname + '.dat','w') as file:
							file.write('E1 , E2\n')
							for i in range(len(vals)):
								file.write('%f , %f \n' %(E1[i],E2[i]))


		soi = ioi.elementSets['SET-1']
		pos = INTEGRATION_POINT

		dir = 'htpp_output/' + odb_name + '/'

		secs = 1
		if len(section) > 0:
			secs = len(section)

		sT = odb.steps['Step-1']
		for f in frame:
			fRa = sT.frames[f]
			odb_data = fRa.fieldOutputs['LE']
			strains = odb_data.getSubset(region=soi,position=pos)
			strings = ['strainxxyyxy',set.lower(),'s'+str(l),'f'+str(f)]
			fname = '_'.join(strings)
			with open(dir + fname + '.dat','w') as file:
				file.write('LExx, LEy, LExy\n')
				for v in strains.values:
					file.write('%f, %f, %f\n' %(v.data[0],v.data[1],v.data[3]))

	elif info != 1:

		soi = ioi.elementSets['SET-1']
		pos = INTEGRATION_POINT

		dir = 'htpp_output/' + odb_name + '/'

		secs = 1
		if len(section) > 0:
			secs = len(section)

		sT = odb.steps['Step-1']
		for f in frame:
			fRa = sT.frames[f]
			odb_data = fRa.fieldOutputs['LE']
			strains = odb_data.getSubset(region=soi,position=pos)
			fname = 'ident_strainxxyyxy'+'param_'+str(param)+'.dat'
			with open(dir + fname + '.dat','w') as file:
				file.write('LExx, LEy, LExy\n')
				for v in strains.values:
					file.write('%f, %f, %f\n' %(v.data[0],v.data[1],v.data[3]))

# ------------------------------------------------------------------------------
