# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
from odbAccess import *
from abaqusConstants import *
import numpy as np
import os
from textRepr import prettyPrint

# --------------------------------------------------------------- STRESS OPTIONS
def htpp_data_stress_option(odb,option):
	step = '\n'.join(s for s in option if 'step=' in s).replace(",", "").replace("\n", "")[6:-1].strip().split(';')
	step = [int(s) for s in step]
	frame = '\n'.join(s for s in option if 'frame=' in s).replace(",", "").replace("\n", "")[7:-1].strip().split(';')
	frame = [int(i) for i in frame]
	elem_set = '\n'.join(s for s in option if 'elemSet=' in s).replace(",", "").replace("\n", "")[9:-1].split(';')
	yls='\n'.join(s for s in option if 'num=' in s).replace(",", "").replace("\n", "")[5:-1].split(';')
	return step,frame,yls,elem_set

# ------------------------------------------------------------------- STRESS AUX
def htpp_data_stress_aux(vals):
	S1 = np.empty(len(vals))
	S2 = np.empty(len(vals))
	SMises = np.empty(len(vals))
	for i in range(len(vals)):
		dim = len(vals[i].data)
		if dim <= 4:
			S1[i] = vals[i].maxInPlanePrincipal
			S2[i] = vals[i].minInPlanePrincipal
			SMises[i] = vals[i].mises
		elif dim > 4:
			stress = np.zeros((3,3))
			stress[0,0] = vals[i].data[0]
			stress[1,1] = vals[i].data[1]
			stress[2,2] = vals[i].data[2]
			stress[0,1] = vals[i].data[3]
			stress[1,0] = stress[0,1]
			stress[0,2] = vals[i].data[4]
			stress[2,0] = stress[0,2]
			stress[1,2] = vals[i].data[5]
			stress[2,1] = stress[1,2]
			w, v = np.linalg.eig(stress)
			S1[i] = max(w[0:2])
			S2[i] = min(w[0:2])

	return S1,S2,SMises

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRESS
# ------------------------------------------------------------------------------
def htpp_data_stress(odb,ioi,option,section,odb_name, info, param):
	step,frame,yls,elem_set = htpp_data_stress_option(odb,option)
	step_names = odb.steps.keys()

	pos = CENTROID

	dir = 'htpp_output/' + odb_name + '/'

	if info == 1 :

		for l in step:
			sT = odb.steps[step_names[l-1]]
			for set in elem_set:
				if set != '':
					soi = ioi.elementSets[set]
					for f in frame:
						fRa = sT.frames[f]
						odb_data = fRa.fieldOutputs['S']
						if len(section) == 0:
							vals = odb_data.getSubset(region=soi, position=pos).values
						else:
							vals = odb_data.getSubset(region=soi,position=pos, sectionPoint=section[-1]).values
						S1,S2,SMises = htpp_data_stress_aux(vals)
						strings = ['stress',set.lower(),'s'+str(l),'f'+str(f)]
						fname = '_'.join(strings)
						with open(dir + fname + '.dat','w') as file:
							file.write('S1 , S2\n')
							for i in range(len(vals)):
								file.write('%f , %f \n' %(S1[i],S2[i]))
		
						strings = ['smises',set.lower(),'s'+str(l),'f'+str(f)]
						fname = '_'.join(strings)
						with open(dir + fname + '.dat','w') as file:
							file.write('SMises\n')
							for i in range(len(vals)):
								file.write('%f \n' %(SMises[i]))


		soi = ioi.elementSets['SET-1']
		pos = INTEGRATION_POINT

		secs = 1
		if len(section) > 0:
			secs = len(section)

		sT = odb.steps['Step-1']
		for f in frame:
			for ls in range(f,f+1):
				fRa = sT.frames[ls]
				odb_data = fRa.fieldOutputs['S']
				strains = odb_data.getSubset(region=soi,position=pos)
				strings = ['stressxxyyxy',set.lower(),'s'+str(l),'f'+str(ls)]
				fname = '_'.join(strings)
				with open(dir + fname + '.dat','w') as file:
					file.write('Sxx, Syy, Sxy\n')
					for v in strains.values:
						file.write('%f, %f, %f\n' %(v.data[0],v.data[1],v.data[3]))
	  
		# soi = ioi.elementSets['SET-1'].elements[0][24768]
		# strings = ['path_peeq','24768']
		# fname = '_'.join(strings)
		# with open(dir + fname + '.dat','w') as file:
		# 	file.write('PEEQ\n')
		# 	for l in step:
		# 		sT = odb.steps[step_names[l-1]]
		# 		for f in range(0,93):
		# 			fRa = sT.frames[f]
		# 			odb_data = fRa.fieldOutputs['SDV1']
		# 			if len(section) == 0:
		# 				vals = odb_data.getSubset(region=soi, position=pos).values
		# 			else:
		# 				vals = odb_data.getSubset(region=soi,position=pos, sectionPoint=section[-1]).values
		# 			#S1,S2,SMises = htpp_data_stress_aux(vals)
		# 			file.write('%f\n' %vals[0].data)

	elif info != 1:

			soi = ioi.elementSets['SET-1']
			pos = INTEGRATION_POINT

			secs = 1
			if len(section) > 0:
				secs = len(section)

			sT = odb.steps['Step-1']
			for f in frame:
				fRa = sT.frames[f]
				odb_data = fRa.fieldOutputs['S']
				strains = odb_data.getSubset(region=soi,position=pos)
				fname = 'ident_stressxxyyxy'+'param_'+str(param)+'.dat'
				with open(dir + fname + '.dat','w') as file:
					file.write('Sxx, Syy, Sxy\n')
					for v in strains.values:
						file.write('%f, %f, %f\n' %(v.data[0],v.data[1],v.data[3]))

# ------------------------------------------------------------------------------
