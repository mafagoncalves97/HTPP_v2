# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import matplotlib.patches as mpt
#from matplotlib.legend_handler import HandlerBase, HandlerPatch
#from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from math import floor, ceil, pi, tan, sin, cos
#from scipy.interpolate import BSpline, splrep
import statistics
import os


# ------------------------------------------------------------. STRAIN PLOT NAME
def htpp_compute_indicator_name(set, step, frame, num, table):
	name = ['INDICATOR [']
	set = set.replace('_', '\_')
	name.append('set: ' + set.upper() + ' /')
	name.append('step: ' + str(step) + ' /')
	name.append('frame: ' + str(frame))

	name.append(']')
	name = ' '.join(name)

	table.append(name)

	return table

# --------------------------------------------------------------- STRAIN OPTIONS


def htpp_compute_indicator_option(option):
	step = '\n'.join(s for s in option if 'step=' in s).replace(
		",", "").replace("\n", "")[6:-1].strip().split(';')
	step = [int(s) for s in step]
	frame = '\n'.join(s for s in option if 'frame=' in s).replace(
		",", "").replace("\n", "")[7:-1].strip().split(';')
	frame = [int(i) for i in frame]
	elem_set = '\n'.join(s for s in option if 'elemSet=' in s).replace(
		",", "").replace("\n", "")[9:-1].split(';')
	lims = '\n'.join(s for s in option if 'lims=' in s).replace(
		"\n", "")[6:-1].split(';')
	tmp = []
	if lims[0] != '':
		for lim in lims:
			l = lim[1:-1].split(',')
			tmp.append([float(l[0]), float(l[1])])
	lims = tmp
	txt = '\n'.join(s for s in option if 'text=' in s).replace(
		"\n", "")[6:-1].split(',')
	if txt[0] == '':
		txt = []
	num = '\n'.join(s for s in option if 'num=' in s).replace(
		",", "").replace("\n", "")[5:-1].split(';')
	num = [float(n) for n in num]

	return step, frame, elem_set, lims, txt, num

# ------------------------------------------------------------------- STRAIN AUX


def htpp_compute_indicator_aux(strain, peeq, evol, ind_terms, odb_name):

	E2E1_ratio = strain[:,1] / strain[:,0]  # MINOR AND MAJOR STRAIN RATIO

	# limiting the range of admissable strain states between [-3,1]
	E2E1_ratio[E2E1_ratio < -15] = -15
	E2E1_ratio[E2E1_ratio > 1] = 1
	E2E1_ratio=np.delete(E2E1_ratio,np.where(np.isnan(E2E1_ratio)==True))

	E2E1_ratio_Max = np.amax(E2E1_ratio)
	E2E1_ratio_Min = np.amin(E2E1_ratio)
	E2E1_ratio_R = E2E1_ratio_Max - E2E1_ratio_Min  # STRAIN STATE RANGE at the end of the test
	E2E1_dispersion = statistics.stdev(E2E1_ratio.tolist())  # STANDARD DEVIATION OF STRAIN STATE (dispersion of E2E1 ratio)
	#print(E2E1_dispersion, E2E1_ratio_Min,E2E1_ratio_Max)
	# ---------------------------------------  EQUIVALENT PLASTIC STRAIN ANALYSIS  ---------------------------------  #

	PEEQ_dispersion = statistics.stdev(peeq.tolist())  # STANDARD DEVIATION OF PEEQ (dispersion of PEEQ) 

	AvPEEQ = np.sum((peeq * evol)) / (np.sum(evol))  # AVERAGE DEFORMATION at the end of the test
	#print(PEEQ_dispersion, AvPEEQ)
	if len(strain) > 15*4:
		PEEQ_test = np.mean(np.sort(peeq)[len(peeq) - 16:-1])  # equiv plastic strain value of the test, at the end of the test. it is the mean of the
		# 15 higher values of PEEQ if there are alot of elements
	else:
		PEEQ_test = np.amax(peeq)

	PEEQ_tens = np.nan_to_num(np.mean(peeq[np.intersect1d(np.where(E2E1_ratio <= -0.26)[0], np.where(E2E1_ratio >= -0.75)[0])]))  # PEEQ mean for the tensile strain state, at the end of the test
	PEEQ_shear = np.nan_to_num(np.mean(peeq[np.intersect1d(np.where(E2E1_ratio <= -0.76)[0], np.where(E2E1_ratio >= -1.5)[0])]))  # PEEQ mean for the shear strain state, at the end of the test
	PEEQ_biaxial = np.nan_to_num(np.mean(peeq[np.where(E2E1_ratio >= 0.75)[0]]))  # PEEQ mean for the biaxial strain state, at the end of the test
	PEEQ_comp = np.nan_to_num(np.mean(peeq[np.intersect1d(np.where(E2E1_ratio <= -1.5)[0], np.where(E2E1_ratio >= -2.5)[0])]))  # PEEQ mean for the compressive strain state,  at the end of the test
	PEEQ_plane = np.nan_to_num(np.mean(peeq[np.intersect1d(np.where(E2E1_ratio <= 0.25)[0], np.where(E2E1_ratio >= -0.25)[0])]))  # PEEQ mean for the plane strain state,  at the end of the test
	PEEQ_max = (PEEQ_test + PEEQ_tens + PEEQ_shear + PEEQ_biaxial + PEEQ_comp + PEEQ_plane) / 6  # STRAIN LEVEL
 
	'''PEEQ_tens = np.nan_to_num(np.sum(peeq[np.intersect1d(np.where(E2E1_ratio <= -0.26)[0], np.where(E2E1_ratio >= -0.75)[0])]*evol[np.intersect1d(np.where(E2E1_ratio <= -0.26)[0], np.where(E2E1_ratio >= -0.75)[0])]))/(len(evol))  # PEEQ mean for the tensile strain state, at the end of the test
	PEEQ_shear = np.nan_to_num(np.sum(peeq[np.intersect1d(np.where(E2E1_ratio <= -0.76)[0], np.where(E2E1_ratio >= -1.5)[0])]*evol[np.intersect1d(np.where(E2E1_ratio <= -0.76)[0], np.where(E2E1_ratio >= -1.5)[0])]))/(len(evol))  # PEEQ mean for the shear strain state, at the end of the test
	PEEQ_biaxial = np.nan_to_num(np.sum(peeq[np.where(E2E1_ratio >= 0.75)[0]]*evol[np.where(E2E1_ratio >= 0.75)[0]]))/(len(evol))  # PEEQ mean for the biaxial strain state, at the end of the test
	PEEQ_comp = np.nan_to_num(np.sum(peeq[np.intersect1d(np.where(E2E1_ratio <= -1.5)[0], np.where(E2E1_ratio >= -2.5)[0])]*evol[np.intersect1d(np.where(E2E1_ratio <= -1.5)[0], np.where(E2E1_ratio >= -2.5)[0])]))/(len(evol))  # PEEQ mean for the compressive strain state,  at the end of the test
	PEEQ_plane = np.nan_to_num(np.sum(peeq[np.intersect1d(np.where(E2E1_ratio <= 0.25)[0], np.where(E2E1_ratio >= -0.25)[0])]*evol[np.intersect1d(np.where(E2E1_ratio <= 0.25)[0], np.where(E2E1_ratio >= -0.25)[0])]))/(len(evol))  # PEEQ mean for the plane strain state,  at the end of the test
	PEEQ_max = (PEEQ_test + PEEQ_tens + PEEQ_shear + PEEQ_biaxial + PEEQ_comp + PEEQ_plane) / 6  # STRAIN LEVEL'''

	#print(PEEQ_test,PEEQ_tens,PEEQ_shear,PEEQ_plane,PEEQ_biaxial,PEEQ_comp,PEEQ_max)
	
	#  -------------------------------------------------------------------------------------------------------------  #
	# --------------------------------------  HETEROGENEITY INDICATOR CALCULATION  ---------------------------------  #
	#  -------------------------------------------------------------------------------------------------------------  #

	# ------------------------------------  LOAD RELATIVE AND ABSOLUTE TERM'S VALUES  ------------------------------  #
	
	Wr1 = ind_terms[0]
	Wr2 = ind_terms[1]
	Wr3 = ind_terms[2]
	Wr4 = ind_terms[3]
	Wr5 = ind_terms[4]
	Wa1 = ind_terms[5]
	Wa2 = ind_terms[6]
	Wa3 = ind_terms[7]
	Wa4 = ind_terms[8]
	Wa5 = ind_terms[9]

	# ------------------------------  INDICATOR TERMS, WITHOUT normalization and relation --------------------------  #

	A1 = E2E1_dispersion
	A2 = sum(E2E1_ratio_R*evol)/(len(evol))
	A3 = PEEQ_dispersion
	A4 = sum(PEEQ_max*evol)/(len(evol))
	A5 = AvPEEQ

	# ----------------------------  TERMS OF THE INDICATOR, WITH normalization and relation ------------------------  #

	I1 = (Wr1 * A1) / Wa1
	I2 = (Wr2 * A2) / Wa2
	I3 = (Wr3 * A3) / Wa3
	I4 = (Wr4 * A4) / Wa4
	I5 = (Wr5 * A5) / Wa5

	# ------------------------------------------------  INDICATOR VALUE ---------------------------------------------  #

	indicator = I1 + I2 + I3 + I4 + I5
	vals=np.empty(12,dtype=float)
	vals[0] = indicator
	vals[1] = A1; vals[2] = A2; vals[3] = A3; vals[4] = A4; vals[5] = A5; 
	vals[6] = sum(PEEQ_tens*evol)/(len(evol)); vals[7] = sum(PEEQ_shear*evol)/(len(evol)); vals[8] = sum(PEEQ_biaxial*evol)/(len(evol)); vals[9] = sum(PEEQ_comp*evol)/(len(evol)); vals[10] = sum(PEEQ_plane*evol)/(len(evol))

	#print(indicator)
	return vals

# ------------------------------------------------------------------- STRAIN AUX


def htpp_compute_strain_ratio(strain, odb_name):
	
	E2E1_ratio = strain[:,1] / strain[:,0]  # MINOR AND MAJOR STRAIN RATIO
	
	#E2E1_ratio=np.delete(E2E1_ratio,np.where(np.isnan(E2E1_ratio)==True))
	
	return E2E1_ratio
	
	
	
# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRAIN
# ------------------------------------------------------------------------------


def htpp_compute_indicator(option,odb_name):
	# plt.rc('text', usetex=True)
	#plt.rc('font', family='serif', size=14)

	dir = 'htpp_output/' + odb_name + '/'

	step, frame, elem_set, lims, txt, num = htpp_compute_indicator_option(option)

	for set in elem_set:
		for l in step:
			for f in frame:
				for n in num:
					strings = ['strain', set.lower(), 's'+str(l), 'f'+str(f)]
					fname = dir + '_'.join(strings) + '.dat'
					strain = np.loadtxt(fname, skiprows=1, delimiter=',')

					strings = ['peeq', set.lower(), 's'+str(l), 'f'+str(f)]
					fname = dir + '_'.join(strings) + '.dat'
					peeq = np.loadtxt(fname)
	 
					fname = dir + 'EVOL' + '.dat'
					evol = np.loadtxt(fname)
	 
					fname = 'ind_param'
					path = os.getcwd()[:-4] + '/'
					fname = path + fname + '.ind'
					ind_terms = np.loadtxt(fname)
	 
					vals = htpp_compute_indicator_aux(strain,peeq,evol,ind_terms,odb_name)
  
					strings = ['indicator',set.lower(),'s'+str(l),'f'+str(f)]
					fname = dir+'_'.join(strings)+'.dat'
					with open(fname,'w') as fi:
						for i in range(len(vals)):
							fi.write('%f\n' %vals[i])

					E2E1_ratio = htpp_compute_strain_ratio(strain, odb_name)
					strings = ['ratio_E2E1',set.lower(),'s'+str(l),'f'+str(f)]
					fname = dir+'_'.join(strings)+'.dat'
					with open(fname,'w') as fi:
						for i in range(len(E2E1_ratio)):
							fi.write('%f\n' %E2E1_ratio[i])
# ------------------------------------------------------------------------------
