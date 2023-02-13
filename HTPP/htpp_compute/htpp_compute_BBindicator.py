# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpt
from matplotlib.legend_handler import HandlerBase, HandlerPatch
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from math import floor, ceil, pi, tan, sin, cos
from scipy.interpolate import BSpline, splrep
import statistics
import os


# ------------------------------------------------------------. STRAIN PLOT NAME
def htpp_compute_BBindicator_name(set, step, frame, num, table):
	name = ['BB_INDICATOR [']
	set = set.replace('_', '\_')
	name.append('set: ' + set.upper() + ' /')
	name.append('step: ' + str(step) + ' /')
	name.append('frame: ' + str(frame))

	name.append(']')
	name = ' '.join(name)

	table.append(name)

	return table

# --------------------------------------------------------------- STRAIN OPTIONS


def htpp_compute_BBindicator_option(option):
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


def htpp_compute_BBindicator_aux(strain, peeq, evol, Smises, odb_name):

	E2E1_ratio = strain[:,1] / strain[:,0]  # MINOR AND MAJOR STRAIN RATIO

	# limiting the range of admissable strain states between [-3,1]
	E2E1_ratio[E2E1_ratio < -3] = -3
	E2E1_ratio[E2E1_ratio > 1] = -3

	# ---------------------------------------  EQUIVALENT PLASTIC STRAIN ANALYSIS  ---------------------------------  #

	if len(strain) > 15 * 4:
		PEEQ_test = np.mean(np.sort(peeq)[len(peeq) - 16:-1])  # equiv plastic strain value of the test, at the end of the test. it is the mean of the
		# 15 higher values of PEEQ if there are alot of elements
	else:
		PEEQ_test = np.amax(peeq)

	PEEQ_tens = np.nan_to_num(np.mean(peeq[np.intersect1d(np.where(E2E1_ratio <= -0.47)[0],
														  np.where(E2E1_ratio >= -0.53)[0])]))  # PEEQ mean for the tensile strain state, at the end of the test
	PEEQ_shear = np.nan_to_num(
		np.mean(peeq[np.intersect1d(np.where(E2E1_ratio <= -0.97)[0], np.where(E2E1_ratio >= -1.03)[0])]))  # PEEQ mean for the shear strain state, at the end of the test
	PEEQ_biaxial = np.nan_to_num(np.mean(peeq[np.where(E2E1_ratio >= 0.97)[0]]))  # PEEQ mean for the biaxial strain state, at the end of the test
	PEEQ_comp = np.nan_to_num(
		np.mean(peeq[np.intersect1d(np.where(E2E1_ratio <= -1.97)[0], np.where(E2E1_ratio >= -2.03)[0])]))  # PEEQ mean for the compressive strain state,  at the end of the test
	PEEQ_plane = np.nan_to_num(
		np.mean(peeq[np.intersect1d(np.where(E2E1_ratio <= 0.03)[0], np.where(E2E1_ratio >= -0.03)[0])]))  # PEEQ mean for the plane strain state,  at the end of the test


	#  -------------------------------------------------------------------------------------------------------------  #
	# --------------------------------------  HETEROGENEITY INDICATOR CALCULATION  ---------------------------------  #
	#  -------------------------------------------------------------------------------------------------------------  #


	SSe = (Smises - (np.mean(Smises))) / (np.mean(Smises))  # identify stress concentrations
	b = 3
	Ze = 1 / (1 + (b * SSe) ** 2)  # stress concentrations penalization
	beta = 20
	s = np.empty([np.size(peeq)])  # compression shear or tension definition of each element
	s[np.where(E2E1_ratio >= -0.75)[0]] = 3  # tension
	s[np.where(E2E1_ratio <= -1.5)[0]] = 1  # compression
	s[np.intersect1d(np.where(E2E1_ratio > -1.5)[0], np.where(E2E1_ratio < -0.75)[0])] = 2  # shear

	delta1 = (1 - (np.tanh(beta * (strain[:,0] + (0.75 * strain[:,1]))))) / 2
	delta2 = ((1 + (np.tanh(beta * (strain[:,0] + (0.75 * strain[:,1]))))) * (1 - (np.tanh(beta * (strain[:,0] + (1.5 * strain[:,1])))))) / 4
	delta3 = (1 + (np.tanh(beta * (strain[:,0] + (1.5 * strain[:,1]))))) / 2

	id1 = (3 / np.sum(evol[np.where(s==1)[0]])) * (np.sum(delta1[np.where(s==1)[0]] * Ze[np.where(s==1)[0]] * evol[np.where(s==1)[0]]))
	id2 = (3 / np.sum(evol[np.where(s==2)[0]])) * (np.sum(delta2[np.where(s==2)[0]] * Ze[np.where(s==2)[0]] * evol[np.where(s==2)[0]]))
	id3 = (3 / np.sum(evol[np.where(s==3)[0]])) * (np.sum(delta3[np.where(s==3)[0]] * Ze[np.where(s==3)[0]] * evol[np.where(s==3)[0]]))

	id1 = sum(id1*evol)/len(evol)
	id2 = sum(id2*evol)/len(evol)
	id3 = sum(id3*evol)/len(evol)
	# ------------------------------------------------  INDICATOR VALUE ---------------------------------------------  #

	indicator = id1 * id2 * id3  # indicator based on stress states
	vals=np.empty(4,dtype=float)
	vals[0] = indicator
	vals[1] = id1
	vals[2] = id2; vals[3] = id3; 
	return vals

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRAIN
# ------------------------------------------------------------------------------


def htpp_compute_BBindicator(option,odb_name):
	# plt.rc('text', usetex=True)
	plt.rc('font', family='serif', size=14)

	dir = 'htpp_output/' + odb_name + '/'

	step, frame, elem_set, lims, txt, num = htpp_compute_BBindicator_option(option)

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
	 
					strings = ['evol', set.lower(), 's'+str(l), 'f'+str(f)]
					fname = dir + 'EVOL' + '.dat'
					evol = np.loadtxt(fname)
	 
					strings = ['smises', set.lower(), 's'+str(l), 'f'+str(f)]
					fname = dir + '_'.join(strings) + '.dat'
					smises = np.loadtxt(fname,skiprows=1)
	 
					vals = htpp_compute_BBindicator_aux(strain,peeq,evol,smises,odb_name)
  
					strings = ['BBindicator',set.lower(),'s'+str(l),'f'+str(f)]
					fname = dir+'_'.join(strings)+'.dat'
					with open(fname,'w') as f:
						for i in range(len(vals)):
							f.write('%f\n' %vals[i])

# ------------------------------------------------------------------------------
