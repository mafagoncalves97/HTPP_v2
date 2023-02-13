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
def htpp_compute_buckling_name(set, step, frame, num, table):
	name = ['BUCKLING [']
	set = set.replace('_', '\_')
	name.append('set: ' + set.upper() + ' /')
	name.append('step: ' + str(step) + ' /')
	name.append('frame: ' + str(frame))

	name.append(']')
	name = ' '.join(name)

	table.append(name)

	return table

# --------------------------------------------------------------- STRAIN OPTIONS


def htpp_compute_buckling_option(option):
	step = '\n'.join(s for s in option if 'step=' in s).replace(
		",", "").replace("\n", "")[6:-1].strip().split(';')
	step = [int(s) for s in step]
	frame = '\n'.join(s for s in option if 'frame=' in s).replace(
		",", "").replace("\n", "")[7:-1].strip().split(';')
	frame = [int(i) for i in frame]
	elem_set = '\n'.join(s for s in option if 'elemSet=' in s).replace(
		",", "").replace("\n", "")[9:-1].split(';')


	return step, frame, elem_set

# ------------------------------------------------------------------- STRAIN AUX


def htpp_compute_buckling_aux(disp, odb_name):
	vals=np.zeros(2)
	U3 = disp[:,2]
	meanU3 = statistics.mean(U3)
	maxU3 = np.max(U3)
	vals[0] = meanU3
	vals[1] = maxU3
	return vals

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRAIN
# ------------------------------------------------------------------------------


def htpp_compute_buckling(option,odb_name):
	# plt.rc('text', usetex=True)
	plt.rc('font', family='serif', size=14)

	dir = 'htpp_output/' + odb_name + '/'

	step, frame, elem_set = htpp_compute_buckling_option(option)

	for set in elem_set:
		for l in step:
			for f in frame:
				strings = 'DISP'
				fname = dir + strings + '.dat'
				disp = np.loadtxt(fname, skiprows=1, delimiter=',')
	
				vals = htpp_compute_buckling_aux(disp,odb_name)

				strings = ['buckling',set.lower(),'s'+str(l),'f'+str(f)]
				fname = dir+'_'.join(strings)+'.dat'
				with open(fname,'w') as f:
					for i in range(0,len(vals)):
						f.write('%.15f\n' %vals[i])

# ------------------------------------------------------------------------------
