# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from math import floor, ceil

# ----------------------------------------------------------- ROTATION PLOT NAME
def htpp_plot_schmitt_name(set,step,table):
	name = ['SCHMITT [']
	set = set.replace('_','\_')
	name.append('set: '  + set.upper() + ' /')
	name.append('step: ' + str(step) + ' /')


	name.append(']')
	name = ' '.join(name)

	table.append(name)

	return table
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------ YIELD LOCUS
# ------------------------------------------------------------------------------
def htpp_plot_schmitt_option(option):
    
	step = '\n'.join(s for s in option if 'step=' in s).replace(",", "").replace("\n", "")[6:-1].strip().split(';')
	step = [int(s) for s in step]
	
	elem_set = '\n'.join(s for s in option if 'elemSet=' in s).replace(",", "").replace("\n", "")[9:-1].split(';')
	

	return step,elem_set

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------ YIELD LOCUS
# ------------------------------------------------------------------------------
def htpp_plot_schmitt_aux(set,l,dir):
    
	plt.rc('font', family='serif',size=10)
	plt.rc('axes', labelsize=12)
 
 
	sc1=np.array([])
	sc2=np.array([])
	sc3=np.array([])
	fig = plt.figure()
	for f in range(0,66):
		strings = ['schmitt',set.lower(),'s'+str(l),'f'+str(f)]
		fname = dir + '_'.join(strings) + '.dat'
		schmitt = np.loadtxt(fname)
		sc1=np.append(sc1,schmitt[4214-1])
		sc2=np.append(sc2,schmitt[6899-1])
		sc3=np.append(sc3,schmitt[5242-1])
		
	plt.plot(np.arange(0,f+1),sc1,'-.b', label='Element 1')
	plt.plot(np.arange(0,f+1),sc2,'-k', label='Element 2')
	plt.plot(np.arange(0,f+1),sc3,'--c', label='Element 3')

	plt.ylabel('Schmitt factor')
	plt.xlabel('Iterations')

	plt.ylim(-1,1)
	plt.xlim(0,93)
	plt.legend(loc='lower center',ncols=3)
	plt.show()
	return fig

# ------------------------------------------------------------------------------
def htpp_plot_schmitt(odb_name,option,pdf,size,table):

	dir = 'htpp_output/' + odb_name + '/'
	step,elem_set= htpp_plot_schmitt_option(option)

	for set in elem_set:
		for l in step:

			fig = htpp_plot_schmitt_aux(set,l,dir)

		

			table = htpp_plot_schmitt_name(set,l,table)

	return pdf,table
# ------------------------------------------------------------------------------
