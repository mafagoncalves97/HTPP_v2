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

# ------------------------------------------------------------. STRAIN PLOT NAME
def htpp_plot_BBindicator_name(set,step,frame,num,table):
	name = ['BBINDICATOR [']
	set = set.replace('_','\_')
	name.append('set: '  + set.upper() + ' /')
	name.append('step: ' + str(step) + ' /')
	name.append('frame: '+ str(frame))

	name.append(']')
	name = ' '.join(name)

	table.append(name)

	return table

# --------------------------------------------------------------- STRAIN OPTIONS
def htpp_plot_BBindicator_option(option):
	step = '\n'.join(s for s in option if 'step=' in s).replace(",", "").replace("\n", "")[6:-1].strip().split(';')
	step = [int(s) for s in step]
	frame = '\n'.join(s for s in option if 'frame=' in s).replace(",", "").replace("\n", "")[7:-1].strip().split(';')
	frame = [int(i) for i in frame]
	elem_set = '\n'.join(s for s in option if 'elemSet=' in s).replace(",", "").replace("\n", "")[9:-1].split(';')
	lims = '\n'.join(s for s in option if 'lims=' in s).replace("\n" , "")[6:-1].split(';')
	tmp = []
	if lims[0] != '':
		for lim in lims:
			l = lim[1:-1].split(',')
			tmp.append([float(l[0]),float(l[1])])
	lims = tmp
	txt = '\n'.join(s for s in option if 'text=' in s).replace("\n","")[6:-1].split(',')
	if txt[0] == '':
		txt = []
	num = '\n'.join(s for s in option if 'num=' in s).replace(",", "").replace("\n", "")[5:-1].split(';')
	num = [float(n) for n in num]

	return step,frame,elem_set,lims,txt,num

# ------------------------------------------------------------------- STRAIN AUX
def htpp_plot_BBindicator_aux(ind_notched,ind_kim,ind_jones,ind_conde,ind_topopt,lims,txt,n,odb_name):

	
	fig = plt.figure(figsize=(5,3))
	barWidth=0.8
	#plt.rc('text', usetex=1)
	plt.rc('font', family='serif',size=8)
	plt.rc('axes', labelsize=8)
	c = ['lightgrey',  'darkgrey', 'dimgrey']
	plt.subplots_adjust(wspace= 0.95, hspace= 0.35)

	sub1 = fig.add_subplot(1,4,1)
	sub1.set_title(r'Tension')
	IT1 = [ind_notched[3], ind_jones[3],  ind_topopt[3]]
	sub1.set_xticks([])
	br1 = np.arange(len(IT1))
	plt.bar(br1, IT1, color = c, width = barWidth,edgecolor ='None')
		

	sub2 = fig.add_subplot(1,4,2)
	sub2.set_title(r'Compression')
	IT2 = [ind_notched[1],  ind_jones[1], ind_topopt[1]]
	sub2.set_xticks([])
	br2 = np.arange(len(IT2))
	plt.bar(br2, IT2, color = c, width = barWidth,edgecolor ='None')

	sub3 = fig.add_subplot(1,4,3)
	sub3.set_title(r'Shear')
	IT3 = [ind_notched[2], ind_jones[2],  ind_topopt[2]]
	sub3.set_xticks([])
	br3= np.arange(len(IT3))
	plt.bar(br3, IT3, color = c, width = barWidth,edgecolor ='None')

	sub4 = fig.add_subplot(1,4,4)
	sub4.set_title(r'$I_b$')
	IT4 = [ind_notched[0],  ind_jones[0],  ind_topopt[0]]
	sub4.set_xticks([])
	br4 = np.arange(len(IT4))
	spec =['Notched',  'D','TopOpt']
	plt.bar(br4, IT4, color = c, width = barWidth,edgecolor ='None', label=spec)



	fig.legend(loc='lower center',ncols=5)
	plt.savefig(r'htpp_output/Figures/BB_indicator_specimens.jpg')
	plt.show()
	return fig

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRAIN
# ------------------------------------------------------------------------------
def htpp_plot_BBindicator(odb_name,option,pdf,size,table):
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif',size=14)

	dir = 'htpp_output/' + odb_name + '/'

	step,frame,elem_set,lims,txt,num = htpp_plot_BBindicator_option(option)

	for set in elem_set:
		for l in step:
			for f in frame:
				for n in num:
		
					strings = [set.lower(),'s'+str(1),'f'+str(53)]
					fname = 'htpp_output\Job_notched45\BBindicator_' + '_'.join(strings) + '.dat'
					ind_notched = np.loadtxt(fname,delimiter=',')
	 
					strings = [set.lower(),'s'+str(1),'f'+str(63)]
					fname = 'htpp_output\Job_kim45\BBindicator_' + '_'.join(strings) + '.dat'
					ind_kim = np.loadtxt(fname,delimiter=',')     

					strings = [set.lower(),'s'+str(1),'f'+str(66)]
					fname = 'htpp_output\Job_jones45\BBindicator_' + '_'.join(strings) + '.dat'
					ind_jones = np.loadtxt(fname,delimiter=',')
	 
					strings = [set.lower(),'s'+str(1),'f'+str(66)]
					fname = 'htpp_output\Job_conde45\BBindicator_' + '_'.join(strings) + '.dat'
					ind_conde = np.loadtxt(fname,delimiter=',')
	 
					strings = [set.lower(),'s'+str(1),'f'+str(93)]
					fname = 'htpp_output\Job_topopt45\BBindicator_' + '_'.join(strings) + '.dat'
					ind_topopt = np.loadtxt(fname,delimiter=',')

					fig = htpp_plot_BBindicator_aux(ind_notched,ind_kim,ind_jones,ind_conde,ind_topopt,lims,txt,n, odb_name)
	 
					pdf.savefig(fig,dpi=500,bbox_inches='tight')

					table = htpp_plot_BBindicator_name(set,l,f,n,table)

	return pdf,table
# ------------------------------------------------------------------------------
