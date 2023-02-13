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
from scipy.stats import norm

# ------------------------------------------------------------. STRAIN PLOT NAME
def htpp_plot_indicator_name(set,step,frame,num,table):
	name = ['INDICATOR [']
	set = set.replace('_','\_')
	name.append('set: '  + set.upper() + ' /')
	name.append('step: ' + str(step) + ' /')
	name.append('frame: '+ str(frame))

	name.append(']')
	name = ' '.join(name)

	table.append(name)

	return table

# --------------------------------------------------------------- STRAIN OPTIONS
def htpp_plot_indicator_option(option):
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
def htpp_plot_indicator_aux(ind_notched,ind_kim,ind_jones,ind_conde,ind_topopt,lims,txt,n,odb_name):

	
	fig = plt.figure()
	barWidth=0.8
	#plt.rc('text', usetex=1)
	plt.rc('font', family='serif',size=8)
	plt.rc('axes', labelsize=8)
	c = ['lightgrey',  'darkgrey', 'dimgrey']
	plt.subplots_adjust(wspace= 1, hspace= 0.4)

	sub1 = fig.add_subplot(2,6,1)
	sub1.set_title(r'Std($\varepsilon_2/\varepsilon_1$)')
	IT = [ind_notched[1],  ind_jones[1],ind_topopt[1]]
	sub1.set_xticks([])
	br1 = np.arange(len(IT))
	plt.bar(br1, IT, color = c, width = barWidth,edgecolor ='None')
		

	sub2 = fig.add_subplot(2,6,2)
	sub2.set_title(r'$(\varepsilon_2/\varepsilon_1)_\mathrm{R}$')
	IT = [ind_notched[2],   ind_jones[2], ind_topopt[2]]
	sub2.set_xticks([])
	br1 = np.arange(len(IT))
	plt.bar(br1, IT, color = c, width = barWidth,edgecolor ='None')

	sub3 = fig.add_subplot(2,6,3)
	sub3.set_title(r'Std($\bar \varepsilon^\mathrm{p}$)')
	IT = [ind_notched[3],   ind_jones[3],  ind_topopt[3]]
	sub3.set_xticks([])
	br1 = np.arange(len(IT))
	plt.bar(br1, IT, color = c, width = barWidth,edgecolor ='None')

	sub4 = fig.add_subplot(2,6,4)
	sub4.set_title(r'$\bar \varepsilon_\mathrm{max}^\mathrm{p}$')
	IT = [ind_notched[4],   ind_jones[4], ind_topopt[4]]
	sub4.set_xticks([])
	br1 = np.arange(len(IT))
	plt.bar(br1, IT, color = c, width = barWidth,edgecolor ='None')

	sub5 = fig.add_subplot(2,6,5)
	sub5.set_title(r'$\bar \varepsilon_\mathrm{tension}^\mathrm{p}$')
	IT = [ind_notched[6],   ind_jones[6], ind_topopt[6]]
	sub5.set_xticks([])
	br1 = np.arange(len(IT))
	plt.bar(br1, IT, color = c, width = barWidth,edgecolor ='None')

	sub6 = fig.add_subplot(2,6,6)
	sub6.set_title(r'$\bar \varepsilon_\mathrm{shear}^\mathrm{p}$')
	IT = [ind_notched[7], ind_jones[7],  ind_topopt[7]]
	sub6.set_xticks([])
	br1 = np.arange(len(IT))
	plt.bar(br1, IT, color = c, width = barWidth, edgecolor ='None')


	sub7 = fig.add_subplot(2,5,6)
	sub7.set_title(r'$\bar \varepsilon_\mathrm{biaxial}^\mathrm{p}$')
	IT = [ind_notched[8],  ind_jones[8], ind_topopt[8]]
	sub7.set_xticks([])
	br1 = np.arange(len(IT))
	plt.bar(br1, IT, color = c, width = barWidth,edgecolor ='None')


	sub8 = fig.add_subplot(2,5,7)
	sub8.set_title(r'$\bar \varepsilon_\mathrm{plane}^\mathrm{p}$')
	IT = [ind_notched[10],   ind_jones[10],  ind_topopt[10]]
	sub8.set_xticks([])
	#sub8.set_ylim(0.0,max(IT))
	br1 = np.arange(len(IT))
	plt.bar(br1, IT, color = c, width = barWidth,edgecolor ='None')


	sub9 = fig.add_subplot(2,5,8)
	sub9.set_title(r'$\bar \varepsilon_\mathrm{compression}^\mathrm{p}$')
	IT = [ind_notched[9], ind_jones[9],  ind_topopt[9]]
	sub9.set_xticks([])
	br1 = np.arange(len(IT))
	plt.bar(br1, IT, color = c, width = barWidth,edgecolor ='None')


	sub10= fig.add_subplot(2,5,9)
	sub10.set_title(r'Av$(\bar \varepsilon^\mathrm{p})$')
	IT = [ind_notched[5],   ind_jones[5],  ind_topopt[5]]
	sub10.set_xticks([])
	br1 = np.arange(len(IT))
	plt.bar(br1, IT, color = c, width = barWidth, edgecolor ='None')


	sub11 = fig.add_subplot(2,5,10)
	sub11.set_title(r'$I_\mathrm{t}$')
	IT = [ind_notched[0],  ind_jones[0],  ind_topopt[0]]
	sub11.set_xticks([])
	#sub11.set_ylim(0.0,max(IT))
	br1 = np.arange(len(IT))
	plt.bar(br1, IT, color = c, width = barWidth,edgecolor ='None', label =['Notched', 'D',  'TopOpt'])



	fig.legend(loc='lower center',ncols=5)
	plt.savefig(r'htpp_output/Figures/nelson_indicator_specimens.jpg')
	plt.show()
	return fig


def http_plot_E2E1_ratio(ratio_jones, ratio_conde,ratio_topopt,lims,txt,n,odb_name):
	
	E2E1_ratio_conde = ratio_conde
	E2E1_ratio_jones = ratio_jones
	E2E1_ratio_topopt = ratio_topopt

	E2E1_ratio_conde[np.isnan(E2E1_ratio_conde)==1]=0
	E2E1_ratio_conde[E2E1_ratio_conde < -15] = -15
	E2E1_ratio_conde[E2E1_ratio_conde > 1] = 1
 
	E2E1_ratio_jones[E2E1_ratio_jones < -15] = -15
	E2E1_ratio_jones[E2E1_ratio_jones > 1] = 1
 
	E2E1_ratio_topopt[E2E1_ratio_topopt < -15] = -15
	E2E1_ratio_topopt[E2E1_ratio_topopt > 1] = 1
 
 
	std_jones = statistics.stdev(E2E1_ratio_jones.tolist())  # STANDARD DEVIATION OF STRAIN STATE (dispersion of E2E1 ratio)
	u_jones = statistics.mean(E2E1_ratio_jones.tolist())  # STANDARD DEVIATION OF STRAIN STATE (dispersion of E2E1 ratio)
 
	std_conde = statistics.stdev(E2E1_ratio_conde.tolist())  # STANDARD DEVIATION OF STRAIN STATE (dispersion of E2E1 ratio)
	u_conde = statistics.mean(E2E1_ratio_conde.tolist())  # STANDARD DEVIATION OF STRAIN STATE (dispersion of E2E1 ratio)
 
	std_topopt = statistics.stdev(E2E1_ratio_topopt.tolist())
	u_topopt = statistics.mean(E2E1_ratio_topopt.tolist())  # STANDARD DEVIATION OF STRAIN STATE (dispersion of E2E1 ratio)
	print(std_jones, std_topopt)

	fig = plt.figure(figsize=(12,6))
	#plt.rc('text', usetex=1)
	plt.rc('font', family='serif',size=14)
	plt.rc('axes', labelsize=14)
	plt.subplots_adjust(wspace= 0.6, hspace= 0.4)

	sub1= fig.add_subplot(1,3,1)
	x = ratio_jones
	numBins = 50
	histj,bins = np.histogram(x, numBins)
	sub1.bar(bins.tolist()[:-1], histj.astype(np.float32)/sum(histj.tolist()),width=(bins.tolist()[-1]-bins.tolist()[-2]),color='grey',alpha=0.8, edgecolor='black', log=True, linewidth=0.2)
	#sub1.plot([-std_jones,-std_jones], [0, 4000], "k--")
	sub1.set_xlabel(r"$\varepsilon_2/\varepsilon_1$")
	sub1.set_ylabel(r"Elements distribution")
	#sub1.set_title(r"D specimen: $\sigma$={:.2f}".format(u_jones, std_jones))
	sub1.set_ylim([0, 1e0])
 
	sub2 = fig.add_subplot(1,3,2)
	xt = ratio_conde
	numBinst = 50
	histt,binst = np.histogram(xt, numBins)
	sub2.bar(binst.tolist()[:-1], histt.astype(np.float32)/sum(histt.tolist()),width=(bins.tolist()[-1]-bins.tolist()[-2]), color='grey',alpha=0.8, edgecolor='black', log=True, linewidth=0.2)
	#sub2.plot([-std_topopt,-std_topopt], [0, 20000], "k--")
	sub2.set_xlabel(r"$\varepsilon_2/\varepsilon_1$")
	sub2.set_ylabel(r"Elements distribution")
	#sub2.set_title(r"Topopt specimen: $\sigma$={:.2f}".format(u_topopt, std_topopt))
	sub2.set_ylim([0, 1e0])
 
	sub3 = fig.add_subplot(1,3,3)
	xt = ratio_topopt
	numBinst = 50
	histt,binst = np.histogram(xt, numBins)
	sub3.bar(binst.tolist()[:-1], histt.astype(np.float32)/sum(histt.tolist()),width=(bins.tolist()[-1]-bins.tolist()[-2]), color='grey',alpha=0.8, edgecolor='black', log=True, linewidth=0.2)
	#sub2.plot([-std_topopt,-std_topopt], [0, 20000], "k--")
	sub3.set_xlabel(r"$\varepsilon_2/\varepsilon_1$")
	sub3.set_ylabel(r"Elements distribution")
	#sub2.set_title(r"Topopt specimen: $\sigma$={:.2f}".format(u_topopt, std_topopt))
	sub3.set_ylim([0, 1e0])
 
	#plt.show()
 
 
# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRAIN
# ------------------------------------------------------------------------------
def htpp_plot_indicator(odb_name,option,pdf,size,table):
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif',size=14)

	dir = 'htpp_output/' + odb_name + '/'

	step,frame,elem_set,lims,txt,num = htpp_plot_indicator_option(option)

	for set in elem_set:
		for l in step:
			for f in frame:
				for n in num:
        
					strings = [set.lower(),'s'+str(1),'f'+str(53)]
					fname = 'htpp_output\Job_notched45\indicator_' + '_'.join(strings) + '.dat'
					ind_notched = np.loadtxt(fname,delimiter=',')
     
					strings = [set.lower(),'s'+str(1),'f'+str(63)]
					fname = 'htpp_output\Job_kim45\indicator_' + '_'.join(strings) + '.dat'
					ind_kim = np.loadtxt(fname,delimiter=',')     

					strings = [set.lower(),'s'+str(1),'f'+str(66)]
					fname = 'htpp_output\Job_jones45\indicator_' + '_'.join(strings) + '.dat'
					ind_jones = np.loadtxt(fname,delimiter=',')
     
					strings = [set.lower(),'s'+str(1),'f'+str(66)]
					fname = 'htpp_output\Job_conde45\indicator_' + '_'.join(strings) + '.dat'
					ind_conde = np.loadtxt(fname,delimiter=',')
     
					strings = [set.lower(),'s'+str(1),'f'+str(121)]
					fname = 'htpp_output\Job_topopt45\indicator_' + '_'.join(strings) + '.dat'
					ind_topopt = np.loadtxt(fname,delimiter=',')

					fig = htpp_plot_indicator_aux(ind_notched,ind_kim,ind_jones,ind_conde,ind_topopt,lims,txt,n, odb_name)
	 
					pdf.savefig(fig,dpi=500,bbox_inches='tight')

					table = htpp_plot_indicator_name(set,l,f,n,table)
     
					strings = [set.lower(),'s'+str(1),'f'+str(66)]
					fname = 'htpp_output\Job_jones45\\ratio_E2E1_' + '_'.join(strings) + '.dat'
					ratio_jones = np.loadtxt(fname,delimiter=',')

     
					strings = [set.lower(),'s'+str(1),'f'+str(66)]
					fname = 'htpp_output\Job_conde45\\ratio_E2E1_' + '_'.join(strings) + '.dat'
					ratio_conde = np.loadtxt(fname,delimiter=',')
     
     
					strings = [set.lower(),'s'+str(1),'f'+str(121)]
					fname = 'htpp_output\Job_topopt45\\ratio_E2E1_' + '_'.join(strings) + '.dat'
					ratio_topopt = np.loadtxt(fname,delimiter=',')

					http_plot_E2E1_ratio(ratio_jones,ratio_conde, ratio_topopt,lims,txt,n,odb_name)
	

	return pdf,table
# ------------------------------------------------------------------------------
