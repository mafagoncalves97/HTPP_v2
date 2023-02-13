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
def htpp_plot_identifiability_name(set,step,frame,num,table):
	name = ['IDENTIFIABILITY [']
	set = set.replace('_','\_')
	name.append('set: '  + set.upper() + ' /')
	name.append('step: ' + str(step) + ' /')
	name.append('frame: '+ str(frame))

	name.append(']')
	name = ' '.join(name)

	table.append(name)

	return table

# --------------------------------------------------------------- STRAIN OPTIONS
def htpp_plot_identifiability_option(option):
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
'''def htpp_plot_identifiability_diff(iden_notched,iden_kim,iden_jones,iden_conde,iden_topopt,lims,txt,n,odb_name):

	# ------------- graph properties ------------- 
	FS = 24
	plt.rcParams['axes.facecolor'] = (1, 1, 1)
	plt.rcParams['figure.facecolor'] = (1, 1, 1)
	plt.rcParams["font.family"] = "sans"
	plt.rcParams["font.serif"] = "Times New Roman"
	plt.rcParams['font.size'] = FS
	params = {"ytick.color" : (0, 0, 0),
		  "xtick.color" : (0, 0, 0),
		  "grid.color" : (0, 0 , 0),
		  "text.color" : (0, 0, 0),
		  "axes.labelcolor" : (0, 0, 0),
		  "axes.edgecolor" : (.15, .15, .15),
		  "text.usetex": False}
	plt.rcParams.update(params)   

	
	print("Plotting difference fields for step: " + str(Step))
# -------------Diff EpsX------------- 
	for i in range(len(Param)):
		fig2 = plt.figure(figsize=(11,7))
		ax = plt.axes([0.01, 0.01, .9, .95])

		plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Strain_xx_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'EpsXX_diff')
	
		plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
		ax.set_xticks([])
		ax.set_yticks([])
		ax.axis('off')
		plt.xlim(-20,80)
		plt.ylim(-20,80)

		cbar = plt.colorbar()

		plt.clim(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
		#sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
		tks = np.linspace(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)
		cbar.solids.set_rasterized(True)
		cbar.set_label(r'$\varepsilon_{xx}$',labelpad=-40,y=1.1, rotation=0)
		#cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
		plt.axis('off')
		cbar.solids.set(alpha=1)
		fig2.savefig(cwd + '\\Results\\Sensitivity_analysis\\EpsXXdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
		plt.close(fig2)
# -------------Diff Epsy------------- 
	for i in range(len(Param)):
		fig3 = plt.figure(figsize=(11,7))
		ax = plt.axes([0.01, 0.01, .9, .95])

		plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Strain_yy_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'EpsYY_diff')
	
		plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
		ax.set_xticks([])
		ax.set_yticks([])
		ax.axis('off')
		plt.xlim(-20,80)
		plt.ylim(-20,80)

		cbar = plt.colorbar()

		plt.clim(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
		#sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
		tks = np.linspace(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)       
		cbar.solids.set_rasterized(True)
		cbar.set_label(r'$\varepsilon_{yy}$',labelpad=-40,y=1.1, rotation=0)
		#cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
		plt.axis('off')
		cbar.solids.set(alpha=1)
		fig3.savefig(cwd + '\\Results\\Sensitivity_analysis\\EpsYYdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
		plt.close(fig3)
# -------------Diff Epsxy------------- 
	for i in range(len(Param)):
		fig4 = plt.figure(figsize=(11,7))
		ax = plt.axes([0.01, 0.01, .9, .95])

		plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Strain_xy_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'EpsYY_diff')
	
		plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
		ax.set_xticks([])
		ax.set_yticks([])
		ax.axis('off')
		plt.xlim(-20,80)
		plt.ylim(-20,80)

		cbar = plt.colorbar()

		plt.clim(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
		#sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
		tks = np.linspace(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)
		cbar.solids.set_rasterized(True)
		cbar.set_label(r'$\varepsilon_{xy}$',labelpad=-40,y=1.1, rotation=0)
		#cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
		plt.axis('off')
		cbar.solids.set(alpha=1)
		fig4.savefig(cwd + '\\Results\\Sensitivity_analysis\\EpsXYdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
		plt.close(fig4)
# -------------Diff SigmaX------------- 
	for i in range(len(Param)):
		fig5 = plt.figure(figsize=(11,7))
		ax = plt.axes([0.01, 0.01, .9, .95])

		plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Sigma_xx_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'SigmaXX_diff')
	
		plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
		ax.set_xticks([])
		ax.set_yticks([])
		ax.axis('off')
		plt.xlim(-20,80)
		plt.ylim(-20,80)

		cbar = plt.colorbar()

		plt.clim(0,max(abs(Sigma_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
		#sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
		tks = np.linspace(0,max(abs(Sigma_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)
		cbar.solids.set_rasterized(True)
		cbar.set_label(r'$\sigma_{xx}$',labelpad=-40,y=1.1, rotation=0)
		#cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
		plt.axis('off')
		cbar.solids.set(alpha=1)
		fig5.savefig(cwd + '\\Results\\Sensitivity_analysis\\SigmaXXdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
		plt.close(fig5)
# -------------Diff SigmaY------------- 
	for i in range(len(Param)):
		fig6 = plt.figure(figsize=(11,7))
		ax = plt.axes([0.01, 0.01, .9, .95])

		plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Sigma_yy_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'SigmaYY_diff')
	
		plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
		ax.set_xticks([])
		ax.set_yticks([])
		ax.axis('off')
		plt.xlim(-20,80)
		plt.ylim(-20,80)

		cbar = plt.colorbar()

		plt.clim(0,max(abs(Sigma_yy_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
		#sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
		tks = np.linspace(0,max(abs(Sigma_yy_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)
		cbar.solids.set_rasterized(True)
		cbar.set_label(r'$\sigma_{yy}$',labelpad=-40,y=1.1, rotation=0)
		#cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
		plt.axis('off')
		cbar.solids.set(alpha=1)
		fig6.savefig(cwd + '\\Results\\Sensitivity_analysis\\SigmaYYdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
		plt.close(fig6)
# -------------Diff SigmaXY------------- 
	for i in range(len(Param)):
		fig7 = plt.figure(figsize=(11,7))
		ax = plt.axes([0.01, 0.01, .9, .95])

		plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Sigma_xy_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'SigmaXY_diff')
	
		plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
		ax.set_xticks([])
		ax.set_yticks([])
		ax.axis('off')
		plt.xlim(-20,80)
		plt.ylim(-20,80)

		cbar = plt.colorbar()

		plt.clim(0,max(abs(Sigma_xy_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
		#sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
		tks = np.linspace(0,max(abs(Sigma_xy_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)
		cbar.solids.set_rasterized(True)
		cbar.set_label(r'$\sigma_{xy}$',labelpad=-40,y=1.1, rotation=0)
		#cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
		plt.axis('off')
		cbar.solids.set(alpha=1)
		fig7.savefig(cwd + '\\Results\\Sensitivity_analysis\\SigmaXYdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
		plt.close(fig7)'''

def htpp_plot_identifiability_index(Swift_notched, yld_notched, Swift_kim, yld_kim, Swift_jones, yld_jones, Swift_conde, yld_conde,Swift_topopt,yld_topopt,lims,txt,n,odb_name):

	fig, ax = plt.subplots(figsize=(12,3))


	# ----------------------------- gamma_K (collinearity index)

	# use LaTeX fonts in the plot
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	Best_solutions = list(reversed(['Notched', 'Sigma', 'D', 'Shape','TopOpt']))
	x_axis = np.arange(len(Best_solutions))
	I_K_swift = list(reversed([Swift_notched[1], Swift_kim[1], Swift_jones[1], Swift_conde[1], Swift_topopt[1]]))
	I_K_yld = list(reversed([yld_notched[1], yld_kim[1], yld_jones[1], yld_conde[1], yld_topopt[1]]))
	gamma_K_swift = list(reversed([Swift_notched[0], Swift_kim[0], Swift_jones[0], Swift_conde[0], Swift_topopt[0]]))
	gamma_K_yld = list(reversed([yld_notched[0], yld_kim[0], yld_jones[0], yld_conde[0], yld_topopt[0]]))

	plt.subplot(1,2,1)

	plt.barh(x_axis-0.2,gamma_K_swift, height=0.4, label = 'Swift', color=['lightgrey'])
	plt.barh(x_axis+0.2,gamma_K_yld, height=0.4, label = 'Yld2000-2d',color=['dimgrey'])
	plt.xlabel(r'Collinearity index ($\gamma_K$)')
	plt.yticks(x_axis, Best_solutions)
	plt.axvline(x=20, color='grey',linestyle='--')
	#plt.gca().xaxis.grid(True)
	plt.xlim([0, 100])
	plt.legend(ncols=2,loc=4,prop={'size': 8})
	plt.tight_layout()


	# ----------------------------- C_K plot (condition number)
	plt.subplot(1,2,2)

	plt.barh(x_axis-0.2,I_K_swift, height=0.4, label = 'Swift', color=['lightgrey'])
	plt.barh(x_axis+0.2,I_K_yld, height=0.4, label = 'Yld2000-2d',color=['dimgrey'])
	#plt.ylabel(r'\textbf{Analysed solutions}', fontsize=11)
	plt.yticks(x_axis, Best_solutions)
	plt.xlabel(r'Identifiability index ($I_K$)')
	plt.axvline(x=2, color='grey',linestyle='--')
	plt.axvline(x=3, color='grey',linestyle='--')
	#plt.gca().xaxis.grid(True)
	plt.xlim([0, 5])
	plt.legend(ncols=2,loc=4,prop={'size': 8})
	plt.tight_layout()

	# ----------------------------- IT1 plot


	plt.show()
	return fig

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRAIN
# ------------------------------------------------------------------------------
def htpp_plot_identifiability(odb_name,option,pdf,size,table):
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif',size=14)

	

	step,frame,elem_set,lims,txt,num = htpp_plot_identifiability_option(option)

	for set in elem_set:
		for l in step:
			for f in frame:
				for n in num:

					#LE and S fields differences
					strings = 'Swift_Collinearity_Identifiability_indexes_Job_notched45_'
					fname = 'htpp_output/' + 'Job_notched45' + '/' + strings + '.dat'
					Swift_notched = np.loadtxt(fname,delimiter=',')
	 
					strings = 'Yld2000_2d_Collinearity_Identifiability_indexes_Job_notched45_'
					fname = 'htpp_output/' + 'Job_notched45' + '/' + strings + '.dat'
					yld_notched = np.loadtxt(fname,delimiter=',')

					strings = 'Swift_Collinearity_Identifiability_indexes_Job_kim45_'
					fname = 'htpp_output/' + 'Job_kim45' + '/' + strings + '.dat'
					Swift_kim = np.loadtxt(fname,delimiter=',')
	 
					strings = 'Yld2000_2d_Collinearity_Identifiability_indexes_Job_kim45_'
					fname = 'htpp_output/' + 'Job_kim45' + '/' + strings + '.dat'
					yld_kim = np.loadtxt(fname,delimiter=',')

					strings = 'Swift_Collinearity_Identifiability_indexes_Job_jones45_'
					fname = 'htpp_output/' + 'Job_jones45' + '/' + strings + '.dat'
					Swift_jones = np.loadtxt(fname,delimiter=',')
	 
					strings = 'Yld2000_2d_Collinearity_Identifiability_indexes_Job_jones45_'
					fname = 'htpp_output/' + 'Job_jones45' + '/' + strings + '.dat'
					yld_jones = np.loadtxt(fname,delimiter=',')
	 
					strings = 'Swift_Collinearity_Identifiability_indexes_Job_conde45_'
					fname = 'htpp_output/' + 'Job_conde45' + '/' + strings + '.dat'
					Swift_conde = np.loadtxt(fname,delimiter=',')
	 
					strings = 'Yld2000_2d_Collinearity_Identifiability_indexes_Job_conde45_'
					fname = 'htpp_output/' + 'Job_conde45' + '/' + strings + '.dat'
					yld_conde = np.loadtxt(fname,delimiter=',')
	 
					strings = 'Swift_Collinearity_Identifiability_indexes_Job_topopt45_'
					fname = 'htpp_output/' + 'Job_topopt45' + '/' + strings + '.dat'
					Swift_topopt = np.loadtxt(fname,delimiter=',')

					strings = 'Yld2000_2d_Collinearity_Identifiability_indexes_Job_topopt45_'
					fname = 'htpp_output/' + 'Job_topopt45' + '/' + strings + '.dat'
					yld_topopt = np.loadtxt(fname,delimiter=',')

					
					fig = htpp_plot_identifiability_index(Swift_notched, yld_notched, Swift_kim, yld_kim, Swift_jones, yld_jones, Swift_conde, yld_conde,Swift_topopt, yld_topopt,lims,txt,n, odb_name)
					plt.savefig(r'htpp_output/Figures/identifiability_index_specimens.jpg')
					pdf.savefig(fig,dpi=500,bbox_inches='tight')

					#fig = htpp_plot_identifiability_diff(ident_diff_notched,ident_diff_notched,ident_diff_notched,ident_diff_notched,ident_diff_notched,lims,txt,n, odb_name)
					#plt.savefig(r'htpp_output/Figures/identifiability_diff_' + odb_name + '.jpg')
					#pdf.savefig(fig,dpi=500,bbox_inches='tight')

					

					table = htpp_plot_identifiability_name(set,l,f,n,table)

	return pdf,table
# ------------------------------------------------------------------------------
