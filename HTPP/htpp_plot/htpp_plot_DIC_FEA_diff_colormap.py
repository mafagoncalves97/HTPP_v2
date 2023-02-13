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

# -------------------------------------------------------------- HANDLER ELLIPSE
class HandlerEllipse(HandlerPatch):
	def create_artists(self,legend,orig_handle,xdescent,ydescent,width,  height,fontsize,trans):
		x, y = 23.5, 3.25
		wid = 15.95
		center = (x,y)
		p = mpt.Ellipse(xy=center,width=wid,height=wid)
		self.update_prop(p,orig_handle,legend)
		p.set_transform(trans)

		return [p]

# ------------------------------------------------------------ HANDLER COLOR MAP
class HandlerColormap(HandlerBase):
	def __init__(self,cmap,num,n,**kw):
		HandlerBase.__init__(self, **kw)
		self.cmap = cmap
		self.num = num
		self.n = n

	def create_artists(self,legend,orig_handle,xdescent,ydescent,width, height,fontsize,trans):
		grey = (0.75,0.75,0.75,1.0)
		stripes = []
		width = 12.0
		x,y = 17.5, -2.8
		if self.n == 0:
			self.num = 1
		for i in range(self.num):
			coord = [x+i*width/self.num,y]
			wid = width/self.num
			if self.n == 0:
				s = mpt.Rectangle(coord,wid,width,transform=trans,fc=grey)
			elif self.n == 1:
				s = mpt.Rectangle(coord,wid,width,transform=trans, fc=self.cmap((2*i+1)/(2*self.num)))

			stripes.append(s)

		return stripes

# ------------------------------------------------------------. STRAIN PLOT NAME
def htpp_plot_DIC_FEA_colormap_name(set,step,frame,table):
	name = ['COLORMAP [']
	set = set.replace('_','\_')
	name.append('set: '  + set.upper() + ' /')
	name.append('step: ' + str(step) + ' /')
	name.append('frame: '+ str(frame))


	name.append(']')
	name = ' '.join(name)

	table.append(name)

	return table

# --------------------------------------------------------------- STRAIN OPTIONS
def htpp_plot_DIC_FEA_colormap_option(option):
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

	return step,frame,elem_set,lims,txt

# ------------------------------------------------------------------- STRAIN AUX
def htpp_plot_DIC_FEA_colormap_aux(coord,FEA_strains,DIC_strains,DIC_data,lims,txt,odb_name):

	'''strains = np.empty((int(len(FEA_strains)/2),3))
	strains[:,0] = abs(FEA_strains[:int(len(FEA_strains)/2),0]-DIC_strains[:int(len(FEA_strains)/2),0])
	strains[:,1] = abs(FEA_strains[:int(len(FEA_strains)/2),1]-DIC_strains[:int(len(FEA_strains)/2),1])
	strains[:,2] = abs(FEA_strains[:int(len(FEA_strains)/2),2]-DIC_strains[:int(len(FEA_strains)/2),2])'''


	strains = np.empty((int(len(FEA_strains)/2),3))
	strains[:,0] = (abs(FEA_strains[:int(len(FEA_strains)/2),0]-DIC_strains[:int(len(FEA_strains)/2),0])/max(FEA_strains[:int(len(FEA_strains)/2),0]))*100
	strains[:,1] = (abs(FEA_strains[:int(len(FEA_strains)/2),1]-DIC_strains[:int(len(FEA_strains)/2),1])/max(FEA_strains[:int(len(FEA_strains)/2),1]))*100
	strains[:,2] = (abs(FEA_strains[:int(len(FEA_strains)/2),2]-DIC_strains[:int(len(FEA_strains)/2),2])/max(abs(FEA_strains[:int(len(FEA_strains)/2),2])))*100

	coords = np.empty((len(coord),2))
	coords[:,0] = coord[:,0]
	coords[:,1] = coord[:,1]

	'''coords = np.empty((len(DIC_data),2))
	coords[:,0] = DIC_data[:,0]
	coords[:,1] = DIC_data[:,1]'''

	rast = True
	fig, axes = plt.subplots(1,3,figsize=(8,3), dpi=150)
	plt.subplots_adjust(wspace=0.2)
 

	cjet = mpl.cm.get_cmap('jet', 12)
	axes[0].scatter(coords[:,0],coords[:,1],c=strains[:,0],s=1.2,cmap=cjet, zorder=100,rasterized=rast)
	axes[0].set_xticks([])
	axes[0].set_yticks([])
	axes[0].axis('off')
	sm1 = plt.cm.ScalarMappable(cmap=cjet, norm=plt.Normalize(vmin=0,vmax=30))
	sm1._A = []
	#tks = np.linspace(min(strains[:,0]),max(strains[:,0]),13)
	tks = [0,5,10,15,20,25,30]
	#clb1 = plt.colorbar(sm1,ticks=tks, ax=axes[0])
	#clb1.solids.set_rasterized(True)
	#clb1.set_label(r'$\varepsilon_{xx}^{FEA}-\varepsilon_{xx}^{DIC}$ [-]',labelpad=-23,y=1.11, rotation=0)
	#clb1.set_ticklabels(['%.2f'%round(i,2) for i in tks])



	axes[1].set_xticks([])
	axes[1].set_yticks([])
	axes[1].axis('off')
	#ax2.set_ylim(-60,60)
	#ax2.set_xlim(-20,20)
	cjet = mpl.cm.get_cmap('jet', 12)
	axes[1].scatter(coords[:,0],coords[:,1],c=strains[:,1],s=1.2,cmap=cjet, zorder=100,rasterized=rast)
	sm2 = plt.cm.ScalarMappable(cmap=cjet, norm=plt.Normalize(vmin=0,vmax=30))
	sm2._A = []
	#tks = np.linspace(min(strains[:,1]),max(strains[:,1]),13)
	tks = [0,5,10,15,20,25,30]
	#clb2 = plt.colorbar(sm2,ticks=tks, ax=axes[1])
	#clb2.solids.set_rasterized(True)
	#clb2.set_label(r'$\varepsilon_{yy}^{FEA}-\varepsilon_{yy}^{DIC}$ [-]',labelpad=-23,y=1.11, rotation=0)
	#clb2.set_ticklabels(['%.2f'%round(i,2) for i in tks])
 
 

	axes[2].set_xticks([])
	axes[2].set_yticks([])
	axes[2].axis('off')
	#ax3.set_ylim(-60,60)
	#ax3.set_xlim(-20,20)
	cjet = mpl.cm.get_cmap('jet', 12)
	axes[2].scatter(coords[:,0],coords[:,1],c=strains[:,2],s=1.2,cmap=cjet, zorder=100,rasterized=rast)
 
	sm3 = plt.cm.ScalarMappable(cmap=cjet, norm=plt.Normalize(vmin=0,vmax=30))
	sm3._A = []
	#tks = np.linspace(min(strains[:,2]),max(strains[:,2]),13)
	tks = [0,5,10,15,20,25,30]
	clb3 = fig.colorbar(sm3,ax=axes,orientation='vertical')
	clb3.solids.set_rasterized(True)
	clb3.set_label(r'Error FEA vs DIC [%]',labelpad=5,y=0.5, rotation=90,size=12)
	clb3.set_ticklabels(['%d'%round(i,6) for i in tks])

	#plt.tight_layout()
	plt.show()
	plt.close()

	err_strains = np.empty((int(len(FEA_strains)/2),3))
	err_strains[:,0] = np.mean((abs(FEA_strains[:int(len(FEA_strains)/2),0]-DIC_strains[:int(len(FEA_strains)/2),0])/max(FEA_strains[:int(len(FEA_strains)/2),0]))*100)
	err_strains[:,1] = np.mean((abs(FEA_strains[:int(len(FEA_strains)/2),1]-DIC_strains[:int(len(FEA_strains)/2),1])/max(FEA_strains[:int(len(FEA_strains)/2),1]))*100)
	err_strains[:,2] = np.mean((abs(FEA_strains[:int(len(FEA_strains)/2),2]-DIC_strains[:int(len(FEA_strains)/2),2])/max(abs(FEA_strains[:int(len(FEA_strains)/2),2])))*100)

	return fig, err_strains

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRAIN
# ------------------------------------------------------------------------------
def htpp_plot_DIC_FEA_diff_colormap(odb_name,option,pdf,size,table):
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif',size=12)

	dir = 'htpp_output/' + odb_name + '/'

	step,frame,elem_set,lims,txt = htpp_plot_DIC_FEA_colormap_option(option)

	for set in elem_set:
		for l in step:
			for f in frame:

				fname = dir + 'COORD_undeformed' + '.dat'
				coords = np.loadtxt(fname,skiprows=1,delimiter=',')

				
				strings = ['strainxxyyxy',set.lower(),'s'+str(l),'f'+str(f)]
				fname = dir + '_'.join(strings) + '.dat'
				FEA_strains = np.loadtxt(fname, skiprows=1, delimiter=',')
	
				strings = ['DIC_in_FEAmesh',set.lower(),'s'+str(l),'f'+str(f)]
				fname = dir + '_'.join(strings) + '.csv'
				strain_dic = np.loadtxt(fname, delimiter=',')

				strings = ['DIC_data',set.lower(),'s'+str(l),'f'+str(f)]
				fname = dir+'_'.join(strings)+'.csv'
				DIC_data = np.loadtxt(fname, delimiter=',')
	 

				fig,error = htpp_plot_DIC_FEA_colormap_aux(coords,FEA_strains,strain_dic,DIC_data,lims,txt,odb_name)

				strings = ['DIC_FEA_error']
				fname = dir+'_'.join(strings)+'.dat'
				with open(fname,'w') as fi:
					fi.write('%f,%f,%f\n' %(error[0,0],error[0,1],error[0,2]))

				pdf.savefig(fig,dpi=500)

				table = htpp_plot_DIC_FEA_colormap_name(set,l,f,table)

	return pdf,table
# ------------------------------------------------------------------------------
