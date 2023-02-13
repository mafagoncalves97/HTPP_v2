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
def htpp_plot_colormap_name(set,step,frame,variable,table):
	name = ['COLORMAP [']
	set = set.replace('_','\_')
	name.append('set: '  + set.upper() + ' /')
	name.append('step: ' + str(step) + ' /')
	name.append('frame: '+ str(frame))
	name.append('variable: '+ str(variable))


	name.append(']')
	name = ' '.join(name)

	table.append(name)

	return table

# --------------------------------------------------------------- STRAIN OPTIONS
def htpp_plot_colormap_option(option):
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
def htpp_plot_colormap_aux(coords,variable,peeq,lims,txt,v,odb_name):


	#s_elas = strain[peeq == 0.0,:]
	#s_plas = strain[peeq != 0.0,:]


	p_elas = variable[peeq == 0.0]
	p_plas = variable[peeq != 0.0]
 
	coords_elas = coords[peeq == 0.0]
	coords_plas = coords[peeq != 0.0]

	if lims:
		xlims,ylims = lims[0],lims[-1]
	else:
		ylims = [0.0,ceil((max(coords[:,0])*10))/10]
		if abs(min(coords[:,1])) > abs(max(coords[:,1])):
			lim_X = floor(min(coords[:,1])*10)/10
			xlims = [lim_X,-lim_X]
		else:
			lim_X = ceil(max(coords[:,1])*10)/10
			xlims = [-lim_X,-lim_X]

	#ylims = [-120,120]
	#xlims = [-170,170]


	fig = plt.figure()
	ax = fig.add_subplot(111)

	blue = (0,0.3,0.85,1.0)
	grey = (0.75,0.75,0.75,1.0)
	black_02 = (0,0,0,0.2)
	al = 'center'
	fnt = 8
	rast = True
	siz = 2

	box = dict(boxstyle='Square,pad=0.05',fc='w',ec='w')
	p = [0.8,0.85,0.85,0.95,0.6]
	if txt:
		p = [float(t) for t in txt]

	
	lbs = ['elastic','plastic']


	cjet = mpl.cm.get_cmap('jet', 12)
	plt.scatter(coords_plas[:,0],coords_plas[:,1],c=p_plas,s=2,cmap=cjet, zorder=100,rasterized=rast,norm=plt.Normalize(vmin=0,vmax=0.26))
	plt.scatter(coords_elas[:,0],coords_elas[:,1],c='grey', s=2, zorder=100,rasterized=rast)

	sm = plt.cm.ScalarMappable(cmap=cjet, norm=plt.Normalize(vmin=0,vmax=0.26))
	sm._A = []
	#tks = np.linspace(0.0,0.26,13)
	#clb = plt.colorbar(sm,ticks=tks, pad=-0.4)
	#clb.solids.set_rasterized(True)
	#clb.set_label(r'$\bar \varepsilon ^\mathrm{p} [-]$', rotation=-90, labelpad=27)
	#clb.set_ticklabels(['%.2f'%round(i,2) for i in tks])

	# cmap = plt.cm.jet
	# hdls = [mpt.Rectangle((0,0),0.5,1,rasterized=True) for _ in lbs]
	# hdl_map = dict(zip(hdls,                     [HandlerColormap(cmap,12,i) for i in range(2)]))
	# leg1 = ax.legend(handles=hdls,labels=lbs,handler_map=hdl_map, fontsize=fnt,loc=4,frameon=False)
	#
	# for text in leg1.get_texts():
	# 	text.set_color("w")
	#
	# c = [mpt.Circle((0.5, 0.5),1.0,fc='none',ec='w',lw=4.0) for _ in lbs]
	# hdl_map = {mpt.Circle:HandlerEllipse()}
	# leg2 = plt.legend(c,lbs,handler_map=hdl_map,loc=4,frameon=False, fontsize=fnt)
	# ax.add_artist(leg1)

	
	plt.ylim(ylims[0],ylims[1])
	plt.xlim(xlims[0],xlims[1])
	plt.yticks([])
	plt.xticks([])
	plt.axis('off')
	plt.tight_layout()
	plt.show()


	#plt.savefig(r'htpp_output/Figures/colormap_peeq_'+odb_name+'.jpg',bbox_inches='tight')
 
	plt.close()

	return fig

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRAIN
# ------------------------------------------------------------------------------
def htpp_plot_colormap(odb_name,option,pdf,size,table):
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif',size=14)

	dir = 'htpp_output/' + odb_name + '/'

	step,frame,elem_set,lims,txt = htpp_plot_colormap_option(option)

	for set in elem_set:
		for l in step:
			for f in frame:
				
				fname = dir + 'COORD_undeformed' + '.dat'
				coords = np.loadtxt(fname,skiprows=1,delimiter=',')

				
				strings = ['peeq',set.lower(),'s'+str(l),'f'+str(f)]
				fname = dir + '_'.join(strings) + '.dat'
				peeq = np.loadtxt(fname)
				
				v = 'peeq'
				strings = [v,set.lower(),'s'+str(l),'f'+str(f)]
				fname = dir + '_'.join(strings) + '.dat'
				peeq = np.loadtxt(fname)
	
				fig = htpp_plot_colormap_aux(coords,peeq,peeq,lims,txt,v, odb_name)

				pdf.savefig(fig,dpi=500,bbox_inches='tight')

				table = htpp_plot_colormap_name(set,l,f,v,table)

	return pdf,table
# ------------------------------------------------------------------------------
