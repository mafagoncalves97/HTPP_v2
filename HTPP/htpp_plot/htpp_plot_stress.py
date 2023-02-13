# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
from re import L
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpt
from matplotlib.legend_handler import HandlerBase, HandlerPatch
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from math import floor, ceil
import math as math

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

# ------------------------------------------------------------- STRESS PLOT NAME
def htpp_plot_stress_name(set,step,frame,num,table):
	name = ['STRESS [']
	set = set.replace('_','\_')
	name.append('set: '  + set.upper() + ' /')
	name.append('step: ' + str(step) + ' /')
	name.append('frame: '+ str(frame))
	if num == 2:
		name.append('/ colormap: ROTATION')

	name.append(']')
	name = ' '.join(name)

	table.append(name)

	return table

# --------------------------------------------------------------- STRESS OPTIONS
def htpp_plot_stress_option(option):
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
	yls=1
	return step,frame,yls,elem_set,lims,txt,num

# ------------------------------------------------------------------- STRESS AUX
def htpp_plot_stress_aux(stress,rotAngle,peeq,lims,txt,n,odb_name,dir,f, peeq_path1, stress_path1,peeq_path2, stress_path2,peeq_path3, stress_path3,peeq_path4, stress_path4):
	s_elas = stress[peeq == 0.0,:]
	s_plas = stress[peeq != 0.0,:]
	p_elas = peeq[peeq == 0.0]
	p_plas = peeq[peeq != 0.0]

	if n == 2 or n == 4:
		rot_elas = rotAngle[peeq == 0.0]
		rot_plas = rotAngle[peeq != 0.0]

	fig1 = plt.figure(num=100)
	ax = fig1.add_subplot(111)



	if lims:
		xlims,ylims = lims[0],lims[-1]
	else:
		ylims = [floor((min(stress[:,0])/100))*100,ceil((max(stress[:,0])/100))*100]
		if abs(min(stress[:,1])) > abs(max(stress[:,1])):
			lim_X = floor(min(stress[:,1])/100)*100
			xlims = [lim_X,-lim_X]
		else:
			lim_X = ceil(max(stress[:,1])/100)*100
			xlims = [-lim_X,lim_X]

	ylims[0] = -1000
	ylims[1] = 1000
	xlims[0] = -1000
	xlims[1] = 1000

	tmpSh = max(ylims[1],abs(xlims[0]))
	shear = [[0,-tmpSh],[0,tmpSh]]
	tmpBi = max(ylims[1],xlims[1])
	eqBiax =[[-tmpBi,tmpBi],[-tmpBi,tmpBi]]

	fig = plt.figure(num=n)
	ax = fig.add_subplot(111)


	blue = (0,0.3,0.85,1.0)
	grey = (0.75,0.75,0.75,1.0)
	black_02 = (0,0,0,0.2)
	al = 'center'
	fnt = 8
	rast = True
	siz = 1.7

	box = dict(boxstyle='Square,pad=0.05',fc='w',ec='w')
	p = [0.85,0.95,0.6]
	if txt:
		p = [float(t) for t in txt]

	name = 'shear'
	plt.plot(shear[0],shear[1],'k-',lw=0.5,zorder=0)
	plt.text(-tmpSh*p[0],tmpSh*p[0],name,ha=al,va=al,fontsize=fnt,bbox=box)

	name = 'uniaxial tension'
	plt.plot([0,0],ylims,'k-',lw=0.5,zorder=0)
	plt.text(0,ylims[1]*p[1],name,ha=al,va=al,ma=al,fontsize=fnt,bbox=box)

	name = 'equibiaxial\ntension'
	plt.plot(eqBiax[0],eqBiax[1],'k-',lw=0.5,zorder=0)
	plt.text(700,700,name,ha=al,va=al,ma=al,fontsize=fnt,bbox=box)

	name = 'equibiaxial\ncompression'
	plt.text(-p[2]*tmpBi,-p[2]*tmpBi,name,ha=al,va=al,ma=al,fontsize=fnt,bbox=box)

	plt.plot([0,0],[ylims[0],ylims[1]],'k-',lw=0.5,zorder=0)
	name = 'uniaxial\ncompression'
	plt.text(-tmpSh*p[0]+30,0.0,name,ha=al,va=al,ma=al,fontsize=fnt,bbox=box)
	plt.plot([-1000,1000],[0,0],'k-',lw=0.5,zorder=0)

	lbs = ['elastic','plastic']
	if n == 1:
		pl = plt.scatter(s_plas[:,1],s_plas[:,0],s=siz,color=blue, edgecolors=black_02,linewidths=0.1,zorder=200,rasterized=rast)
		el = plt.scatter(s_elas[:,1],s_elas[:,0],s=0.5,color=grey, edgecolors='none',linewidths=0.2,zorder=100,rasterized=rast)

		leg = ax.legend([el,pl],lbs,loc=4,frameon=False, fontsize=fnt,markerscale=3.0,handletextpad=0.1)
		leg.legendHandles[0].set_edgecolor(grey)
		leg.legendHandles[1].set_edgecolor(blue)

	elif n == 2:
		cjet = mpl.cm.get_cmap('jet',12)
		pl = plt.scatter(s_plas[:,1],s_plas[:,0],s=siz,c=rot_plas,cmap=cjet, edgecolors=black_02,linewidths=0.1,zorder=200,rasterized=rast, norm=mpl.colors.Normalize(vmin=0.0,vmax=90.0))
		tks = np.linspace(0.0,90.0,13)
		clb = plt.colorbar(ticks=tks)
		clb.solids.set_rasterized(True)
		clb.set_label(r'$\gamma$ [$^\circ$]', labelpad=-25, y=1.11, rotation=0)
		clb.set_ticklabels(['%.1f'%i for i in tks])

		el = plt.scatter(s_elas[:,1],s_elas[:,0],s=0.5,color=grey, edgecolors='none',linewidths=0.2,zorder=100,rasterized=rast)

		cmap = plt.cm.jet
		hdls = [mpt.Rectangle((0,0),0.5,1,rasterized=True) for _ in lbs]
		hdl_map = dict(zip(hdls,                     [HandlerColormap(cmap,12,i) for i in range(2)]))
		leg1 = ax.legend(handles=hdls,labels=lbs,handler_map=hdl_map, fontsize=fnt,loc=4,frameon=False)

		for text in leg1.get_texts():
			text.set_color("w")

		c = [mpt.Circle((0.5, 0.5),1.0,fc='none',ec='w',lw=4.0) for _ in lbs]
		hdl_map = {mpt.Circle:HandlerEllipse()}
		leg2 = plt.legend(c,lbs,handler_map=hdl_map,loc=4,frameon=False, fontsize=fnt)
		ax.add_artist(leg1)

	elif n == 3:
		size = np.linspace(0,5,12)

		cjet = mpl.cm.get_cmap('jet', 12)
		pl = plt.scatter(s_plas[:,1],s_plas[:,0],s=2.5,c=p_plas,cmap=cjet,zorder=200,rasterized=rast, norm=plt.Normalize(vmin=0,vmax=0.26))

		sm = plt.cm.ScalarMappable(cmap=cjet, norm=plt.Normalize(vmin=0,vmax=0.26))
		sm._A = []
		tks = np.linspace(0.0,0.26,13)

		clb = plt.colorbar(sm,ticks=tks, orientation='vertical')
		clb.solids.set_rasterized(True)
		clb.set_label(r'$\bar \epsilon ^\mathrm{p}$ [-]',labelpad=-23,y=1.11, rotation=0)
		clb.set_ticklabels(['%.2f'%round(i,2) for i in tks])

		el = plt.scatter(s_elas[:,1],s_elas[:,0],s=3,color=grey, edgecolors='none',linewidths=0.2,zorder=100,rasterized=rast)

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

	elif n == 4:
		A = np.array([s_plas[:,0],s_plas[:,1],rot_plas,p_plas])
		B = A[:,np.argsort(A[0,:])]

		c0 = []
		for i in range(len(c1)):
			c0.append(int(np.sum(c1[0:i])))

		cjet = mpl.cm.get_cmap('jet', 12)
		tot = 0
		stk = np.linspace(0.1,8.0,12)
		for i in range(0,len(c0)-1):
			x = B[1,c0[i]:c0[i+1]]
			y = B[0,c0[i]:c0[i+1]]
			pl = plt.scatter(x,y,s=stk,c=B[2,c0[i]:c0[i+1]], cmap=cjet,edgecolors=black_02,linewidths=0.2, zorder=200,rasterized=True, norm=mpl.colors.Normalize(vmin=0.0,vmax=90.0))

		tks = np.linspace(0.0,90.0,13)
		clb = plt.colorbar(ticks=tks)
		clb.solids.set_rasterized(rast)
		clb.set_label(r'$\gamma$ [$^\circ$]', labelpad=-25, y=1.11, rotation=0)
		clb.set_ticklabels(['%.1f'%i for i in tks])

		el = plt.scatter(s_elas[:,1],s_elas[:,0],s=stk[0],color=grey, edgecolors='none',linewidths=0.2,zorder=100,rasterized=rast)

		cmap = plt.cm.jet
		hdls = [mpt.Rectangle((0,0),0.5,1,rasterized=True) for _ in lbs]
		hdl_map = dict(zip(hdls,                     [HandlerColormap(cmap,12,i) for i in range(2)]))
		leg1 = ax.legend(handles=hdls,labels=lbs,handler_map=hdl_map, fontsize=fnt,loc=4,frameon=False)

		for text in leg1.get_texts():
			text.set_color("w")

		c = [mpt.Circle((0.5, 0.5),1.0,fc='none',ec='w',lw=4.0) for _ in lbs]
		hdl_map = {mpt.Circle:HandlerEllipse()}
		leg2 = plt.legend(c,lbs,handler_map=hdl_map,loc=4,frameon=False, fontsize=fnt)
		ax.add_artist(leg1)

	for k, spine in ax.spines.items():
		spine.set_zorder(1000)
  
  
	####################################################################################################### PLOT YIELD LOCUS
	
	x1=np.arange(-1000,1100, 50)
	x2=np.arange(-1000,1100, 50)
	xy=np.zeros((len(x1),1))
	
	yyi=np.zeros((len(x1),len(x1)))
	yyf=np.zeros((len(x1),len(x1)))
	h=0
	l=0
	for i in range(len(x1)):
		for j in range(len(x1)):
			yyi[i,j],sigma_yi = Yld2000(x1[i], x2[j],xy[h],np.zeros(len(peeq)))
			yyf[i,j], sigma_yf = Yld2000(x1[i], x2[j],xy[h],peeq)
			l=l + 1
		h = h+1

	x1m, x2m = np.meshgrid(x1, x2)
	plt.contour(x1m,x2m,yyi, 0,linestyles='-',linewidths = 0.7, colors = 'k')
	plt.text(120,-350,r'$\sigma_\mathrm{y0}=%d\ \mathrm{MPa}$'%sigma_yi,fontsize=8,rotation=42)
	plt.contour(x1m,x2m,yyf, 0,linestyles='-',linewidths = 0.7, colors = 'k')
	plt.text(320,-640,r'$\sigma_\mathrm{y}=%d\ \mathrm{MPa}$'%sigma_yf,fontsize=8,rotation=42)
 
	plt.xlabel('$\sigma_2$ [MPa]')
	plt.ylabel('$\sigma_1$ [MPa]')
	plt.ylim(-1000,1000)
	plt.xlim(-1000,1000)
	
	# plt.yticks(np.arange(ylims[0],ylims[1]+200.0,200.0))
	# plt.xticks(np.arange(xlims[0],xlims[1]+200.0,200.0))
	#plt.tight_layout()
	#plt.gcf().set_size_inches(7,5)
	
	#plt.savefig(r'htpp_output/Figures/stress_diagram_bar_'+odb_name+'.jpg',bbox_inches='tight')
	plt.show()
	#clb.remove()
	#plt.savefig(r'htpp_output/Figures/stress_diagram_'+odb_name+'.jpg',bbox_inches='tight')
	plt.close()

	####################################################################################################### STRAIN PATH
 
	s_elas1 = stress_path1[peeq_path1 == 0.0,:]
	s_plas1 = stress_path1[peeq_path1 != 0.0,:]
	p_elas1 = peeq_path1[peeq_path1 == 0.0]
	p_plas1 = peeq_path1[peeq_path1 != 0.0]
 
	s_elas2 = stress_path2[peeq_path2 == 0.0,:]
	s_plas2 = stress_path2[peeq_path2 != 0.0,:]
	p_elas2 = peeq_path2[peeq_path2 == 0.0]
	p_plas2 = peeq_path2[peeq_path2 != 0.0]
 
	s_elas3 = stress_path3[peeq_path3 == 0.0,:]
	s_plas3 = stress_path3[peeq_path3 != 0.0,:]
	p_elas3 = peeq_path3[peeq_path3 == 0.0]
	p_plas3 = peeq_path3[peeq_path3 != 0.0]
 
	s_elas4 = stress_path4[peeq_path4 == 0.0,:]
	s_plas4 = stress_path4[peeq_path4 != 0.0,:]
	p_elas4 = peeq_path4[peeq_path4 == 0.0]
	p_plas4 = peeq_path4[peeq_path4 != 0.0]


	fig2 = plt.figure(figsize=(8,5))
	#ax1 = fig2.add_subplot(111)


	ylims[0] = -650
	ylims[1] = 650
	xlims[0] = -650
	xlims[1] = 650

	tmpSh = max(ylims[1],abs(xlims[0]))
	shear = [[0,-tmpSh],[0,tmpSh]]
	tmpBi = max(ylims[1],xlims[1])
	eqBiax =[[-tmpBi,tmpBi],[-tmpBi,tmpBi]]



	blue = (0,0.3,0.85,1.0)
	grey = (0.75,0.75,0.75,1.0)
	black_02 = (0,0,0,0.2)
	al = 'center'
	fnt = 8
	rast = True
	siz = 1.7

	box = dict(boxstyle='Square,pad=0.05',fc='w',ec='w')
	p = [0.85,0.95,0.6]
	if txt:
		p = [float(t) for t in txt]

	name = 'shear'
	plt.plot(shear[0],shear[1],'k-',lw=0.5,zorder=0)
	plt.text(-tmpSh*p[0],tmpSh*p[0],name,ha=al,va=al,fontsize=fnt,bbox=box)

	name = 'uniaxial tension'
	plt.plot([0,0],ylims,'k-',lw=0.5,zorder=0)
	plt.text(0,ylims[1]*p[1],name,ha=al,va=al,ma=al,fontsize=fnt,bbox=box)

	name = 'equibiaxial\ntension'
	plt.plot(eqBiax[0],eqBiax[1],'k-',lw=0.5,zorder=0)
	plt.text(520,600,name,ha=al,va=al,ma=al,fontsize=fnt,bbox=box)

	name = 'equibiaxial\ncompression'
	plt.text(-p[2]*tmpBi-50,-p[2]*tmpBi,name,ha=al,va=al,ma=al,fontsize=fnt,bbox=box)

	plt.plot([0,0],[ylims[0],ylims[1]],'k-',lw=0.5,zorder=0)
	name = 'uniaxial\ncompression'
	plt.text(-tmpSh*p[0]+30,0.0,name,ha=al,va=al,ma=al,fontsize=fnt,bbox=box)
	plt.plot([-1000,1000],[0,0],'k-',lw=0.5,zorder=0)

	lbs = ['elastic','plastic']
 
	cjet = mpl.cm.get_cmap('jet',12)
	#el1 = plt.scatter(s_elas1[:,1],s_elas1[:,0],s=15,color=grey, marker='o',edgecolors='none',zorder=100,rasterized=rast)
	#pl1 = plt.scatter(s_plas1[:,1],s_plas1[:,0],s=5,c=p_plas1,marker='o',cmap=cjet, edgecolors=black_02,zorder=200,rasterized=rast, norm=mpl.colors.Normalize(vmin=0.0,vmax=0.2))
	el2 = plt.scatter(s_elas2[:,1],s_elas2[:,0],s=10,color=grey, marker='^',edgecolors='none',zorder=100,rasterized=rast)
	pl2 = plt.scatter(s_plas2[:,1],s_plas2[:,0],s=10,c=p_plas2, marker='^',cmap=cjet,zorder=200,rasterized=rast, norm=mpl.colors.Normalize(vmin=0.0,vmax=0.1))
	el3 = plt.scatter(s_elas3[:,1],s_elas3[:,0],s=10,color=grey, marker='+', edgecolors='none',zorder=100,rasterized=rast)
	pl3 = plt.scatter(s_plas3[:,1],s_plas3[:,0],s=10,c=p_plas3, marker='+',cmap=cjet,zorder=200,rasterized=rast, norm=mpl.colors.Normalize(vmin=0.0,vmax=0.1))
	el4 = plt.scatter(s_elas4[:,1],s_elas4[:,0],s=10,color=grey, marker='x', edgecolors='none',zorder=100,rasterized=rast)
	pl4 = plt.scatter(s_plas4[:,1],s_plas4[:,0],s=10,c=p_plas4, marker='x',cmap=cjet,zorder=200,rasterized=rast, norm=mpl.colors.Normalize(vmin=0.0,vmax=0.1))
 
	tks = np.linspace(0.0,0.1,13)
	clb = plt.colorbar(ticks=tks)
	clb.solids.set_rasterized(True)
	clb.set_label(r'$\bar \varepsilon_\mathrm{p}$ [-]', labelpad=-25, y=1.11, rotation=0)
	clb.set_ticklabels(['%.2f'%i for i in tks])

	plt.contour(x1m,x2m,yyi, 0,linestyles='-',linewidths = 0.7, colors = 'k')
	plt.text(120,-350,r'$\sigma_\mathrm{y0}=%d\ \mathrm{MPa}$'%sigma_yi,fontsize=8,rotation=42)
	#plt.contour(x1m,x2m,yyf, 0,linestyles='-',linewidths = 0.7, colors = 'k')
	#plt.text(320,-640,r'$\sigma_\mathrm{y}=%d\ \mathrm{MPa}$'%sigma_yf,fontsize=8,rotation=42)
 
	plt.xlabel('$\sigma_2$ [MPa]')
	plt.ylabel('$\sigma_1$ [MPa]')
	plt.ylim(-650,650)
	plt.xlim(-650,650)
	


	
	plt.close()
 
 
 
 
 
 

	######################################################################################################## SIGMAXX VS SIGMAYY 

	
	'''fig2 = plt.figure()

	s_plas_1 = np.array([])
	s_plas_2 = np.array([])
	p_plas = np.array([])
	for ls in range(f,f+1):
		strings = ['stressxxyyxy','set-1','s'+str(1),'f'+str(ls)]
		fname_stress_ls = dir+'_'.join(strings)
		strings = ['peeq','set-1','s'+str(1),'f'+str(ls)]
		fname_peeq_ls = dir+'_'.join(strings)
		stress_ls = np.loadtxt(fname_stress_ls+'.dat',skiprows=1,delimiter=',')
		print((len(stress_ls)-1)/2)
		stress_ls = stress_ls[0:7852]
		peeq_ls = np.loadtxt(fname_peeq_ls+'.dat',delimiter=',')
  
		#np.append(p_elas, peeq_ls[peeq_ls == 0.0])
		p_plas=np.append(p_plas,peeq_ls[peeq_ls != 0.0])
  
		s_plas_1=np.append(s_plas_1,stress_ls[peeq_ls != 0.0,0])
		s_plas_2=np.append(s_plas_2,stress_ls[peeq_ls != 0.0,1])

	cjet = mpl.cm.get_cmap('jet', 12)
	pl = plt.scatter(s_plas_2,s_plas_1,s=2.5,c=p_plas,cmap=cjet,rasterized=rast)

		
		


		#el = plt.scatter(s_elas[:,1],s_elas[:,0],s=3,color=grey, edgecolors='none',linewidths=0.2,zorder=100,rasterized=rast)
	tks = np.linspace(0.0,0.26,13)
	sm = plt.cm.ScalarMappable(cmap=cjet, norm=plt.Normalize(vmin=0,vmax=0.26))
	sm._A = []
	clb = plt.colorbar(sm,ticks=tks, orientation='vertical')
	clb.solids.set_rasterized(True)
	clb.set_label(r'$\bar \epsilon ^\mathrm{p}$ [-]',labelpad=-23,y=1.11, rotation=0)
	clb.set_ticklabels(['%.2f'%round(i,2) for i in tks])

	#plt.savefig(r'htpp_output/Figures/stress_diagram_total_'+odb_name+'.jpg',bbox_inches='tight')
	
	
	x1=np.arange(-1000,1100, 50)
	x2=np.arange(-1000,1100, 50)
	xy=0.1+np.zeros((len(x1),1))
	
	yyi=np.zeros((len(x1),len(x1)))
	yyf=np.zeros((len(x1),len(x1)))
	h=0
	l=0
	for i in range(len(x1)):
		for j in range(len(x1)):
			yyi[i,j],sigma_yi = Yld2000(x1[i], x2[j],xy[h],np.zeros(len(peeq)))
			yyf[i,j], sigma_yf = Yld2000(x1[i], x2[j],xy[h],peeq)
			l=l + 1
		h = h+1
	x1m, x2m = np.meshgrid(x1, x2)
	plt.contour(x1m,x2m,yyi, 0,linestyles='-',linewidths = 0.7, colors = 'k')
	plt.text(120,-350,r'$\sigma_\mathrm{y0}=%d\ \mathrm{MPa}$'%sigma_yi,fontsize=8,rotation=42)
	plt.contour(x1m,x2m,yyf, 0,linestyles='-',linewidths = 0.7, colors = 'k')
	plt.text(300,-550,r'$\sigma_\mathrm{y}=%d\ \mathrm{MPa}$'%sigma_yf,fontsize=8,rotation=42)
 
	plt.xlabel('$\sigma_2$ [MPa]')
	plt.ylabel('$\sigma_1$ [MPa]')
	plt.ylim(-1000,1000)
	plt.xlim(-1000,1000)
	plt.show()
	plt.close()'''

 
	return fig1


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------ YIELD LOCUS 
# ------------------------------------------------------------------------------
def Yld2000(sig_XX, sig_YY, xy,peeq):
	am = np.zeros((2,3,3))
	a = [1.011,0.964,1.191,0.995, 1.010, 1.018, 0.977, 0.935]
	em = 6
	s = [sig_XX, sig_YY,xy]
	am[0,0,0] =  2.0*a[0]
	am[0,0,1] = -1.0*a[0]
	am[0,1,0] = -1.0*a[1]
	am[0,1,1] =  2.0*a[1]
	am[0,2,2] =  3.0*a[6]

	am[1,0,0] = -2.0*a[2]+2.0*a[3]+8.0*a[4]-2.0*a[5]
	am[1,0,1] =      a[2]-4.0*a[3]-4.0*a[4]+4.0*a[5]
	am[1,1,0] =  4.0*a[2]-4.0*a[3]-4.0*a[4]+    a[5]
	am[1,1,1] = -2.0*a[2]+8.0*a[3]+2.0*a[4]-2.0*a[5]
	am[1,2,2] =  9.0*a[7]

	am[0] = am[0]/3.0
	am[1] = am[1]/9.0

	p = [1.0,-1.0]

	y = np.zeros((2,3))
	for n in range(2):
		y[n] = am[n].dot(s)

	x = np.zeros((2,2))
	for n in range(2):
		a = math.sqrt((y[n,0]-y[n,1])**2 + 4*y[n,2]**2)
		for i in range(2):
			x[n,i] = 0.5*(y[n,0]+y[n,1]+p[i]*a)

	phi = np.zeros(2)
	phi[0] = abs(x[0,0]-x[0,1])**em
	phi[1] = abs(2*x[1,1]+x[1,0])**em + abs(2*x[1,0]+x[1,1])**em

	q = phi[0] + phi[1]
	if q <= 0.0:
		q = 0.0
	if np.sum(peeq)>0:
		sigma_y = 979.46*(max(peeq))**0.194
		se = ((0.5*q)**(1.0/em)) - sigma_y
	else:
		sigma_y = 979.46*(0.00535)**0.194
		se = ((0.5*q)**(1.0/em)) - sigma_y


	return se, sigma_y


# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRESS
# ------------------------------------------------------------------------------
def htpp_plot_stress(odb_name,option,pdf,size,table):

	size = [[[-0.1,-0.2],[6.0,4.6]],[[-0.2,-0.1],[5.725,4.675]],              [[-0.2,-0.1],[6.625,4.675]],[[-0.2,-0.1],[5.725,4.675]]]

	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif',size=14)

	dir = 'htpp_output/' + odb_name + '/'

	step,frame,yls,elem_set,lims,txt,num = htpp_plot_stress_option(option)

	for set in elem_set:
		for l in step:
			for f in frame:
				for n in num:
					strings = ['stress',set.lower(),'s'+str(l),'f'+str(f)]
					fname = dir + '_'.join(strings) + '.dat'
					stress = np.loadtxt(fname,skiprows=1,delimiter=',')

					strings = ['peeq',set.lower(),'s'+str(l),'f'+str(f)]
					fname = dir + '_'.join(strings) + '.dat'
					peeq = np.loadtxt(fname)				

					###### PEEQ and STRESS PATHS
					strings = 'path_peeq_1565'
					fname = dir + strings + '.dat'
					peeq_path1 = np.loadtxt(fname,skiprows=1,delimiter=',')

					strings = 'path_stress_1565'
					fname = dir + strings + '.dat'
					stress_path1 = np.loadtxt(fname, skiprows=1,delimiter=',')
     
					strings = 'path_peeq_24768'
					fname = dir + strings + '.dat'
					peeq_path2 = np.loadtxt(fname,skiprows=1,delimiter=',')

					strings = 'path_stress_24768'
					fname = dir + strings + '.dat'
					stress_path2 = np.loadtxt(fname, skiprows=1,delimiter=',')
     
					strings = 'path_peeq_44050'
					fname = dir + strings + '.dat'
					peeq_path3 = np.loadtxt(fname,skiprows=1,delimiter=',')

					strings = 'path_stress_44050'
					fname = dir + strings + '.dat'
					stress_path3 = np.loadtxt(fname, skiprows=1,delimiter=',')
     
					strings = 'path_peeq_36235'
					fname = dir + strings + '.dat'
					peeq_path4 = np.loadtxt(fname,skiprows=1,delimiter=',')

					strings = 'path_stress_36235'
					fname = dir + strings + '.dat'
					stress_path4 = np.loadtxt(fname, skiprows=1,delimiter=',')
     
     
					rotAngle = None
					if n == 2 or n == 4:
						strings = ['rotation',set.lower(),'s'+str(l),'f'+str(f)]
						fname = dir + '_'.join(strings) + '.dat'
						rotAngle = np.loadtxt(fname)

					lims=[[-1000,1000],[-1000,1000]]

					fig1= htpp_plot_stress_aux(stress,rotAngle,peeq,lims,txt,n,odb_name,dir,f, peeq_path1, stress_path1,peeq_path2, stress_path2,peeq_path3, stress_path3,peeq_path4, stress_path4)

					#pdf.savefig(fig1,dpi=500,bbox_inches='tight')
					#pdf.savefig(fig2,dpi=500,bbox_inches='tight')

					table = htpp_plot_stress_name(set,l,f,n,table)

	return pdf,table
# ------------------------------------------------------------------------------
