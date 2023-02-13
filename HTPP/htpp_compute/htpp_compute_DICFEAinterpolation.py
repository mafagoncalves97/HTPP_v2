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
import csv
import numpy
import scipy.interpolate as inp 
from pandas import *
import math
from more_itertools import locate


# --------------------------------------------------------------- STRAIN OPTIONS


def htpp_compute_DICFEAinterpolation_option(option):
	step = '\n'.join(s for s in option if 'step=' in s).replace(",", "").replace("\n", "")[6:-1].strip().split(';')
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
	
	return step, frame, elem_set, lims, txt

# ------------------------------------------------------------------- STRAIN AUX


def htpp_compute_DICFEAinterpolation_aux(coord_fea, strain_fea, strain_dic, coord_dic):

	DIC_strains = np.empty((len(strain_dic),6))
	DIC_strains[:,0] = strain_dic[:,2] #exx
	DIC_strains[:,1] = strain_dic[:,3] #eyy
	DIC_strains[:,2] = strain_dic[:,4] #exy
	DIC_strains[:,3] = strain_dic[:,0] #X
	DIC_strains[:,4] = strain_dic[:,1] #Y
	
	DIC_coords = np.empty((len(coord_dic),4))
	DIC_coords[:,0] = (coord_dic[:,2]-coord_dic[:,4]) * 0.057845 #X 
	DIC_coords[:,1] = (coord_dic[:,3]-coord_dic[:,5]) * 0.057845 #Y
	DIC_coords[:,2] = coord_dic[:,0] #x_pix
	DIC_coords[:,3] = coord_dic[:,1] #y_pic
	print(max(DIC_coords[:,1]),min(DIC_coords[:,1]))
	#indexC = [np.any(i) for i in np.isnan(DIC_coords)]
	#DIC_coords = np.delete(DIC_coords,indexC,axis=0)
	#These two introduced for zones after replacing ;;; for zeros in excel
	#indexZ = [np.any(i) for i in DIC_coords==0]
	#indexS = [np.any(i) for i in np.isnan(DIC_strains)]
	#DIC_strains = np.delete(DIC_strains,indexS,axis=0)
	#DIC_coords = np.delete(DIC_coords,indexZ,axis=0)
	#DIC_coords = np.delete(DIC_coords,indexS,axis=0)
	#DIC_coords = DIC_coords[:96492]
	#DIC_coords[:,1] =  - DIC_coords[:,1]
	#indexS = [np.any(i) for i in np.isnan(DIC_strains)]
	#DIC_strains = np.delete(DIC_strains,indexS,axis=0)
	#indexSS = [np.any(i) for i in np.isnan(DIC_coords)]
	#DIC_coords = np.delete(DIC_coords,indexSS,axis=0)
	#DIC_strains = np.delete(DIC_strains,indexC,axis=0)



	
 
	table_var = np.empty((len(DIC_coords),3))
	table_var[:,0:2] = DIC_coords[:,2:4]
	table_var[:,2] = range(0,len(DIC_coords))

	table_strains_var = np.empty((len(DIC_strains),3))
	table_strains_var[:,0:2] = DIC_strains[:,3:5]
	table_strains_var[:,2] = range(0,len(DIC_strains))


	print(len(table_strains_var))
	print(len(table_var))
	for i in range(0,len(table_strains_var)):
		xx = np.where(table_var[:,0]==[table_strains_var[i,0]])
		yy = np.where(table_var[:,1]==[table_strains_var[i,1]])
		idx = np.intersect1d(xx,yy,i)
		if len(idx)==0:
			print(xx,yy,i,idx)
		table_strains_var[i,2]=table_var[idx[0],2]
		DIC_strains[i,5]=table_var[idx[0],2]


	#print(table_var)
	#print(table_strains_var)
	aux_Sxx = DIC_strains[:,0]
	aux_Syy = DIC_strains[:,1]
	aux_Sxy = DIC_strains[:,2]
	aux_X = DIC_strains[:,3]
	aux_Y = DIC_strains[:,4]

 
	arr1inds = DIC_strains[:,5].argsort()
	aux_Sxx = aux_Sxx[arr1inds]
	aux_Syy = aux_Syy[arr1inds]
	aux_Sxy = aux_Sxy[arr1inds]
	aux_X = aux_X[arr1inds]
	aux_Y = aux_Y[arr1inds]

	DIC_strains[:,0] = aux_Sxx
	DIC_strains[:,1] = aux_Syy
	DIC_strains[:,2] = aux_Sxy
	DIC_strains[:,3] = aux_X
	DIC_strains[:,4] = aux_Y
	#DIC_strains = np.delete(DIC_strains,indexZ,axis=0)
	print(DIC_strains,DIC_coords)

	DIC_data = np.empty((len(DIC_coords),5))
	DIC_data[:,0] = DIC_coords[:,0] -39.045
	DIC_data[:,1] = (-DIC_coords[:,1]) + 32.58
	DIC_data[:,2] = DIC_strains[:,0]
	DIC_data[:,3] = DIC_strains[:,1]
	DIC_data[:,4] = 2*DIC_strains[:,2]

	print(max(DIC_data[:,3]),min(DIC_data[:,3]))

	rmesh = coord_fea[:,0:2]


	#Determine the number of experimental increments/stage
	ndis = len(coord_dic)
	tsteps = 1
	
	#Allocate the receiver displacement mesh
	nrenod=len(rmesh)
	ruxuy=np.zeros((nrenod*tsteps,3))
	
 
	#Determine the pointers for receiver mesh - smoothing spatial variables
	#Mesh ordered by y-coordinate
	fcoord=rmesh[0,1]
	for i in range(0,len(rmesh)):
		if int(rmesh[i,1]) == int((fcoord)):
			colp = i - 1
			break
	rowp = int(nrenod/(colp))
	icolp = 0
	ecolp = colp

		
	#Determines the number of points in each stage

	ndonod = ndis

	#Initialize the indeces for receiver mesh
	dmeshX = DIC_data[:,0]
	dmeshY =  DIC_data[:,1]
	uxdo = DIC_data[:,2]
	uydo = DIC_data[:,3]
	uxydo = DIC_data[:,4]
	

	#Interpolation
	ruxuy[:,0]=inp.griddata((dmeshX,dmeshY), uxdo,(rmesh[:,0], rmesh[:,1]),'nearest')
	ruxuy[:,1]=inp.griddata((dmeshX,dmeshY), uydo,(rmesh[:,0], rmesh[:,1]),'nearest')
	ruxuy[:,2]=inp.griddata((dmeshX,dmeshY), uxydo,(rmesh[:,0], rmesh[:,1]),'nearest')



	#Apply spatial smoothing
	
	for k in range(0,rowp):
		ruxuy[icolp:ecolp,0] = np.convolve(ruxuy[icolp:ecolp,0],1,'valid')
		ruxuy[icolp:ecolp,1] = np.convolve(ruxuy[icolp:ecolp,1],1,'valid')
		ruxuy[icolp:ecolp,2] = np.convolve(ruxuy[icolp:ecolp,2],1,'valid')

		icolp = icolp+colp
		ecolp = ecolp+colp
	

	return ruxuy, DIC_data

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- STRAIN
# ------------------------------------------------------------------------------


def htpp_compute_DICFEAinterpolation(option,odb_name):
	

	dir = 'htpp_output/' + odb_name + '/'

	step, frame, elem_set, lims, txt = htpp_compute_DICFEAinterpolation_option(option)

	for set in elem_set:
		for l in step:
			for f in frame:
				strings = ['strainxxyyxy', set.lower(), 's'+str(l), 'f'+str(f)]
				fname = dir + '_'.join(strings) + '.dat'
				strain = np.loadtxt(fname, skiprows=1, delimiter=',')

				strings = ['COORD_undeformed']
				fname = dir + '_'.join(strings) + '.dat'
				coord = np.loadtxt(fname, skiprows=1, delimiter=',')
	
				fname = dir + 'Coords_U_DIC_Conde45_TimeStep_66.csv'
				coord_dic = np.loadtxt(fname, skiprows=1, delimiter=';')
	
				fname = dir + 'Strains_DIC_Conde45_TimeStep_66.csv'
				strain_dic = np.loadtxt(fname, skiprows=1, delimiter=';')
	
				vals, DIC_data = htpp_compute_DICFEAinterpolation_aux(coord, strain, strain_dic, coord_dic)
	
				strings = ['DIC_in_FEAmesh',set.lower(),'s'+str(l),'f'+str(f)]
				fname = dir+'_'.join(strings)+'.csv'
				with open(fname,'w',newline='') as file:
					csvwriter = csv.writer(file)
					csvwriter.writerows(vals) 

				strings = ['DIC_data',set.lower(),'s'+str(l),'f'+str(f)]
				fname = dir+'_'.join(strings)+'.csv'
				with open(fname,'w',newline='') as file:
					csvwriter = csv.writer(file)
					csvwriter.writerows(DIC_data) 

# ------------------------------------------------------------------------------
