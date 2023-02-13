# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
import os
import shutil
import sys
from math import sin, cos, pi, degrees
import numpy as np
from scipy.optimize import fsolve

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------- IMPORT FILES
# -------------------------------------------------------------- Yield Functions
from htpp_yfunc_mises import htpp_yfunc_mises, htpp_yfunc_mises_param
from htpp_yfunc_hill48 import htpp_yfunc_ylocus_hill48, htpp_yfunc_hill48_param
from htpp_yfunc_yld2004 import htpp_yfunc_ylocus_yld2004, htpp_yfunc_yld2004_param
from htpp_yfunc_yld2000_2d import htpp_yfunc_ylocus_yld2000_2d, htpp_yfunc_yld2000_2d_param

# -------------------------------------------------------------- YIELD LOCUS AUX
def htpp_yfunc_yieldlocus_aux_ori(ang,yldstress,s12):
	xx = yldstress*cos(ang)
	yy = yldstress*sin(ang)
	xy = yldstress*s12
	s = [xx, yy, 0, xy, 0, 0]

	return s,xx,yy

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------ YIELD LOCUS
# ------------------------------------------------------------------------------
def htpp_yfunc_yieldlocus_ori(option,odb_name):
	yldname = '\n'.join(s for s in option if 'yldname=' in s).replace(",", "").replace("\n", "")[9:-1].strip()
	shear = '\n'.join(s for s in option if 'shear=' in s).replace(",", "").replace("\n", "")[7:-1].split(';')
	shear = [float(s) for s in shear]

	dir = 'htpp_output/' + odb_name + '/'

	npts = 500
	angle = np.linspace(0.0,2*pi,npts)
	alpha = 1

	for s12 in shear:
		xYld,yYld,xyYld = [],[],[]
#	---------------------------------------------------------------------- mises
		if yldname == 'mises':
			s0 = htpp_yfunc_mises_param(odb_name)
			for i in range(0,len(angle)):
					s,xx,yy = htpp_yfunc_yieldlocus_aux_ori(angle[i],s0,s12)
					alpha = fsolve(htpp_yfunc_mises,1.0,args=(s,s0))
					xYld.append(xx*alpha)
					yYld.append(yy*alpha)

#	--------------------------------------------------------------------- hill48
		elif yldname == 'hill48':
			c,s0 = htpp_yfunc_hill48_param(odb_name)
			for i in range(0,len(angle)):
					s,xx,yy = htpp_yfunc_yieldlocus_aux_ori(angle[i],s0,s12)
					alpha = fsolve(htpp_yfunc_ylocus_hill48,1.0,args=(s,c,s0))
					xYld.append(xx*alpha)
					yYld.append(yy*alpha)

#	---------------------------------------------------------------- yld2004-18p
		elif yldname == 'yld2004-18p':
			cp1,cp2,a,s0 = htpp_yfunc_yld2004_param(odb_name)
			for i in range(0,len(angle)):
					s,xx,yy = htpp_yfunc_yieldlocus_aux_ori(angle[i],s0,s12)
					alpha = fsolve(htpp_yfunc_ylocus_yld2004,1.0,args=(s,cp1,cp2,a,s0))
					xYld.append(xx*alpha)
					yYld.append(yy*alpha)

#	----------------------------------------------------------------- yld2000-2d
		elif yldname == 'yld2000-2d':
			a,em,s0 = htpp_yfunc_yld2000_2d_param(odb_name)
			for i in range(0,len(angle)):
				s,xx,yy = htpp_yfunc_yieldlocus_aux_ori(angle[i],s0,s12)
				alpha = fsolve(htpp_yfunc_ylocus_yld2000_2d,1.0,args=(s,a,em,s0))
				xYld.append(xx*alpha)
				yYld.append(yy*alpha)
		else:
			return

		xYld_n = [x/s0 for x in xYld]
		yYld_n = [y/s0 for y in yYld]
		strings = ['yieldlocus',yldname,str(s12)]
		fname = '_'.join(strings)
		with open(dir + fname + '.dat','w') as f:
			f.write('Sxx , Syy, Sxx/So , Syy/So \n')
			for i in range(len(xYld)):
				f.write('%f , %f , %f , %f\n' %(xYld[i],yYld[i],xYld_n[i], yYld_n[i]))
# ------------------------------------------------------------------------------
