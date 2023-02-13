# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
import os
import shutil
import sys
import math
import numpy as np

# ------------------------------------------------------------------------------
# --------------------------------------------------------------------- ADD PATH
sys.path.insert(0, os.getcwd() + '/htpp_compute')

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------- IMPORT FILES
from htpp_compute_indicator import htpp_compute_indicator
from htpp_compute_BBindicator import htpp_compute_BBindicator
from htpp_compute_DICFEAinterpolation import htpp_compute_DICFEAinterpolation

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------- MAIN
# ------------------------------------------------------------------------------
def htpp_compute():
	# ------------------------------------------------------ FIND SETTINGS FILES
	path = os.getcwd()[:-4]
	files = []
	for file in os.listdir(path):
		if file.endswith('.htpp'):
			files.append(file)

	for file in files:
		# ------------------------------------------ OPEN AND READ SETTINGS FILE
		data_settings = []
		try:
			path = os.getcwd()[:-4]
			with open(path + file, 'r') as f:
				for line in f:
					if len(line.split()) > 0:
						if line.split()[0][0:2][0:2] == '\\*':
							continue
						else:
							data_settings.append(line.split(", "))
		except Exception as e:
			print('Error: ' + file)
			return

		i = 1
		for line in data_settings:
			if '*ODB' in line[0].upper():
				odb_file = line[-1][6:-2]
				name = odb_file
			elif '*COMPUTE' in line[0].upper():
				data_begin = i
			elif '*ENDCOMPUTE' in line[0].upper():
				data_end = i
			i+=1
		# ------------------------------------------------- EXTRACT DATA OPTIONS
		data_options = data_settings[data_begin:data_end]
		all_data = []
		for line in data_options:
			if '**INDICATOR' in line[0].upper():
				all_data.append(line)
			elif '**BBINDICATOR' in line[0].upper():
				all_data.append(line)
			elif '**INTERPOLATION' in line[0].upper():
				all_data.append(line)

		# ------------------------------------------------ CALL DATA SUBROUTINES
		if len(all_data) > 0:
			tlen = 20.0
			st = int(round(tlen/len(all_data)))
			stri = "> HTPP COMPUTE "
			sys.stdout.write(stri+"[%s]    %d%s  (%s)"%(" "*int(tlen), 0,'%',name))
			sys.stdout.flush()
			n = 1
			for option in all_data:
				opt = option[0]
				if opt == '**INDICATOR':
					htpp_compute_indicator(option,name)
				elif opt == '**BBINDICATOR':
					htpp_compute_BBindicator(option,name)
				elif opt == '**INTERPOLATION':
					htpp_compute_DICFEAinterpolation(option,name)

				per = 100.0*float(st*n)/float(tlen)
				bars = st*n
				form = "[%s%s] ~ %d%s  (%s)"
				if per >= 100.0:
					per = 100.0
					bars = int(tlen)
					form = "[%s%s]  %d%s  (%s)"
				sys.stdout.write("\r"+stri+form %("-"*bars," "*(int(tlen)-bars),per,'%',name))
				sys.stdout.flush()
				n+=1

			sys.stdout.write("\r"+stri+"[%s]  %d%s  (%s)\n"%("-"*int(tlen), 100,'%',name))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':
	htpp_compute()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
