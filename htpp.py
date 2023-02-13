# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
import os
import sys

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------- MAIN
# ------------------------------------------------------------------------------
def main():

	path = os.getcwd()
	os.chdir(path + '/HTPP')
	#os.system('python htpp_run.py')
	#os.system('abaqus cae noGUI=htpp_data.py')
	#os.system('python htpp_compute.py')
	#os.system('python htpp_yfunc.py')
	os.system('python htpp_plot.py')
	#os.system('python htpp_pdfsave.py')

if __name__ == '__main__':
	main()
# -------------------------------------------------------------------------- END
