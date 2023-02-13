# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
#from odbAccess import *
#from abaqusConstants import *
#import numpy as np
import os
#from textRepr import prettyPrint



def htpp_run_simulation_option(odb,option):

	input_file = '\n'.join(s for s in option if 'name=' in s).replace(",", "").replace("\n", "")[6:-1]

	return input_file



# ------------------------------------------------------------------------------
# ----------------------------------------------------------------------- RUN
# ------------------------------------------------------------------------------
def htpp_run_simulation(option,odb_name):


	input_file = htpp_run_simulation_option(odb_name,option)
	#When asked for a Input file, write ..\inputfilename
	os.system('abaqus job='+input_file+' user=..\\\\UMMDp_FLC.f interactive')
	
	


# ------------------------------------------------------------------------------
