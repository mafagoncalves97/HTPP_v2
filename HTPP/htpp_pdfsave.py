from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import os
import sys
import math
import cv2

# ------------------------------------------------------------------------------
# --------------------------------------------------------------------- ADD PATH
sys.path.insert(0, os.getcwd() + '/htpp_plot')
sys.path.insert(0, os.getcwd() + '/htpp_output')
sys.path.insert(0, os.getcwd() + '/htpp_graphic')



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------- MAIN
# ------------------------------------------------------------------------------
def htpp_pdfsave():
	# ------------------------------------------------------ FIND SETTINGS FILES
	path = os.getcwd()[:-4]
	

	size = [[[-0.1,-0.2],[6.0,4.6]],[[0.05,-0.1],[7.725,4.875]],[[0.05,-0.1],[7.725,4.875]]]

	# ------------------------------------------------ CALL DATA SUBROUTINES
	
	odb_name='Job_topopt45'
	strain = cv2.imread(r'htpp_output/Figures/strain_diagram_'+odb_name+'.jpg',cv2.COLOR_BGR2RGB)[...,::-1]
	stress = cv2.imread(r'htpp_output/Figures/stress_diagram_'+odb_name+'.jpg',cv2.COLOR_BGR2RGB)[...,::-1]
	peeq = cv2.imread(r'htpp_output/Figures/colormap_peeq_'+odb_name+'.jpg',cv2.COLOR_BGR2RGB)[...,::-1]

	fig = plt.figure(figsize=(24, 8))
	
	# setting values to rows and column variables
	rows = 1
	columns =3
	
	# reading images
	Image1 = strain
	Image2 = stress
	Image3 = peeq
	
	# Adds a subplot at the 1st position
	fig.add_subplot(rows, columns, 1)
	
	# showing image
	plt.imshow(Image1, aspect='auto')
	plt.axis('off')
	
	# Adds a subplot at the 2nd position
	fig.add_subplot(rows, columns, 2)
	
	# showing image
	plt.imshow(Image2, aspect='auto')
	plt.axis('off')
	
	# Adds a subplot at the 3rd position
	fig.add_subplot(rows, columns, 3)
	
	# showing image
	plt.imshow(Image3[0:4000,600:3300,:], aspect='auto')
	plt.axis('off')

	
	plt.tight_layout(h_pad=1)
	

	plt.savefig(r'htpp_output/Figures/strain__stress_peeq_diagram_'+odb_name+'.jpg')
	plt.show()
	plt.close()




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':
	htpp_pdfsave()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
