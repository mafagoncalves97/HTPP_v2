import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# use LaTeX fonts in the plot
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')

########################## IDENTIFIABILITY INDEX ###########################
Best_solutions = list(reversed(['Notched', 'Sigma', 'D', 'Shape','TopOpt']))
x_axis = np.arange(len(Best_solutions))
S_r_plastic = list(reversed([3.19,4.37,4.56,3.80,0]))
S_r_total = list(reversed([1.84,0.45,0.99,0.63,0]))

fig, ax = plt.subplots(figsize =(6, 3))

plt.barh(x_axis+0.2,S_r_total, height=0.4, label = 'Swift law',color=['dimgrey'])
plt.barh(x_axis-0.2,S_r_plastic, height=0.4, label = 'Yld2000-2d', color=['lightgrey'])
plt.xlabel(r'Identifiability index $I_\mathrm{k}$')
plt.yticks(x_axis, Best_solutions)
#plt.gca().xaxis.grid(True)
plt.axvline(x=2, color='grey',linestyle='--')
plt.axvline(x=3, color='grey',linestyle='--')
plt.xlim([0, 5])
plt.legend(ncols=2,loc=1,prop={'size': 8})

# Add annotation to bars
for i in ax.patches:
    plt.text(i.get_width()+0.1, i.get_y()+0.1,
             str(round((i.get_width()), 2)),
             fontsize = 8,
             color ='black')

plt.tight_layout()

#plt.show()
plt.close()

########################## COLLINEARITY INDEX ###########################
Best_solutions = list(reversed(['Notched', 'Sigma', 'D', 'Shape','TopOpt']))
x_axis = np.arange(len(Best_solutions))
S_r_plastic = list(reversed([16.98,76.25,94.83,41.03,0]))
S_r_total = list(reversed([4.99,1.3,2.32,1.54,0]))

fig, ax = plt.subplots(figsize =(7, 3))

plt.barh(x_axis+0.2,S_r_total, height=0.4, label = 'Swift law',color=['dimgrey'])
plt.barh(x_axis-0.2,S_r_plastic, height=0.4, label = 'Yld2000-2d', color=['lightgrey'])
plt.xlabel(r'Collinearity index $\gamma_\mathrm{k}$')
plt.yticks(x_axis, Best_solutions)
plt.axvline(x=20, color='grey',linestyle='--')
#plt.gca().xaxis.grid(True)
plt.xlim([0, 100])
plt.legend(ncols=2,loc=1,prop={'size': 8})

# Add annotation to bars
for i in ax.patches:
    plt.text(i.get_width()+0.2, i.get_y()+0.1,
             str(round((i.get_width()), 2)),
             fontsize = 8,
             color ='black')

plt.tight_layout()

#plt.show()
plt.close()

########################## DETERMINANT MEASURE ###########################
Best_solutions = list(reversed(['Notched', 'Sigma', 'D', 'Shape','TopOpt']))
x_axis = np.arange(len(Best_solutions))
S_r_plastic = list(reversed([44.04,125.97,125.14,61.68,0]))
S_r_total = list(reversed([32.23,49.37,13.96,46.47,0]))

fig, ax = plt.subplots(figsize =(6, 3))

plt.barh(x_axis+0.2,S_r_total, height=0.4, label = 'Swift law',color=['dimgrey'])
plt.barh(x_axis-0.2,S_r_plastic, height=0.4, label = 'Yld2000-2d', color=['lightgrey'])
plt.xlabel(r'Determinant measure $\rho_\mathrm{k}$')
plt.yticks(x_axis, Best_solutions)
#plt.gca().xaxis.grid(True)
plt.xlim([0, 150])
plt.legend(ncols=2,loc=1,prop={'size': 8})

# Add annotation to bars
for i in ax.patches:
    plt.text(i.get_width()+0.2, i.get_y()+0.1,
             str(round((i.get_width()), 2)),
             fontsize = 8,
             color ='black')

plt.tight_layout()

#plt.show()
plt.close()


########################## SENSITIVITY MEASURE ###########################

Notched=list([11.09,10.55,21.23,36.33,27.69,17.86,58.12,21.69,6.91,11.60,107.8])
Kim=list([705.38,475.16,14.75,41.34,468.74,12.78,459.43,15.94,2.07,2.21,1140.82])
Jones=list([269.28,16.11,268.21,47.53,44.44,15.59,300.72,90.92,2.28,2.95,25.96])
Conde=list([13.68,12.45,21.15,57.19,51.35,18.63,457.78,31.20,2.06,3.35,808.77])
TopOpt=list([1,1,1,1,1,1,1,1,1,1,1])


fig, ax = plt.subplots(figsize =(8,4))
c = ['lightgrey',  'darkgrey', 'dimgrey','black']

Best_solutions = list([r'$\alpha_1$',r'$\alpha_2$',r'$\alpha_3$',r'$\alpha_4$',r'$\alpha_5$',r'$\alpha_6$',r'$\alpha_7$',r'$\alpha_8$',r'$K$',r'$\varepsilon_\mathrm{0}$',r'$n$'])
x_axis = 3*np.arange(len(Best_solutions))


plt.bar(x_axis-1,Notched, width=0.5, label = 'Notched',color=['lightgrey'])
plt.bar(x_axis-0.5,Kim, width=0.5, label = 'Kim', color=['darkgrey'])
plt.bar(x_axis,Jones, width=0.5, label = 'D',color=['dimgrey'])
plt.bar(x_axis+0.5,Conde, width=0.5, label = 'Shape', color=['black'])

plt.bar(x_axis+1,TopOpt, width=0.5, label = 'TopOpt',color=['dimgrey'])

plt.ylabel(r'Sensitivity measure $\delta_\mathrm{k}$')

plt.xticks(x_axis, Best_solutions)
#plt.gca().xaxis.grid(True)
plt.ylim([0, 1200])
plt.legend(ncols=5,loc=1,prop={'size': 8})

# # Add annotation to bars
# for i in ax.patches:
#     plt.text(i.get_width()+0.2, i.get_y()+0.1,
#              str(round((i.get_width()), 2)),
#              fontsize = 8,
#              color ='black')

plt.tight_layout()

plt.show()
plt.close()



'''Best_solutions = list(reversed(['Notched', 'Sigma', 'D', 'Shape','TopOpt']))
x_axis = np.arange(len(Best_solutions))
S_r_plastic = list(reversed([23.3, 44.9, 45.3, 35.7, 64.6]))
S_r_total = list(reversed([34.7,49.3,55.4,51.5,64.6]))

fig, ax = plt.subplots(figsize =(5, 3))

plt.barh(x_axis+0.2,S_r_total, height=0.4, label = 'Total',color=['dimgrey'])
plt.barh(x_axis-0.2,S_r_plastic, height=0.4, label = 'Plastic', color=['lightgrey'])
plt.xlabel(r'Stress triaxiality and Lode angle indicator ($S_\mathrm{LT}$)')
plt.yticks(x_axis, Best_solutions)
#plt.gca().xaxis.grid(True)
plt.xlim([0, 100])
plt.legend(ncols=2,loc=1,prop={'size': 8})

# Add annotation to bars
for i in ax.patches:
    plt.text(i.get_width()+0.2, i.get_y()+0.1,
             str(round((i.get_width()), 2)),
             fontsize = 8,
             color ='black')

plt.tight_layout()

plt.show()'''