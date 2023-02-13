# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2019 replay file
# Internal Version: 2018_09_24-19.41.51 157541
# Run by Utilizador on Fri Feb 10 14:47:00 2023
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(1.11979, 1.1169), width=164.833, 
    height=110.796)
session.viewports['Viewport: 1'].makeCurrent()
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
execfile('htpp_data.py', __main__.__dict__)
#: Model: C:/Users/Utilizador/OneDrive - Universidade de Aveiro/Doutoramento/3ºano/Specimens_comparation/HTPP/Job_topopt45.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       5
#: Number of Node Sets:          5
#: Number of Steps:              1
#: > HTPP DATA  [                    ]    0%  (Job_topopt45)> HTPP DATA  [--                  ] ~ 10%  (Job_topopt45)> HTPP DATA  [----                ] ~ 20%  (Job_topopt45)> HTPP DATA  [------              ] ~ 30%  (Job_topopt45)> HTPP DATA  [--------            ] ~ 40%  (Job_topopt45)> HTPP DATA  [----------          ] ~ 50%  (Job_topopt45)> HTPP DATA  [------------        ] ~ 60%  (Job_topopt45)> HTPP DATA  [--------------      ] ~ 70%  (Job_topopt45)> HTPP DATA  [----------------    ] ~ 80%  (Job_topopt45)> HTPP DATA  [------------------  ] ~ 90%  (Job_topopt45)> HTPP DATA  [--------------------]  100%  (Job_topopt45)> HTPP DATA  [--------------------]  100%  (Job_topopt45)
print 'RT script done'
#: RT script done
