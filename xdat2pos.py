import ase.io.vasp as vasp
import numpy as np
import sys
import os
import shutil
import re

print("This is a script to convert selected frame XDATCAR to POSCAR.")
print("Select frame is "+' '.join(str(x) for x in sys.argv[1:]))

INCAR_wf = '''System=Pt111
ISTART=1
ICHARG=1
ISMEAR=1
SIGMA=0.2
NWRITE=0
GGA=RP
ENCUT=400
NSW=0
EDIFF=1E-6
LREAL=Auto
NELM=500
NELMIN=8
ALGO=N
EDIFFG=-0.05
POTIM=0.5
ISIF=2
LWAVE=FALSE
LCHARG=T
IBRION=-1
NCORE=10
ISYM=0
LVHAR =.TRUE.
IVDW=11'''

print("Warning: INCAR is set to default. You can check in source code.")

frame = [int(x) for x in sys.argv[1:]]

current_folder = os.getcwd()

def xdat2pos(xdatcar,step):

    traj = vasp.read_vasp_xdatcar(xdatcar,index = step)

    folder = os.path.basename(os.getcwd()) + '_' + str(step)

    os.mkdir(folder)
    
    vasp.write_vasp(current_folder+'/'+folder+'/'+'POSCAR',traj,direct=True, sort=False)

    if os.path.exists(current_folder+'/'+'KPOINTS'):
        shutil.copy(current_folder+'/'+'KPOINTS',current_folder+'/'+folder)
    else:
        print("Warning: KPOINTS not found!")

    if os.path.exists(current_folder+'/'+'POTCAR'):
        shutil.copy(current_folder+'/'+'POTCAR',current_folder+'/'+folder)
    else:
        print("Warning: POTCAR not found!")
    
    if os.path.exists(current_folder+'/'+'sub.sh'):
        shutil.copy(current_folder+'/'+'sub.sh',current_folder+'/'+folder)
    else:
        print("Warning: sub.sh not found!")

    if os.path.exists(current_folder+'/'+'SELECTED_ATOMS_LIST'):
        shutil.copy(current_folder+'/'+'SELECTED_ATOMS_LIST',current_folder+'/'+folder)

    with open(current_folder+'/'+folder+'/'+'INCAR', 'w') as file:
        file.write(INCAR_wf)

for i in frame:

    print("Converting frame "+str(i))

    xdat2pos(current_folder+'/'+'XDATCAR',i)

print("Done")
