#!/usr/bin/env python
"""
A script which averages a CHGCAR or LOCPOT file in one direction to make a 1D curve.
User must specify filename and direction on command line。
Example: vtotav.py LOCPOT z
Depends on ase
"""

import os
import sys
import numpy as np
import math
import string
import datetime
import time
import re
from ase.calculators.vasp import VaspChargeDensity

starttime = time.perf_counter()
print("Starting calculation at", end=' ')
print(time.strftime("%H:%M:%S on %a %d %b %Y"))

if len(sys.argv) != 3:
    print("\n** ERROR: Must specify name of file and direction on command line.")
    print("eg. vtotav.py z 1e-5.")
    sys.exit(0)

# Read information from command line
# First specify location of LOCPOT
#LOCPOTfile = sys.argv[1].lstrip()
LOCPOTfile = "LOCPOT"
CHGCARfile = "CHGCAR"
OUTCARfile = "OUTCAR"

# Next the direction to make average in
# input should be x y z, or X Y Z. Default is Z.
allowed = "xyzXYZ"
direction = sys.argv[1].lstrip()
cutoff = float(sys.argv[2].lstrip())
if allowed.find(direction) == -1 or len(direction)!=1 :
    print("** WARNING: The direction was input incorrectly.")
    print("** Setting to z-direction by default.")
if direction.islower():
    direction = direction.upper()
filesuffix = "_%s" % direction

# Open geometry and density class objects
#-----------------------------------------
vasp_charge = VaspChargeDensity(filename = CHGCARfile)
vasp_potl = VaspChargeDensity(filename = LOCPOTfile)
e_femi = float(re.findall(r"E-fermi\s*:\s*([\d.-]+)", open(OUTCARfile).read())[0])
charge = vasp_charge.chg[-1]
potl = vasp_potl.chg[-1]
atoms = vasp_charge.atoms[-1]
del vasp_charge
del vasp_potl

potl=potl*atoms.get_volume()

print("\nReading file: %s, %s, %s" % (LOCPOTfile, CHGCARfile, OUTCARfile))
print("Performing average in %s direction" % direction)
print("Cutoff density is %10.6f" % cutoff)

# Read in lattice parameters and scale factor
#---------------------------------------------
cell = atoms.cell

# Find length of lattice vectors
#--------------------------------
latticelength = np.dot(cell, cell.T).diagonal()
latticelength = latticelength**0.5

print("Lattice parameters:")
print("a = %10.6f Ang" % latticelength[0])
print("b = %10.6f Ang" % latticelength[1])
print("c = %10.6f Ang" % latticelength[2])

# Read in potential data
#------------------------
ngridpts = np.array(potl.shape)
totgridpts = ngridpts.prod()
print("Potential stored on a %dx%dx%d grid" % (ngridpts[0],ngridpts[1],ngridpts[2]))
print("Total number of points is %d" % totgridpts)
print("Reading potential data from file...", end=' ')
sys.stdout.flush()
print("done.")

# Perform average
#-----------------
if direction=="X":
    idir = 0
    a = 1
    b = 2
elif direction=="Y":
    a = 0
    idir = 1
    b = 2
else:
    a = 0
    b = 1
    idir = 2
a = (idir+1)%3
b = (idir+2)%3
# At each point, sum over other two indices
average_potl = np.zeros(ngridpts[idir],np.float64)
for ipt in range(ngridpts[idir]):
    if direction=="X":
        average_potl[ipt] = potl[ipt,:,:].sum()
    elif direction=="Y":
        average_potl[ipt] = potl[:,ipt,:].sum()
    else:
        average_potl[ipt] = potl[:,:,ipt].sum()

average_charge = np.zeros(ngridpts[idir],np.float64)
for ipt in range(ngridpts[idir]):
    if direction=="X":
        average_charge[ipt] = charge[ipt,:,:].sum()
    elif direction=="Y":
        average_charge[ipt] = charge[:,ipt,:].sum()
    else:
        average_charge[ipt] = charge[:,:,ipt].sum()

# Scale by number of grid points in the plane.
# The resulting unit will be eV.
average_potl /= ngridpts[a]*ngridpts[b]

# Scale by size of area element in the plane,
# gives unit e/Ang. I.e. integrating the resulting
# CHG_dir file should give the total charge.
area = np.linalg.det([ (cell[a,a], cell[a,b] ),
                        (cell[b,a], cell[b,b])])
dA = area/(ngridpts[a]*ngridpts[b])
average_charge *= dA

for i,charge in enumerate(average_charge):
    if charge < cutoff:
        vac_pos = i #*latticelength[idir]/float(ngridpts[idir]-1)
        print("Vacuum distance is %10.6f Ang" % (vac_pos*latticelength[idir]/float(ngridpts[idir]-1)))
        break

e_vac = average_potl[vac_pos]

print("E_vac is %10.6f eV" % e_vac)
print("E_fermi is %10.6f eV" % e_femi)
print("Work function  = E_vac - E_fermi is %10.6f eV" % (e_vac-e_femi))
print("HER potential = E_vac - E_fermi - 4.44 eV is %10.6f eV" % (e_vac-e_femi-4.44))

# Print out average
#-------------------

averagepotl = LOCPOTfile + filesuffix
averagecharge = CHGCARfile + filesuffix
print("Writing averaged data to file %s and %s..." % (averagepotl,averagecharge), end=' ')
sys.stdout.flush()

outputfile = open(averagepotl,"w")
outputfile.write("#  Distance(Ang)     Potential(eV)\n")

outputfile1 = open(averagecharge,"w")
outputfile1.write("#  Distance(Ang)     Chg. density (e/Ang)\n")

xdiff = latticelength[idir]/float(ngridpts[idir]-1)
for i in range(ngridpts[idir]):
    x = i*xdiff
    outputfile.write("%15.8g %15.8g\n" % (x,average_potl[i]))
outputfile.close()

for i in range(ngridpts[idir]):
    x = i*xdiff
    outputfile1.write("%15.8g %15.8g\n" % (x,average_charge[i]))
outputfile1.close()

print("done.")

endtime = time.perf_counter()
runtime = endtime-starttime
print("\nEnd of calculation.")
print("Program was running for %.2f seconds." % runtime)