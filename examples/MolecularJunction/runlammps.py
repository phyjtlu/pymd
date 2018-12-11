import sys
import string
import os
import time
import glob
import numpy as N
import numpy.linalg as LA

from md import *
from phbath import *
from ebath import *
from lammps import *
from matrix import *
from myio import *

# -------------------------------------------------------------------------------------
# temperature
T = 1000
nrep = 1
# time = 0.658fs #time unit
dt = 0.25/0.658
# number of md steps
nmd = 10**2
# transiesta run dir,
# where to find default settings, and relaxed structure *.XV
# SDir="../CGrun/"
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# initialise lammps run
args = "-screen none"
lmp = lammps(infile='in.langevin', cmdargs=args.split())
#
#lmp.f0 = lmp.f0*0.0
#print lmp.els
# forces...
# q is 1d array made from
# the displacement from equilibrium in unit of lmp.conv * 0.06466 Ang.,
# which is the internal unit of md
q = N.zeros(len(lmp.xyz))
print (q)
print (lmp.f0)
print (lmp.force(q))

# -------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------
print "initialise md"

# Initialize axyz:
# with the format:
# [["C",0,0,0],["C",0,0,1.0]]
#
axyz = []
for i, a in enumerate(lmp.els):
    axyz.append([get_atomname(a), lmp.xyz[i*3],
                 lmp.xyz[i*3+1], lmp.xyz[i*3+2]])
print ("axyz:", axyz)

fileObject = open('axyz.txt', 'w+')  
for ip in axyz:  
    fileObject.write(str(ip))
    fileObject.write('\n')  
fileObject.close() 

# we fix the 1st atom
# constraint is a list of vectors.
# the force on along each vector is zerofied in md run
constraint = []

fixatoms = range(0*3, (71+1)*3)
fixatoms.extend(range(409*3, (480+1)*3))

for i in fixatoms:
    tmp = N.zeros(len(lmp.xyz))
    tmp[i] = 1.0
    constraint.append(tmp)

print ("constraint:", constraint)

# Molecular Junction atom indices
slist = range(198, 282+1)
# -----------------------------------------------------------------------
# atom indices that are connecting to debyge bath
ecatsl = range(72, 197+1)
ecatsr = range(283, 408+1)

dynamicatoms = slist+ecatsl+ecatsr
dynamicatoms.sort()
print "the following atoms are dynamic:\n"
print dynamicatoms
print len(dynamicatoms)

# if slist is not given, md will initialize it using xyz
mdrun = md(dt, nmd, T, syslist=None, axyz=axyz,
           nrep=nrep, npie=1, constr=constraint)
# attache lammps driver to md
mdrun.AddLMPint(lmp)
# --------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------
# debye bath
# number of dynamical atoms
gamma = 1.316423628402082*10**-2
ndl = len(ecatsl)
ndr = len(ecatsr)
etal = gamma*N.identity(3*ndl, N.float)
etar = gamma*N.identity(3*ndr, N.float)
#print("eta:", eta)
# --------------------------------------------------------------------------------------
ebl = ebath(ecatsl, T, mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etal)
mdrun.AddBath(ebl)

ebr = ebath(ecatsr, T, mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etar)
mdrun.AddBath(ebr)
# ----------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# MD
# ------------------------------------------------------------------------------
mdrun.Run()
# ------------------------------------------------------------------------------
# close lammps instant
lmp.quit()
# ------------------------------------------------------------------------------
