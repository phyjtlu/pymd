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
T = 300
nrep = 1
# time = 0.658fs #time unit
dt = 0.5
# number of md steps
nmd = 100
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
print q
print(lmp.force(q))

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
#print ("axyz:", axyz)

# we fix the 1st atom
# constraint is a list of vectors.
# the force on along each vector is zerofied in md run
constraint = []
for i in range(3*2):
    tmp = N.zeros(len(lmp.xyz))
    tmp[i] = 1.0
    constraint.append(tmp)
#print ("constraint:", constraint)

# if slist is not given, md will initialize it using xyz
mdrun = md(dt, nmd, T, syslist=None, axyz=axyz, nrep=nrep, npie=1)
# attache lammps driver to md
mdrun.AddLMPint(lmp)
# --------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------
# debye bath
# number of dynamical atoms
gamma = 10**-5
nd = len(lmp.xyz)
eta = gamma*N.identity(nd, N.float)
print("eta:", eta)
# --------------------------------------------------------------------------------------

# -----------------------------------------------------------------------
# atom indices that are connecting to debyge bath
ecats = range(nd/3)
ebl = ebath(ecats, T, mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=eta)
mdrun.AddBath(ebl)
#
ebr = ebath(ecats, T, mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=eta)
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
