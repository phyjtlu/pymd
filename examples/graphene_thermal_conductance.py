#!/usr/bin/env python
import time

from ebath import *
from lammpsdriver import *
from matrix import *
from md import *
from myio import *
from phbath import *
from postprocessing import *

lammpsinfile = [
    #"log none",
    "units metal ",
    "dimension 3 ",
    "boundary f p p",
    "atom_style full",
    "read_data graphene.data ",
    "pair_style rebo ",
    "pair_coeff * * CH.rebo C",
    "min_style cg",
    "minimize 1e-25 1e-25 5000 10000",
]
# temperature
T = 300
delta = 0.2
Thot = T*(1+delta/2)
Tcold = T*(1-delta/2)
nstart = 0
nstop = 5
# time = 0.658fs #time unit
dt = 0.25/0.658
# number of md steps
nmd = 2**10
# initialise lammps run
lmp = lammpsdriver(infile=lammpsinfile)
time_start = time.time()

print("initialise md")
fixatoms = list(range(0*3, (7+1)*3))
fixatoms.extend(list(range(88*3, (95+1)*3)))

# print(("constraint:",constraint))
# Molecular Junction atom indices
slist = list(range(24*3, (71+1)*3))
cutslist = [list(range(24*3, (39+1)*3)),
            list(range(40*3, (55+1)*3)), list(range(56*3, (71+1)*3))]
# atom indices that are connecting to debyge bath
ecatsl = list(range(8*3, (23+1)*3))
ecatsr = list(range(72*3, (87+1)*3))
dynamicatoms = slist+ecatsl+ecatsr
dynamicatoms.sort()
print("the following atoms are dynamic:\n")
print(dynamicatoms)
print(len(dynamicatoms))
# if slist is not given, md will initialize it using xyz
mdrun = md(dt, nmd, T, ecatsl=ecatsl, ecatsr=ecatsr, slist=slist, cutslist=None, syslist=None, axyz=lmp.axyz, writepq=True, rmnc=True,
           nstart=nstart, nstop=nstop, npie=1, constr=fixatoms, nstep=100)
# attache lammps driver to md
mdrun.AddLMPint(lmp)
# unit in 0.658211814201041 fs
damp = 25/0.658211814201041

etal = (1.0/damp)*N.identity(len(ecatsl), N.float)
etar = (1.0/damp)*N.identity(len(ecatsr), N.float)
# atom indices that are connecting to bath
ebl = ebath(ecatsl, Thot, mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etal, classical=False, zpmotion=False)
mdrun.AddBath(ebl)

ebr = ebath(ecatsr, Tcold, mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etar, classical=False, zpmotion=False)
mdrun.AddBath(ebr)
# MD
mdrun.Run()
# close lammps instant
lmp.quit()
CalHF()
CalTC(delta=delta, dlist=0)
time_end = time.time()
print('time cost', time_end-time_start, 's')

# ----------------
