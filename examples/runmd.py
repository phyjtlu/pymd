import time
import numpy as np
from baths import ebath
from tools import calHF, calTC
from lammpsdriver import lammpsdriver
from md import md

lammpsinfile = [
    #'log none'
    'atom_style full',
    'units metal',
    'boundary f p p',
    'read_data structure.data',
    'pair_style rebo',
    'pair_coeff * * CH.rebo C H',
]
# temperature
T = 300
delta = 0.1
nstart = 0
nstop = 5
# time = 0.658fs #time unit
dt = 0.25/0.658
# number of md steps
nmd = 2**12
# initialise lammps run
lmp = lammpsdriver(infile=lammpsinfile)
time_start = time.time()

print("initialise md")
fixatoms = [range(0*3, (19+1)*3), range(181*3, (200+1)*3)]

# Molecular Junction atom indices
slist = list(range(70*3, (130+1)*3))
cutslist = [list(range(70*3, (89+1)*3)),
            list(range(90*3, (109+1)*3)), list(range(110*3, (130+1)*3))]
# atom indices that are connecting to debyge bath
ecatsl = list(range(20*3, (69+1)*3))
ecatsr = list(range(131*3, (180+1)*3))
# if slist is not given, md will initialize it using xyz
mdrun = md(dt, nmd, T, syslist=None, axyz=lmp.axyz, writepq=True, rmnc=False,
           nstart=nstart, nstop=nstop, npie=1, constr=fixatoms, nstep=100)
# attache lammps driver to md
mdrun.AddPotential(lmp)
# unit in 0.658211814201041 fs
damp = 100/0.658211814201041

etal = (1.0/damp)*np.identity(len(ecatsl), np.float)
etar = (1.0/damp)*np.identity(len(ecatsr), np.float)
# atom indices that are connecting to bath
ebl = ebath(ecatsl, T*(1+delta/2), mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etal, classical=False, zpmotion=True)
mdrun.AddBath(ebl)
ebr = ebath(ecatsr, T*(1-delta/2), mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etar, classical=False, zpmotion=True)
mdrun.AddBath(ebr)

mdrun.AddPowerSection([ecatsl, slist, ecatsr])
mdrun.Run()
lmp.quit()
calHF()
calTC(delta=delta)
time_end = time.time()
print('time cost', time_end-time_start, 's.')
