# Semi-classical Langevin Molecular Dynamics

This is a set of scripts to run semi-classical Langevin molecular dynamics on junctions connecting to several electron or phonon baths. The details of the method are presented in Prog. Surf. Sci. [<https://doi.org/10.1016/j.progsurf.2018.07.002]>.

To do the molecular dynamics, we need a force driver. Currently, we implemented Siesta, Brenner and LAMMPS potentials.

The important feature of this script is that, the statistics of the bath degrees of freeom is quantum. For example, the zero point fluctuations are included; they fulfill the quantum-mechanical Bose-Einstein distribution. Moreover, the electron bath is allowed to be in a nonequilibrium steady state (non-thermal).

## Install

Install LAMMPS and required python packages, which can be done by following script.

```bash
#!/bin/sh
# install.sh Compile & install LAMMPS shared library & required packages in python
# Usage: sh install.sh
# -------------------------------------------------------
if [ -d 'lammps' ]
then
    echo 'LAMMPS already exists' 
else
    echo 'Download LAMMPS source code'
    git clone https://github.com/lammps/lammps.git -b stable
fi
echo 'Activate Intel compilation environment'
source /opt/intel/parallel_studio_xe_2020/psxevars.sh
cd lammps/src
echo 'Install select packages'
make yes-body yes-class2 yes-manybody yes-molecule yes-kspace yes-user-reaxc yes-user-phonon
echo 'Compile LAMMPS shared library'
make -j 8 mode=shlib intel_cpu
echo 'Install LAMMPS shared library in python'
make install-python
echo 'Compile the LAMMPS executable file'
make -j 8 intel_cpu
echo 'Install required python packages'
python -m pip install -U numpy netCDF4 tqdm
```

## Input file

### QTB-MD thermal conductance calcuation

```python
import time
import numpy as np
from baths import ebath
from tools import calHF, calTC
from lammpsdriver import lammpsdriver
from md import md

lammpsinfile = [
    #'log none',
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
nstop = 3
# time = 0.658fs #time unit
dt = 0.25/0.658
# number of md steps
nmd = 2**12
# initialise lammps run
lmp = lammpsdriver(infile=lammpsinfile)
time_start = time.time()

print('initialise md')
fixatoms = [range(0*3, (19+1)*3), range(181*3, (200+1)*3)]

# Molecular Junction atom indices
slist = range(70*3, (130+1)*3)
cutslist = [range(70*3, (89+1)*3),
            range(90*3, (109+1)*3), range(110*3, (130+1)*3)]
# atom indices that are connecting to debyge bath
ecatsl = range(20*3, (69+1)*3)
ecatsr = range(131*3, (180+1)*3)
# if slist is not given, md will initialize it using xyz
mdrun = md(dt, nmd, T, axyz=lmp.axyz,
           nstart=nstart, nstop=nstop)
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

mdrun.AddConstr(fixatoms)

mdrun.CalPowerSpec()
mdrun.AddPowerSection([ecatsl, slist, ecatsr])
mdrun.CalAveStruct()
mdrun.SaveTraj()
mdrun.SaveAll()
mdrun.RemoveNC()

mdrun.Run()
lmp.quit()
calHF()
calTC(delta=delta)
time_end = time.time()
print('time cost', time_end-time_start, 's.')
```

### NEGF thermal conductance calcuation

```python
import time
import numpy as np
from negf import bpt
from matplotlib import pyplot as plt
lammpsinfile = [
    'atom_style full',
    'units metal',
    'boundary f p p',
    'read_data structure.data',
    'pair_style rebo',
    'pair_coeff * * CH.rebo C H',
]
time_start = time.time()
atomfixed = [range(0*3, (19+1)*3), range(181*3, (200+1)*3)]
atomofbath = [range(20*3, (69+1)*3), range(131*3, (180+1)*3)]
mybpt = bpt(lammpsinfile, 0.25, 0.1, atomofbath, atomfixed, 100)
mybpt.plotresult()
# T_H/C = T*(1±delta/2)
T = [100, 200, 300, 400, 500, 600, 700,
     800, 900, 1000]
delta = 0.1
thermalconductance = []
for temp in T:
    thermalconductance.append([temp, mybpt.thermalconductance(temp, delta)])
np.savetxt('thermalconductance.dat', thermalconductance)
plt.figure(5)
plt.plot(np.array(thermalconductance)[
    :, 0], np.array(thermalconductance)[:, 1])
plt.xlabel('Temperature(K)')
plt.ylabel('Thermal Conductance(nW/K)')
plt.savefig('thermalconductance.png')
time_end = time.time()
print('time cost', time_end-time_start, 's.')
```
