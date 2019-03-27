from md import *
from phbath import *
from ebath import *
from lammpsdriver import *
from matrix import *
from myio import *

lammpsinfile=[
"units         metal",
"atom_style    full",
"dimension     3",
"boundary      p  p  p ",
"read_data  test.data",
"pair_style   rebo",
"pair_coeff * * /opt/lammps/potentials/CH.airebo C H ",
]
#-------------------------------------------------------------------------------------
#temperature
T = 300
nrep = 2
#time = 0.658fs #time unit
dt = 0.25/0.658
#number of md steps
nmd = 10**2
#transiesta run dir,
#where to find default settings, and relaxed structure *.XV
#SDir="../CGrun/"
#-------------------------------------------------------------------------------------

#initialise lammps run
lmp = lammpsdriver(infile=lammpsinfile)
print "initialise md"

constraint = []

fixatoms = range(0*3, (71+1)*3)
fixatoms.extend(range(697*3, (768+1)*3))

for i in fixatoms:
    tmp = N.zeros(len(lmp.xyz))
    tmp[i] = 1.0
    constraint.append(tmp)

# Molecular Junction atom indices
slist = range(342, 426+1)
# atom indices that are connecting to debyge bath
ecatsl = range(72, 341+1)
ecatsr = range(427, 696+1)

dynamicatoms = slist+ecatsl+ecatsr
dynamicatoms.sort()
print "the following atoms are dynamic:\n"
print dynamicatoms
print len(dynamicatoms)

# if slist is not given, md will initialize it using xyz
mdrun = md(dt, nmd, T, syslist=None, axyz=lmp.axyz,writepq=False,
           nrep=nrep, npie=1, constr=constraint,nstep=10**5)
# attache lammps driver to md
mdrun.AddLMPint(lmp)

gamma = 0.0001*2
#gamma = 1.316423628402082*10**-2
ndl = len(ecatsl)
ndr = len(ecatsr)
etal = gamma*N.identity(3*ndl, N.float)
etar = gamma*N.identity(3*ndr, N.float)

ebl = ebath(ecatsl, T*1.1, mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etal,zpmotion=False)
mdrun.AddBath(ebl)

ebr = ebath(ecatsr, T*0.9, mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etar,zpmotion=False)
mdrun.AddBath(ebr)

mdrun.Run()

lmp.quit()