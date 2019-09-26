from md import *
from phbath import *
from ebath import *
from lammpsdriver import *
from matrix import *
from myio import *
from postprocessing import *

lammpsinfile=[
"units metal ",
"dimension 3 ",
"boundary p p p",
"atom_style full",
"read_data test.dat ",
"pair_style rebo ",
"pair_coeff * * CH.airebo C H",
]
#-------------------------------------------------------------------------------------
#temperature
T = 300
delta=0.1
nrep = 2
#time = 0.658fs #time unit
dt = 0.25/0.658
#number of md steps
nmd = 2**10
#transiesta run dir,
#where to find default settings, and relaxed structure *.XV
#SDir="../CGrun/"
#-------------------------------------------------------------------------------------

#initialise lammps run
lmp = lammpsdriver(infile=lammpsinfile)
#
lmp.f0 = lmp.f0*0.0
#print lmp.els
#forces...
#q is 1d array made from 
#the displacement from equilibrium in unit of lmp.conv * 0.06466 Ang., 
#which is the internal unit of md
q=N.zeros(len(lmp.xyz))
#print q
lmp.force(q)

print "initialise md"

constraint = []

fixatoms = range(0*3, (71+1)*3)
fixatoms.extend(range(697*3, (768+1)*3))

for i in fixatoms:
    tmp = N.zeros(len(lmp.xyz))
    tmp[i] = 1.0
    constraint.append(tmp)
print ("constraint:",constraint)
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
mdrun = md(dt, nmd, T, syslist=None, axyz=lmp.axyz,nrep=nrep, npie=1,constr=constraint)
# attache lammps driver to md
mdrun.AddLMPint(lmp)
#--------------------------------------------------------------------------------------

gamma = 10**-4
ndl = len(ecatsl)
ndr = len(ecatsr)
etal = gamma*N.identity(3*ndl, N.float)
etar = gamma*N.identity(3*ndr, N.float)
#--------------------------------------------------------------------------------------

#-----------------------------------------------------------------------
#atom indices that are connecting to bath
ebl = ebath(ecatsl, T*(1+delta), mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etal,zpmotion=False)
mdrun.AddBath(ebl)

ebr = ebath(ecatsr, T*(1-delta), mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etar,zpmotion=False)
mdrun.AddBath(ebr)
#----------------------------------------------------------------------

#---------------------------------------------------------------------------------
#MD
#------------------------------------------------------------------------------
mdrun.Run()
#------------------------------------------------------------------------------
#close lammps instant
lmp.quit()
CTC(delta=delta)
#----------------
