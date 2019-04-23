from md import *
from phbath import *
from ebath import *
from lammpsdriver import *
from matrix import *
from myio import *

lammpsinfile=[
"units metal ",
"dimension 3 ",
"boundary p p p",
"atom_style full",
"read_data test.dat ",
"pair_style rebo ",
"pair_coeff * * CH.airebo C",
]
#-------------------------------------------------------------------------------------
#temperature
T = 4.2
nrep = 1
#time = 0.658fs #time unit
dt = 0.5
#number of md steps
nmd = 2**14
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
print q
lmp.force(q)
print "initialise md"
#we fix the 1st atom
#constraint is a list of vectors.
#the force on along each vector is zerofied in md run
constraint = []
for i in range(3*2):
    tmp = N.zeros(len(lmp.xyz))
    tmp[i]=1.0
    constraint.append(tmp)
print ("constraint:",constraint)

# if slist is not given, md will initialize it using xyz
mdrun = md(dt, nmd, T, syslist=None, axyz=lmp.axyz,writepq=True,nrep=nrep, npie=1)
# attache lammps driver to md
mdrun.AddLMPint(lmp)
#--------------------------------------------------------------------------------------

gamma = 10**-5
nd=len(lmp.xyz)
eta=gamma*N.identity(nd,N.float)
print("eta:",eta)
#--------------------------------------------------------------------------------------

#-----------------------------------------------------------------------
#atom indices that are connecting to debyge bath
ecats=range(nd/3)
eb = ebath(ecats,T,mdrun.dt,mdrun.nmd,wmax=1.,nw=500,bias=0.0,efric=eta)
#mdrun.AddBath(eb)
#----------------------------------------------------------------------

#---------------------------------------------------------------------------------
#MD
#------------------------------------------------------------------------------
mdrun.Run()
#------------------------------------------------------------------------------
#close lammps instant
lmp.quit()
#----------------
