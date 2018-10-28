#!/applications/mbrsoft/bin/python
#import matplotlib.pyplot as PP
import sys, string, os, time, glob
import Scientific.IO.NetCDF as nc
import numpy as N
import numpy.linalg as LA

from Inelastica import SiestaIO as SIO
from md import *
from phbath import *
from ebath import *
from siesta import *
from matrix import *
from io import *


#-------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
T = 4.2
nrep = 16
#time = 0.658fs #time unit
dt = 0.15
nmd = 2**16
#transiesta run dir,
#where to find default settings, and relaxed structure *.XV
SDir="../CGrun/"

nperlayer=1    #number of atoms per layer
# We need to ensure correct hamiltonian, so no cut
cut = 0        #how many layers to cut from each side

mesh=200.
dmtol=0.0001
#-------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
import Inelastica.SiestaIO as SIO

os.system("cp "+SDir+"/*.vps .")
print "read in structure information"
fn=glob.glob(SDir+"/*.XV")
print fn

# Cut layers
geom=cutlayers(fn[0],nalayer=nperlayer,nl=cut,nr=cut,outfile="STRUCT.fdf")
cell, xyzonly, snr, anr = geom.pbc,geom.xyz,geom.snr,geom.anr
print anr

#xyz with element labels
#for example ['Au',0,0,0]
xyz = [[PeriodicTable[a],b[0],b[1],b[2]] for (a,b) in zip(anr,xyzonly)]
print xyz


#-------------------------------------------------------------------------------------
#initialise siesta run
#def __init__(self,label,xyz,cell,mesh=100.,dmtol=0.001,constraints=[],lunit="Ang",eunit="eV"):
srun = siesta("mdrun",xyz,cell,mesh=mesh,dmtol=dmtol,tdir=SDir)
print srun.els

srun.start()

#forces...
#q is 1d array made from 
#the displacement from equilibrium in unit of 0.06466 Ang., 
#which is the internal unit of md
q=N.zeros(len(xyz))
print q
srun.force(q)

srun.quit()
