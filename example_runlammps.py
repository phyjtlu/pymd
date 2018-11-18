#!/home/jtlu/anaconda2/bin/python
#import matplotlib.pyplot as PP
import sys, string, os, time, glob
import numpy as N
import numpy.linalg as LA

from md import *
from phbath import *
from ebath import *
from lammps import *
from matrix import *
from myio import *


#-------------------------------------------------------------------------------------
#temperature
T = 4.2
nrep = 1
#time = 0.658fs #time unit
dt = 0.15
#number of md steps
nmd = 2**10
#transiesta run dir,
#where to find default settings, and relaxed structure *.XV
#SDir="../CGrun/"
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#initialise lammps run
lmp = lammps(infile='in.test')
print lmp.els
#forces...
#q is 1d array made from 
#the displacement from equilibrium in unit of 0.06466 Ang., 
#which is the internal unit of md
q=N.zeros(len(lmp.xyz))
print q
lmp.force(q)
#-------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
print "initialise md"

#Initialize axyz:
#with the format:
# [["C",0,0,0],["C",0,0,1.0]]
#
axyz = []
for i,a in enumerate(lmp.els):
    axyz.append([get_atomname(a),lmp.xyz[i*3],lmp.xyz[i*3+1],lmp.xyz[i*3+2]])
print ("axyz:",axyz)



#if slist is not given, md will initialize it using xyz
mdrun = md(dt,nmd,T,syslist=None,axyz=axyz,nrep=nrep,npie=1)
#attache lammps driver to md
mdrun.AddLMPint(lmp)
#--------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------
#debye bath
#number of dynamical atoms
gamma = 10**-5
nd=len(lmp.xyz)
eta=gamma*N.identity(nd,N.float)
print("eta:",eta)
#--------------------------------------------------------------------------------------

#-----------------------------------------------------------------------
#atom indices that are connecting to debyge bath
ecats=range(nd/3)
eb = ebath(ecats,T,mdrun.dt,mdrun.nmd,wmax=1.,nw=500,bias=0.0,efric=eta)
mdrun.AddBath(eb)
#----------------------------------------------------------------------

#---------------------------------------------------------------------------------
#MD
#------------------------------------------------------------------------------
mdrun.Run()
#------------------------------------------------------------------------------
#close lammps instant
lmp.quit()
#------------------------------------------------------------------------------
