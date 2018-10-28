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
T=0.
#curcof = 243414.
nrep = 8
#time = 0.658fs #time unit
dt = 2
nmd = 2**14
#transiesta run dir,
#where to find default settings, and relaxed structure *.XV
SDir="../"

nperlayer=16    #number of atoms per layer
nl=6            #number of layers left
nr=6            #number of layers right
cut=nl-3        #how many layers to cut from each side


#friction coefficient
gamma=0.0001
mesh=150.
dmtol=0.0001
#-------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
import Inelastica.SiestaIO as SIO
#read in the relaxed structure
#cell   3x3 unit vectors
#snr    species numbers defined in fdf
#anr    atomic numbers 
#xyzonly    xyz coordinates
#natoms number of atoms

#copy 
os.system("cp "+SDir+"/*.vps "+SDir+"/*.psf .")
print "read in structure information"
fn=glob.glob(SDir+"/*.XV")
print fn


#def cutlayers(infile,nalayer,nl,nr,outfile,ord=None):
#    """
#    cut down some layers for md simulation
#    infile  input STRUCT.fdf file
#    nalayer number of atoms per layer
#    nl      nl layers from left
#    nr      nr layers from left
#    ord     atom lists in new order
#    """
geom=cutlayers(fn[0],nperlayer,cut,cut,"STRUCT.fdf")
cell, xyzonly, snr, anr = geom.pbc,geom.xyz,geom.snr,geom.anr

#xyz with element labels
#for example ['Au',0,0,0]
xyz = [[PeriodicTable[a],b[0],b[1],b[2]] for (a,b) in zip(anr,xyzonly)]

#slist are list of dynamical atoms(system) in python index
slist=N.array(range(len(xyz)))
slist=slist[nperlayer*(nl-cut):-nperlayer*(nr-cut)]
print "the following atoms are dynamic:\n"
print slist
print len(slist)

#dynxyz=N.array([xyz[i] for i in slist])



#--------------------------------------------------------------------------------------
#readin dynamical matrix and reorder it
#fn=glob.glob("PHrun/Dev*.nc")
#dyn,U,hw=ReadDynmat(fn[0],Orderw0)

#--------------------------------------------------------------------------------------
#debye bath
#number of dynamical atoms
nd=len(slist)
eta=gamma*N.identity(3*nd,N.float)

#--------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------
#def __init__(self,dt,nmd,T,syslist=None,xyz=None,harmonic=False,dyn=None,savepq=True):
print "initialise md"
#print md.__doc__
mdrun = md(dt,nmd,T,slist,xyz,nrep=nrep,npie=128)


#mdrun.setDyn(eph.DynMat)
#mdrun.SetSyslist(syslist)
#mdrun.SetHarm(True)
#mdrun.SetMD(dt,nmd)
#mdrun.SetT(T)

#-------------------------------------------------------------------------------------
#initialise siesta run
#def __init__(self,label,xyz,cell,mesh=100.,dmtol=0.001,constraints=[],lunit="Ang",eunit="eV"):
srun = siesta("mdrun",xyz,cell,mesh=mesh,dmtol=dmtol,tdir=SDir)
print srun.els
mdrun.AddSint(srun)
mdrun.sint.start()

#print potforce(None,N.zeros(3*len(slist)))
#mdrun.sint.quit()
#sys.exit()



#----------------------------------------------------
# electron bath
#def __init__(self,cats,T,dt,nmd,wmax=None,nw=None,bias=0.,efric=None,exim=None,exip=None):
#----------------------------------------------------
#class ebath:
#    """
#    cats   index of dynamical atoms connecting the bath 
#    T      equlibrium temperature of electrons
#    bias    applied bias
#    wmax    cutoff energy of electron bath
#    nw      number of energy points
#    efric   friction matrix
#    exim,exip
#    zeta1,zeta2
#
#    dt      md time step
#    nmd     md steps
#
#    noise   noise series
#    kernel  friction kernel 
#    ml      length of kernel(=1)
#
#    ebath   include ebath or not
#    """
#    def __init__(self,cats,T,dt,nmd,wmax=None,nw=None,bias=0.,\
#            efric=None,exim=None,exip=None,zeta1=None,zeta2=None):
#----------------------------------------------------
#atom indices that are connecting to electron bath
ecats=range(nd)
eb = ebath(ecats,T,mdrun.dt,mdrun.nmd,wmax=1.,nw=500,bias=0.0,efric=eta)
mdrun.AddBath(eb)

print "add electron bath"
#eb.genoi()
#PP.plot(eb.noise[:,0])
#PP.show()
#----------------------------------------------------


#---------------------------------------------------------------------------------
#MD
#------------------------------------------------------------------------------
mdrun.Run()
#------------------------------------------------------------------------------
#quit siesta
mdrun.sint.quit()
