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
#user input parameters
#print 'PhononNetCDF : Calculate MAMA using HSSigma and Heph NetCDF file'
args = sys.argv[1:]
if len(args) != 3:
    print "Usage : /usr/bin/python %s <T> <dT> <eV>"%(sys.argv[0])
    #print 'PhononNetCDF : Calculate MAMA using HSSigma and Heph NetCDF file'
    sys.exit()
else:
    T=float(args[0])
    dT=float(args[1])
    eV=float(args[2])
    print "Average Temperature: %s\n"%T
    print "Temperature Difference: %s\n"%(dT)
    print "Applied bias (muL-muR): %s\n"%(eV)

#--------------------------------------------------------------------------------------
curcof = 243414.
nrep = 8
dt = 8
nmd = 2**16
hwcut=0.03
#-------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
#readin the phonon self-energies
eph=ReadSig("Sig.nc")
#--------------------------------------------------------------------------------------
import Inelastica.SiestaIO as SIO
#read in the relaxed structure
#cell   3x3 unit vectors
#snr    species numbers defined in fdf
#anr    atomic numbers 
#xyzonly    xyz coordinates
#natoms number of atoms

print "read in structure information"
fn=glob.glob("CGrun/*.XV")
print fn

#fullstru="SIESTASTRUCT.fdf"
#geom2geom fn fullstru 1 1 1
nperlayer=25
nl=6
nr=7
cut=nl-2
#siesta order
#Five atom chain,generated from Mads mathematica code
Orderw = [5, 2, 4, 3, 1, 6, 7, 8, 13, 10, 11, 12, 9]
Orderw =[nperlayer*nl+i for i in Orderw]
geom=cutlayers(fn[0],nperlayer,cut,cut,Orderw,"STRUCT.fdf")
cell, xyzonly, snr, anr = geom.pbc,geom.xyz,geom.snr,geom.anr

#xyz with element labels
#for example ['Au',0,0,0]
xyz = [[PeriodicTable[a],b[0],b[1],b[2]] for (a,b) in zip(anr,xyzonly)]

#slist are list of dynamical atoms(system) in python index
slist=N.array(range(len(xyz)))
slist=slist[nperlayer*(nl-cut):-nperlayer*(nr-cut)-1]
print "the following atoms are dynamic:\n"
print slist
print len(slist)

#--------------------------------------------------------------------------------------
#def __init__(self,dt,nmd,T,syslist=None,xyz=None,harmonic=False,dyn=None,savepq=True):
print "initialise md"
#print md.__doc__
fn=glob.glob("PHrun/Dev*.nc")
dyn=ReadDynmat(fn[0])
mdrun = md(dt,nmd,T,slist,xyz,harmonic=True,dyn=dyn,nrep=nrep)
print mdrun.dyn


#mdrun.setDyn(eph.DynMat)
#mdrun.SetSyslist(syslist)
#mdrun.SetHarm(True)
#mdrun.SetMD(dt,nmd)
#mdrun.SetT(T)

#-------------------------------------------------------------------------------------
##initialise siesta run
##def __init__(self,label,xyz,cell,mesh=100.,dmtol=0.001,constraints=[],lunit="Ang",eunit="eV"):
#constraints=[[1,25],[-75,-1]]    
#srun = siesta("mdrun",xyz,cell,constraints=constraints)
#print srun.els
#mdrun.AddSint(srun)
#mdrun.sint.start()
#
##print potforce(None,N.zeros(3*len(slist)))
##mdrun.sint.quit()
##sys.exit()

#---------------------------------------------------------------------------------
#debye bath
#---------------------------------------------------------------------------------
#def __init__(self,T,cats,debye,nw,dt=None,nmd=None,ml=None,mcof=2.0,\
        #        gamma=None,gwl=None,K00=None,K01=None,V01=None):
# initialise a phonon bath with debye damping
# phl = phbath(300,[0,1,2,3,4,5,6,7,8,9,10,11,12],0.022,100,mdrun.dt,mdrun.nmd)
# phl.gnoi()
# phl.gmem()
# PP.plot(phl.noise[:,0])
# PP.show()
# ##add the bath 
# mdrun.AddBath(phl)
# print "memory length: %s\n"%mdrun.ml
# 
#---------------------------------------------------------------------------------


#----------------------------------------------------
# electron bath
#def __init__(self,cats,T,dt,nmd,wmax=None,nw=None,bias=0.,efric=None,exim=None,exip=None):
#----------------------------------------------------
###ecats=range(13)
###eb = ebath(ecats,T,mdrun.dt,mdrun.nmd,wmax=1.0,nw=500,bias=eV,efric=eph.efric,exim=eph.xim,exip=eph.xip)
###mdrun.AddBath(eb)
#eb.genoi()
#PP.plot(eb.noise[:,0])
#PP.show()
#----------------------------------------------------



#---------------------------------------------------------------------------------
# left phonon bath
#---------------------------------------------------------------------------------
#    def __init__(self,T,cats,debye,nw,dt=None,nmd=None,ml=None,mcof=2.0,gamma=None,gwl=None,K00=None,K01=None,V01=None):
    #        self.T,self.debye,self.cats = T,debye,N.array(cats,dtype='int')
# initialise a phonon bath with debye damping
#cats=range(16)
#phl = phbath(T+dT/2,cats,hwcut,200,mdrun.dt,mdrun.nmd,ml=100)
#phl.gmem()
#mdrun.AddBath(phl)

#cats: python index of system atoms that are connecting to the left phonon bath
cats=range(mdrun.syslist[0],mdrun.syslist[0]+5)
cats=[0,1,2,3,4]
print cats
phl = phbath(T+dT/2,cats,hwcut,200,mdrun.dt,mdrun.nmd,ml=100,sig=eph.SigL,gwl=eph.wl)
phl.gmem()
mdrun.AddBath(phl)

#print phl.kernel.shape
#PP.plot(phl.noise[:,0])
#PP.savefig("noise_left.pdf",format="pdf")
#PP.close()
#PP.plot([N.trace(a) for a in phl.kernel])
#PP.savefig("friction_left.pdf",format="pdf")
#PP.close()

#print "memory length: %s\n"%mdrun.ml
#print phl.ml
#---------------------------------------------------------------------------------

#----------------------------------------------------
# right phonon bath
#----------------------------------------------------
#    def __init__(self,T,cats,debye,nw,dt=None,nmd=None,ml=None,mcof=2.0,gamma=None,gwl=None,K00=None,K01=None,V01=None):
    #        self.T,self.debye,self.cats = T,debye,N.array(cats,dtype='int')
# initialise a phonon bath with debye damping
#phl = phbath(300,[0,1,2,3,4],0.022,200,mdrun.dt,mdrun.nmd,ml=100,gamma=gaml,gwl=eph.wl)
###phr = phbath(T-dT/2,[8,9,10,11,12],hwcut,200,mdrun.dt,mdrun.nmd,ml=100,sig=eph.SigR,gwl=eph.wl)
###phr.gmem()
###mdrun.AddBath(phr)
###
#print phr.kernel.shape

#cats: python index of system atoms that are connecting to the right phonon bath
cats=range(mdrun.syslist[-1]-4,mdrun.syslist[-1]+1)
cats=[8,9,10,11,12]
print cats
phr = phbath(T-dT/2,cats,hwcut,200,mdrun.dt,mdrun.nmd,ml=100,sig=eph.SigR,gwl=eph.wl)
phr.gmem()
mdrun.AddBath(phr)

#PP.plot(phr.noise[:,0])
#PP.savefig("noise_right.pdf",format="pdf")
#PP.close()
#PP.plot([N.trace(a) for a in phr.kernel])
#PP.savefig("friction_right.pdf",format="pdf")
#PP.close()
#add the bath 
#mdrun.AddBath(phr)
#print "memory length: %s\n"%mdrun.ml
#print phr.ml
#---------------------------------------------------------------------------------



#---------------------------------------------------------------------------------
#MD
#---------------------------------------------------------------------------------
mdrun.Run()
#------------------------------------------------------------------------------
#quit siesta
#mdrun.sint.quit()
