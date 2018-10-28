#!/applications/mbrsoft/bin/python
#import matplotlib.pyplot as PP
import sys, string, os, time, glob
import Scientific.IO.NetCDF as nc
from md import *
from phbath import *
from ebath import *
from siesta import *


#-------------------------------------------------------------------------------------
#user input parameters
#print 'Usage : /usr/bin/python thermal.py  <T> <dT>'
#print 'PhononNetCDF : Calculate MAMA using HSSigma and Heph NetCDF file'
args = sys.argv[1:]
if len(args) != 3:
    print 'Usage : /usr/bin/python thermal.py  <T> <dT>'
    print 'PhononNetCDF : Calculate MAMA using HSSigma and Heph NetCDF file'
    sys.exit()
else:
    T=float(args[0])
    dT=float(args[1])
    eV=float(args[2])
    print "Average Temperature: %s\n"%T
    print "Temperature Difference: %s\n"%(dT)
    print "Applied bias: %s\n"%(eV)

#--------------------------------------------------------------------------------------
curcof = 243414.
nrep = 1
dt = 8
nmd = 2**12
hwcut=0.03
#-------------------------------------------------------------------------------------

def ReadEPHNCFile(filename):
    """
    Reads a NetCDF file that describes dynamical matrix, self-energies
    """
    class eph:
        pass

    file = nc.NetCDFFile(filename,'r')
    print 'Reading from %s' % filename

    # General attributes
    eph.filename = filename
    eph.wl= N.array(file.variables['Wlist'])
    eph.hw= N.array(file.variables['hw'])
    eph.U= N.array(file.variables['U'])
    eph.DynMat= N.array(file.variables['DynMat'])
    eph.SigL= N.array(file.variables['ReSigL'])+1j*N.array(file.variables['ImSigL'])
    eph.SigR= N.array(file.variables['ReSigR'])+1j*N.array(file.variables['ImSigR'])
    eph.efric=N.array(file.variables['Friction'])
    eph.xim=N.array(file.variables['NC'])
    eph.xip=N.array(file.variables['NCP'])

    file.close()

    return eph


def ReadMDNCFile(filename):
    """
    Reads a NetCDF file 
    """
    class mdmath:
        pass

    file = nc.NetCDFFile(filename,'r')
    print 'Reading from %s' % filename

    # General attributes
    mdmath.filename = filename
    mdmath.cell= N.array(file.variables['UnitCell'])
    mdmath.xyz = N.array(file.variables['XYZ'])
    mdmath.dynatom = N.array(file.variables['DynamicAtoms'])
    mdmath.atomlist= N.array(file.variables['AtomList'])

    file.close()

    return mdmath 

#-------------------------------------------------------------------------------------
#eph=ReadEPHNCFile("./EPH.nc")
#mdmath = ReadMDNCFile("./MD1.nc")
#xyz = [["Au",a[0],a[1],a[2]] for a in mdmath.xyz]          
#slist=N.array(range(13))+32 #python indices
#print slist


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
import Inelastica.SiestaIO as SIO
cell, xyz, snr, anr, natoms=SIO.ReadFDFFile("STRUCT.fdf")

xyz = [[PeriodicTable[a],b[0],b[1],b[2]] for (a,b) in zip(anr,xyz)]          
print xyz
slist=N.array(range(len(xyz)))
print slist
#--------------------------------------------------------------------------------------
#    def __init__(self,dt,nmd,T,syslist=None,xyz=None,harmonic=False,dyn=None,savepq=True):
print "initialise md"
print md.__doc__
test = md(dt,nmd,T,slist,xyz,harmonic=False,dyn=N.zeros((3*natoms,3*natoms)),nrep=1)
print test.xyz


#test.setDyn(eph.DynMat)
#syslist=range(len(eph.DynMat)/3)
#test.SetSyslist(syslist)
#test.SetHarm(True)
#test.SetMD(dt,nmd)
#test.SetT(T)

#-------------------------------------------------------------------------------------
#initialise siesta run
#def __init__(self,label,xyz,cell,mesh=100.,dmtol=0.001,constraints=[],lunit="Ang",eunit="eV"):
srun = siesta("test",xyz,cell,constraints=[])
print srun.els
test.AddSint(srun)
test.sint.start()

#print potforce(None,N.zeros(3*len(slist)))
#test.sint.quit()
#sys.exit()

#---------------------------------------------------------------------------------
#debye bath
#---------------------------------------------------------------------------------
#def __init__(self,T,cats,debye,nw,dt=None,nmd=None,ml=None,mcof=2.0,\
        #        gamma=None,gwl=None,K00=None,K01=None,V01=None):
# initialise a phonon bath with debye damping
# phl = phbath(300,[0,1,2,3,4,5,6,7,8,9,10,11,12],0.022,100,test.dt,test.nmd)
# phl.gnoi()
# phl.gmem()
# PP.plot(phl.noise[:,0])
# PP.show()
# ##add the bath 
# test.AddBath(phl)
# print "memory length: %s\n"%test.ml
# 
#---------------------------------------------------------------------------------


#----------------------------------------------------
# electron bath
#def __init__(self,cats,T,dt,nmd,wmax=None,nw=None,bias=0.,efric=None,exim=None,exip=None):
#----------------------------------------------------
###ecats=range(13)
###eb = ebath(ecats,T,test.dt,test.nmd,wmax=1.0,nw=500,bias=eV,efric=eph.efric,exim=eph.xim,exip=eph.xip)
###test.AddBath(eb)
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
cats=range(16)
phl = phbath(T+dT/2,cats,hwcut,200,test.dt,test.nmd,ml=100)
phl.gmem()
test.AddBath(phl)

#print phl.kernel.shape
#PP.plot(phl.noise[:,0])
#PP.savefig("noise_left.pdf",format="pdf")
#PP.close()
#PP.plot([N.trace(a) for a in phl.kernel])
#PP.savefig("friction_left.pdf",format="pdf")
#PP.close()

#print "memory length: %s\n"%test.ml
#print phl.ml
#---------------------------------------------------------------------------------

#----------------------------------------------------
# right phonon bath
#----------------------------------------------------
#    def __init__(self,T,cats,debye,nw,dt=None,nmd=None,ml=None,mcof=2.0,gamma=None,gwl=None,K00=None,K01=None,V01=None):
    #        self.T,self.debye,self.cats = T,debye,N.array(cats,dtype='int')
# initialise a phonon bath with debye damping
#phl = phbath(300,[0,1,2,3,4],0.022,200,test.dt,test.nmd,ml=100,gamma=gaml,gwl=eph.wl)
###phr = phbath(T-dT/2,[8,9,10,11,12],hwcut,200,test.dt,test.nmd,ml=100,sig=eph.SigR,gwl=eph.wl)
###phr.gmem()
###test.AddBath(phr)
###
#print phr.kernel.shape

#PP.plot(phr.noise[:,0])
#PP.savefig("noise_right.pdf",format="pdf")
#PP.close()
#PP.plot([N.trace(a) for a in phr.kernel])
#PP.savefig("friction_right.pdf",format="pdf")
#PP.close()
#add the bath 
#test.AddBath(phr)
#print "memory length: %s\n"%test.ml
#print phr.ml
#---------------------------------------------------------------------------------



#---------------------------------------------------------------------------------
#MD
#---------------------------------------------------------------------------------
test.Run()
#------------------------------------------------------------------------------
#quit siesta
test.sint.quit()
