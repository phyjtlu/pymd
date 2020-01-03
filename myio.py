#!/applications/mbrsoft/bin/python

from netCDF4 import Dataset
import numpy as N
import numpy.linalg as LA
from Inelastica import MakeGeom as MG
import sys,os
from matrix import *
from functions import *


def cutlayers(infile,nalayer,nl,nr,outfile,ord=None):
    """
    cut down some layers for md simulation
    infile  input STRUCT.fdf file
    nalayer number of atoms per layer
    nl      nl layers from left
    nr      nr layers from left
    ord     atom lists in new order
    """
    print("reading ",infile)
    geom=MG.Geom(infile)
    xyz=geom.xyz
    snr=geom.snr
    anr=geom.anr
    pbc=geom.pbc
    

    if ord is not None:
        anr,xyz=reordxyz(anr,xyz,ord)

    zs = [a[2] for a in xyz]
    olen = max(zs)-min(zs)
    
    
    na=len(xyz)
    nal = nl*nalayer
    nar = nr*nalayer
    if(nal+nar >= na):
        print("Cuting too many atoms")
        sys.exit(0)
    
    nna = int(na-nal-nar)
    nxyz = [xyz[nal+i] for i in range(nna)]
    nsnr = [snr[nal+i] for i in range(nna)]
    nanr = [anr[nal+i] for i in range(nna)]
    
    nzs = [a[2] for a in nxyz]
    nlen = max(nzs)-min(nzs)
    
    pbc[2][2]=pbc[2][2]-(olen-nlen)

    geom.xyz=nxyz
    geom.pbc=pbc
    geom.snr=nsnr
    geom.anr=nanr
    geom.natoms=nna
    
    geom.writeFDF(outfile)
    geom.writeXYZ(os.path.splitext(outfile)[0]+".xyz")

    return geom

def reordxyz(anr,xyz,ord):
    """
    anr list of atom numbers
    xyz old xyz
    ord atom list that needs to reord,siesta index
    """
    old=sorted(ord)
    nl=list(range(old[0]-1))+[i-1 for i in ord]+list(range(old[-1],len(xyz)))
    if len(nl)!=len(anr):
        print("reordxyz:length error")
        sys.exit()

    return [anr[i] for i in nl],[xyz[i] for i in nl]

def ReadEPHNCFile(filename):
    """
    Reads a NetCDF file that describes dynamical matrix, self-energies
    """
    class eph:
        pass

    file = Dataset(filename,'r')
    print('Reading from %s' % filename)

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


def ReadNewEPHNCFile(filename):
    """
    Reads a NetCDF file that describes dynamical matrix, self-energies
    """
    class eph:
        pass

    file = Dataset(filename,'r')
    print('Reading from %s' % filename)

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
    eph.zeta1=N.array(file.variables['zeta1'])
    eph.zeta2=N.array(file.variables['zeta2'])
    file.close()
    return eph

def WriteEPHNCfile(filename,wl,hw,U,DynMat,SigL,SigR,Friction,NC,NCP,zeta1,zeta2):
    """
    Write a NetCDF file contains information for harmonic analysis
    """
    fn=Dataset(filename,'w')
    print('Writing to %s' %filename)
    dhw=len(hw)
    dwl=len(wl)
    dsl=len(SigL[0])
    dsr=len(SigR[0])

    fn.createDimension('NPh',len(hw))
    fn.createDimension('NWl',len(wl))
    fn.createDimension('Nsl',len(SigL[0]))
    fn.createDimension('Nsr',len(SigR[0]))
    Write2NetCDFFile(fn,wl,'Wlist',('NWl',),units='eV')
    Write2NetCDFFile(fn,hw,'hw',('NPh',),units='eV')
    Write2NetCDFFile(fn,U,'U',('NPh','NPh',),units='None')
    Write2NetCDFFile(fn,DynMat,'DynMat',('NPh','NPh',),units='eV**2')
    Write2NetCDFFile(fn,SigL.real,'ReSigL',('NWl','Nsl','Nsl',),units='eV**2')
    Write2NetCDFFile(fn,SigL.imag,'ImSigL',('NWl','Nsl','Nsl',),units='eV**2')
    Write2NetCDFFile(fn,SigR.real,'ReSigR',('NWl','Nsl','Nsl',),units='eV**2')
    Write2NetCDFFile(fn,SigR.imag,'ImSigR',('NWl','Nsl','Nsl',),units='eV**2')
    Write2NetCDFFile(fn,Friction,'Friction',('NPh','NPh',),units='eV**2')
    Write2NetCDFFile(fn,NC,'NC',('NPh','NPh',),units='eV**2')
    Write2NetCDFFile(fn,NCP,'NCP',('NPh','NPh',),units='eV**2')
    Write2NetCDFFile(fn,zeta1,'zeta1',('NPh','NPh',),units='eV**2')
    Write2NetCDFFile(fn,zeta2,'zeta2',('NPh','NPh',),units='eV**2')
    print('Finished writing.')
    fn.close()


def Write2NetCDFFile(file,var,varLabel,dimensions,units=None,description=None):
    print('Write2NetCDFFile:', varLabel, dimensions)
    tmp = file.createVariable(varLabel,'d',dimensions)
    tmp[:] = var
    if units: tmp.units = units
    if description: tmp.description = description

def ReadNetCDFVar(file,var):
    print("ReadNetCDFFile: reading "+ var)
    f = Dataset(file,'r')
    vv=N.array(f.variables[var])
    f.close()
    return vv

def ReadMDNCFile(filename):
    """
    Reads a NetCDF file 
    """
    class mdmath:
        pass

    file = Dataset(filename,'r')
    print('Reading from %s' % filename)

    # General attributes
    mdmath.filename = filename
    mdmath.cell= N.array(file.variables['UnitCell'])
    mdmath.xyz = N.array(file.variables['XYZ'])
    mdmath.dynatom = N.array(file.variables['DynamicAtoms'])
    mdmath.atomlist= N.array(file.variables['AtomList'])

    file.close()

    return mdmath 

def ReadDynmat(filename,order=None):
    """
    Reads the NetCDF file output from PHrun, Dev*.nc
    return the dynamical matrix in real space
    
    order:  new order of atoms in siesta index
    """
    file = Dataset(filename,'r')
    print('Reading from %s' % filename)

    hw= N.array(file.variables['hw'])
    U= N.array(file.variables['U'])
    file.close()

    #reorder
    if order is not None:
        n3=3*len(order)
        if n3 != len(hw):
            print("ReadDynmat: length of order error!")
            sys.exit()
        idx=ord2idx(order)
        nU=0.0*U
        for i in range(len(idx)):
            nU[:,i]=U[:,idx[i]]
    else:
        nU=U

    dyn=mm(nU.T,N.diag(hw**2),nU)
    dyn=0.5*(dyn+dyn.T)
    return dyn,nU,hw

def ord2idx(order):
    n3=3*len(order)
    idx=N.zeros(n3,N.int)
    for i in range(len(order)):
        for j in range(3):
            idx[3*i+j]=3*(order[i]-1)+j
    return idx

def ReadSig(filename):
    """
    Reads a NetCDF file that describes dynamical matrix, self-energies
    """
    file = Dataset(filename,'r')
    print('Reading from %s' % filename)

    class eph:
        pass
    # General attributes
    eph.wl= N.array(file.variables['Wlist'])
    eph.SigL= N.array(file.variables['ReSigL'])+1j*N.array(file.variables['ImSigL'])
    eph.SigR= N.array(file.variables['ReSigR'])+1j*N.array(file.variables['ImSigR'])

    return eph



def ReadwbLambda(filename,order=None):
    """
    Reads a NetCDF file that describes dynamical matrix, self-energies
    """
    file = Dataset(filename,'r')
    print('Reading from %s' % filename)

    mus=N.array(file.variables['muLR'])
    bias=mus[0]-mus[1]
    # General attributes
    eta= N.array(file.variables['eta'])
    xim= N.array(file.variables['xim'])
    xip= N.array(file.variables['xip'])
    zeta1= N.array(file.variables['zeta1'])
    zeta2= N.array(file.variables['zeta2'])

    return bias,eta,xim,xip,zeta1,zeta2


def ReadLambda(filename,w0,order=None):
    """
    Reads a NetCDF file that describes dynamical matrix, self-energies
    """
    file = Dataset(filename,'r')
    print('Reading from %s' % filename)

    # General attributes
    wl = N.array(file.variables['wl'])
    mus=N.array(file.variables['muLR'])
    bias=mus[0]-mus[1]
    print("applied bias in Lambda.nc:",bias)
    id=nearest(w0,wl)
    print("using energy point:", wl[id])

    w00=wl[id]

    eta0= N.array(file.variables['ImPir2'][id])
    eta= -(eta0+N.transpose(eta0))/2/w00
    zeta2= -(eta0-N.transpose(eta0))/2/w00/bias
    xim0= N.array(file.variables['RePir2'][id])
    xim= -(xim0-N.transpose(xim0))/2/bias
    zeta1= (xim0+N.transpose(xim0))/2/bias

    xip= N.array(file.variables['ReLamLR'][id])
    xip= -N.pi*(xip+N.transpose(xip))/2/w00

    return bias,eta,xim,xip,zeta1,zeta2

