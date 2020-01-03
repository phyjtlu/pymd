#!/usr/bin/python

import sys, string, os, time, glob
import Scientific.IO.NetCDF as nc
import numpy as N
import numpy.linalg as LA

#from Inelastica import SiestaIO as SIO
from md import *
from phbath import *
from ebath import *
#from siesta import *
from matrix import *
from myio import *

import numpy as np
import units as U
from functions import *
from myfft import *
from optparse import OptionParser, OptionGroup

from brenner import *
#import py_brenner_module


#--------------------------------------------------------------------------------------
#misc driver routines
def Write2NetCDFFile(file,var,varLabel,dimensions,units=None,description=None):
    #print 'Write2NetCDFFile:', varLabel, dimensions
    tmp = file.createVariable(varLabel,'d',dimensions)
    tmp[:] = var
    if units: tmp.units = units
    if description: tmp.description = description

def ReadNetCDFVar(file,var):
    print("ReadNetCDFFile: reading "+ var)
    f = nc.NetCDFFile(file,'r')
    vv=N.array(f.variables[var])
    f.close()
    return vv

def sameq(q1,q2):
    if(len(q1) != len(q2)):
        #print "sameq: qq1 and qq2 not same length"
        return False 
    qq1 = N.array(q1)
    qq2 = N.array(q2)
    dif = qq1-qq2
    if(max(abs(dif)) < 10e-3):
        return True
    else:
        return False

def ReadEPHNCFile(filename):
    """
    Reads a NetCDF file that describes dynamical matrix, self-energies
    """
    class eph:
        pass

    file = nc.NetCDFFile(filename,'r')
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


def ReadMatlabConfigFile(filename):
    """
    Reads a NetCDF file that describes dynamical matrix, self-energies
    """
    class config:
        pass

    file = nc.NetCDFFile(filename,'r')
    print('Reading from %s' % filename)

    # General attributes[:,0]
    config.filename = filename
    config.R= np.array(file.variables['R'][:,0]) # In Ang
    config.R0= np.array(file.variables['R0'][:,0]) # In Ang
    config.typevec= np.array(file.variables['typevec'][:,0],int)
    config.D00ev2= np.array(file.variables['D00ev2'][:])
    config.Rfix= np.array(file.variables['Rfix'][:,0]) # In Ang
    config.typefix= np.array(file.variables['typefix'][:,0],int)
    config.cell= np.array(file.variables['cell'][:,0]) # In Ang
    file.close()

    return config

def ReadMatlabPHNCFile(filename):
    """
    Reads a NetCDF file that describes dynamical matrix, self-energies
    """
    class eph:
        pass

    file = nc.NetCDFFile(filename,'r')
    print('Reading from %s' % filename)

    # General attributes
    eph.filename = filename
    eph.T= float(file.variables['T'][:,0]) #np.array(file.variables['T'][:,0])
    eph.dT= float(file.variables['dT'][:,0]) #np.array(file.variables['dT'][:,0])
    eph.hwcut= float(file.variables['hwcut'][:,0]) #np.array(file.variables['hwcut'][:,0])
    eph.dt= float(file.variables['dt'][:,0]) #np.array(file.variables['dt'][:,0])
    eph.nw= int(file.variables['nw'][:,0]) #np.array(file.variables['nw'][:,0],int)
    eph.nmd= int(file.variables['nmd'][:,0]) #np.array(file.variables['nmd'][:,0],int)
    eph.nmemL= int(file.variables['nmemL'][:,0]) #np.array(file.variables['nmemL'][:,0],int)
    eph.nmemR= int(file.variables['nmemR'][:,0]) #np.array(file.variables['nmemR'][:,0],int)
    eph.NA_L= np.array(file.variables['NA_L'][:,0],int)
    eph.NA_R= np.array(file.variables['NA_R'][:,0],int)
    eph.idL= np.array(file.variables['idLn'][:,0],int)
    eph.idR= np.array(file.variables['idRn'][:,0],int)
    eph.E= np.array(file.variables['E'][:,0])
    eph.SigL= np.array(file.variables['SigL_Re'][:])+1j*np.array(file.variables['SigL_Im'][:])
    eph.SigR= np.array(file.variables['SigR_Re'][:])+1j*np.array(file.variables['SigR_Im'][:])
    #for iw in range(eph.nw):
    #eph.SigL[iw,1,1]=eph.SigL2[1,1,iw]
    #eph.SigR[iw,1,1]=eph.SigR2[1,1,iw]
    #N=3*eph.NA_L
    #eph.SigL=np.zeros((N,N,eph.nw),np.complex)
    #for iw in range(eph.nw):
    #    for id1 in range(N):
    #        for id2 in range(N):
        #eph.SigL[iw,:,:]= np.array(file.variables['SigL_Re'][:])+1j*np.array(file.variables['SigL_Im'][:])
        #eph.SigR[iw,:,:]= np.array(file.variables['SigR_Re'][:])+1j*np.array(file.variables['SigR_Im'][:])
    #print "was here"
    file.close()

    return eph

def WriteEPHNCfile(filename,E,nw,NA_L,NA_R,fL,fR,KernelL,KernelR,gammaL,gammaR):
    """
    Write a NetCDF file contains information for harmonic analysis
    """
    fn=nc.NetCDFFile(filename,'w')
    print('Writing to %s' %filename)
    
    #fn.createDimension('NPh',len(hw))
    fn.createDimension('NWl',len(E))
    fn.createDimension('Nsl',len(fL[:,0]))
    fn.createDimension('Nsr',len(fR[:,0]))
    fn.createDimension('Nsl2',len(fL[0,:]))
    fn.createDimension('Nsr2',len(fR[0,:]))

    fn.createDimension('Nkl',len(KernelL[:,0,0]))
    fn.createDimension('Nkr',len(KernelR[:,0,0]))
    fn.createDimension('Nkl2',len(KernelL[0,:,0]))
    fn.createDimension('Nkr2',len(KernelR[0,:,0]))
    fn.createDimension('Nkl3',len(KernelL[0,0,:]))
    fn.createDimension('Nkr3',len(KernelR[0,0,:]))

    fn.createDimension('Ngl1',len(gammaL[:,0,0]))
    fn.createDimension('Ngr1',len(gammaR[:,0,0]))
    fn.createDimension('Ngl2',len(gammaL[0,:,0]))
    fn.createDimension('Ngr2',len(gammaR[0,:,0]))
    fn.createDimension('Ngl3',len(gammaL[0,0,:]))
    fn.createDimension('Ngr3',len(gammaR[0,0,:]))

    Write2NetCDFFile(fn,fL.real,'fL',('Nsl','Nsl2',),units='dunno')
    Write2NetCDFFile(fn,fR.real,'fR',('Nsr','Nsr2',),units='dunno')
    Write2NetCDFFile(fn,KernelL.real,'KernelL',('Nkl','Nkl2','Nkl3',),units='dunno')
    Write2NetCDFFile(fn,KernelR.real,'KernelR',('Nkr','Nkr2','Nkr3',),units='dunno')
    Write2NetCDFFile(fn,gammaL.real,'gammaL',('Ngl1','Ngl2','Ngl3',),units='eV')
    Write2NetCDFFile(fn,gammaR.real,'gammaR',('Ngr1','Ngr2','Ngr3',),units='eV')

    print('Finished writing.')
    fn.close()

def setupParameters():
    # User defined input:
    global general
    usage = "usage: %prog [options] DestinationDirectory" #""#
    description = ""
    parser = OptionParser(usage,description=description)
    Paths = OptionGroup(parser, "Paths for noise generation")
    Paths.add_option("-f", "--file2load", dest='file2load', default=None,type='string',
                          help="Path to file holding MD-data.")
    Paths.add_option("-s", "--file2save", dest='file2save', default='out.nc',type='string',
                          help="Path to file where noise is saved.")
    parser.add_option_group(Paths)
    (general, args) = parser.parse_args()
#--------------------------------------------------------------------------------------
#-------------------------------Program Starts-----------------------------------------
#--------------------------------------------------------------------------------------
setupParameters()
print(general.file2load) #SystemInput
#--------------------------------------------------------------------------------------
# Load Input:
#T=300                     # Temperature (K)
#dT=60                     # Temperature difference (K)
#hwcut=0.03                # Lower cutoff frequency (eV)
#dt=8                      # Time step
#nw=50                     # Number of energy points Sigma was calculated for
#nmd=100                   # Number of md steps for which noise should be generated
#nmem=50                   # Length of memory (in # of timesteps)
#NA=2                      # Number of lead atoms
#dim=3                     # 3 dimensions
#N=NA*dim                  # Number of degree of freedom
#idL=np.zeros((NA));       # Index of central system atoms connecting to the bath (python index)
#E=np.linspace(hwcut,2,nw) # Energy points where selfenergy is calculated
#Sig=tensor[iE,idim,idim]  # Tensor with selfenergy (energy points,degree of freedom,degree of freedom)
#filename="GraphenePH.nc"
eph=ReadMatlabPHNCFile(general.file2load)
dim=3
NL=eph.NA_L*dim
NR=eph.NA_R*dim

# Config and dynamical matrix:
#filename="GrapheneConfig.nc"
config=ReadMatlabConfigFile(general.file2load)
print(config.typevec.shape,config.R.shape,config.R0.shape,config.D00ev2.shape,eph.SigL.shape,eph.E.shape)
#--------------------------------------------------------------------------------------
# Noise generation:
#eph.nmd=200 #200
eph.nmd=2**6
#print 'TgTest',eph.nmemL,eph.nmemR
#eph.nmemL=150
#eph.nmemR=150
#print 'TgTest',eph.nmemL,eph.nmemR
#eph.hwcut=0.22
##eph.hwcut=0.5
eph.T=500#3500
eph.dT=0#.1*2*eph.T #60
#print 'dt',eph.dt
#eph.dt=0.5/(0.658211814201041) # from fs to JT units
eph.dt=eph.dt # 0.5->0.05 fs due to hydrogen!
#print 'dt',eph.dt
eph.eta_ad=0.005#.01#0.000000001#

phl = phbath(eph.T+eph.dT/2,list(range(eph.NA_L)),eph.hwcut/2,eph.nw,eph.dt,eph.nmd,ml=eph.nmemL,sig=eph.SigL,gwl=eph.E,eta_ad=eph.eta_ad,classical=False,zpmotion=True)##sig=eph.SigL,gwl=eph.wl)
phl.gmem()
#phl.gnoi()
phr = phbath(eph.T-eph.dT/2,list(range(eph.NA_R))-eph.NA_R,eph.hwcut/2,eph.nw,eph.dt,eph.nmd,ml=eph.nmemR,sig=eph.SigR,gwl=eph.E,eta_ad=eph.eta_ad, classical=False,zpmotion=True)##sig=eph.SigL,gwl=eph.wl)
phr.gmem()
#phr.gnoi()
#classical=False: to use 2*nu*kB*T or full,zpmotion=False:To include QM +1/2 or neglect.

#--------------------------------------------------------------------------------------
# MDrun:
# What is xyz? I input actual positions but maybe it should be displacement=config.R-config.R0?
# config.typevec = types of atoms in system = vector with ones here.
# config.D00ev2 = dynamical matrix in eV^2.
curcof = 243414.
nrep = 4
dt = eph.dt
nmd = eph.nmd
hwcut = eph.hwcut
T=eph.T
#def __init__(self,dt,nmd,T,syslist=None,xyz=None,harmonic=False,\
#             dyn=None,savepq=True,nrep=1,npie=8):
#mdrun = md(dt,nmd,T,syslist=config.typevec,xyz=config.R,harmonic=True,dyn=config.D00ev2,nrep=nrep,npie=1)
Harmonic=0
if Harmonic==1:
   mdrun = md(dt,nmd,T,harmonic=True,dyn=config.D00ev2,nrep=nrep,npie=1)
else:
   #   PeriodicTable = {'H':1,1:'H','D':1001,1001:'D','He':2,2:'He','Li':3,3:'Li','Be':4,4:'Be','B':5,5:'B','C':6,6:'C','N':7,7:'N','O':8,8:'O','F':9,9:'F','Ne':10,10:'Ne','Na':11,11:'Na','Mg':12,12:'Mg','Al':13,13:'Al','Si':14,14:'Si','P':15,15:'P','S':16,16:'S','Cl':17,17:'Cl','Ar':18,18:'Ar','K':19,19:'K','Ca':20,20:'Ca','Sc':21,21:'Sc','Ti':22,22:'Ti','V':23,23:'V','Cr':24,24:'Cr','Mn':25,25:'Mn','Fe':26,26:'Fe','Co':27,27:'Co','Ni':28,28:'Ni','Cu':29,29:'Cu','Zn':30,30:'Zn','Ga':31,31:'Ga','Ge':32,32:'Ge','As':33,33:'As','Se':34,34:'Se','Br':35,35:'Br','Kr':36,36:'Kr','Rb':37,37:'Rb','Sr':38,38:'Sr','Y':39,39:'Y','Zr':40,40:'Zr','Nb':41,41:'Nb','Mo':42,42:'Mo','Tc':43,43:'Tc','Ru':44,44:'Ru','Rh':45,45:'Rh','Pd':46,46:'Pd','Ag':47,47:'Ag','Cd':48,48:'Cd','In':49,49:'In','Sn':50,50:'Sn','Sb':51,51:'Sb','Te':52,52:'Te','I':53,53:'I','Xe':54,54:'Xe','Cs':55,55:'Cs','Ba':56,56:'Ba','La':57,57:'La','Ce':58,58:'Ce','Pr':59,59:'Pr','Nd':60,60:'Nd','Pm':61,61:'Pm','Sm':62,62:'Sm','Eu':63,63:'Eu','Gd':64,64:'Gd','Tb':65,65:'Tb','Dy':66,66:'Dy','Ho':67,67:'Ho','Er':68,68:'Er','Tm':69,69:'Tm','Yb':70,70:'Yb','Lu':71,71:'Lu','Hf':72,72:'Hf','Ta':73,73:'Ta','W':74,74:'W','Re':75,75:'Re','Os':76,76:'Os','Ir':77,77:'Ir','Pt':78,78:'Pt','Au':79,79:'Au','Hg':80,80:'Hg','Tl':81,81:'Tl','Pb':82,82:'Pb','Bi':83,83:'Bi','Po':84,84:'Po','At':85,85:'At','Rn':86,86:'Rn','Fr':87,87:'Fr','Ra':88,88:'Ra','Ac':89,89:'Ac','Th':90,90:'Th','Pa':91,91:'Pa','U':92,92:'U','Np':93,93:'Np','Pu':94,94:'Pu','Am':95,95:'Am','Cm':96,96:'Cm','Bk':97,97:'Bk','Cf':98,98:'Cf','Es':99,99:'Es','Fm':100,100:'Fm','Md':101,101:'Md','No':102,102:'No'}
   typearray=[6,1]
   typevec2=[]
   R0=[]
   R=[]
   for it in range(len(config.typevec)):
       if config.typevec[it] == 1:
           typevec2.append(6)
           for j in range(3):
              R0.append(config.R0[3*it+j])
              R.append(config.R[3*it+j])
	
   #anr=np.array(config.typevec*6,dtype='i') # This should be a numpy array for vec*6 to work!
   anr=np.array(typevec2,dtype='i') # This should be a numpy array for vec*6 to work!
   config.R0=np.array(R0)
   config.R=np.array(R)
   config.D00ev2=np.identity(len(typevec2)*3)
   print(config.R.shape,config.R0.shape,anr.shape)


   #the new structure without H
   f = open("STRUCT.noH.xyz","w")
   f.write("%s\n\n"%(len(anr)))
   for i in range(len(anr)):
       f.write("C  ")
       for j in range(3):
           f.write("  %f"%(config.R[3*i+j]))
       f.write(" \n")
   f.close()



   #Rtemp=np.reshape(config.R,[3,int(len(config.R)/3)])
   #els=[[PeriodicTable[a]] for (a) in (anr)]
   #print 'Check R: ',len(config.R0)
   #cell=np.array([3,3,100],dtype='d')		# Unitcell dimensions.
   #R=np.array([0,0,0,0,1,0,1,0,0],dtype='d')	# Position vector.
   #typevec=np.array([6,6,6],dtype='i')		# type of atoms (atomic numbers).
   #print config.Rfix,config.typefix,config.cell
   #brennerrun=brenner(xyz=R,anr=typevec,cell=cell,devicedir=0)#,constrained,anr_constrained)
   #config.cell[2]=39.3522
   brennerrun=brenner(xyz=config.R0,anr=anr,cell=config.cell,devicedir=3,constrained=config.Rfix,anr_constrained=config.typefix*6)
   brennerrun.initforce()
   print('First Force: ',brennerrun.f0)#,brennerrun.f0[0:6]#config.R[0:6]
   brennerrun.f0=N.zeros(len(config.R0))
   mdrun = md(dt,nmd,T,harmonic=False,dyn=config.D00ev2,nrep=nrep,npie=1,savepq=True)#savepq=False
   mdrun.brennerrun=brennerrun

mdrun.AddBath(phl)
mdrun.AddBath(phr)
#mdrun.savepq=1
mdrun.initialise()
mdrun.Run()


# # Check force directly:
#config.cell[0]=config.cell[0]+100
#newR=N.concatenate((config.R,config.Rfix))
#anr_b=N.concatenate((anr,np.array(config.typefix*6,dtype='i')))
# #u1,f,u3=py_brenner_module.brennerf(config.cell,anr,config.R)
#u1,f,u3=py_brenner_module.brennerf(config.cell,anr_b[:],newR[:])
#print config.cell,anr[0],f[0:6]
#print 'R: ', newR[1194::]#[0:6]#[1195::]

# Output:
# kappa.dat: IL,IR (has to devide by DeltaT to get ThermalConductance and manually average over all MD sequences)
# MDx.nc:
#        double power(nmd, two) ; vv-correlation (averaged over all MD runs-> only need the last .nc file).
#        double p(nph) ; initial positions
#        double q(nph) ; initial velocities
#        double t(one) ; times
#        double ipie(one) ;
#        double phis(mem, nph) ; history of positions
#        double qhis(mem, nph) ; history of velocities
#        double ps(nmd, nph) ; position trajectories
#        double qs(nmd, nph) ; velocity trajectories
#        double noise0(nmd, n0) ; 0-> Left noise
#        double noise1(nmd, n1) ; 1-> Right noise
