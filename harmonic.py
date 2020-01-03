#!/applications/mbrsoft/bin/python

import sys,time
import os.path
import numpy as N
from numpy import linalg as LA
import Scientific.IO.NetCDF as nc

import units as U
from matrix import *
from functions import *
from myfft import *
from noise import *
from spectrum import *

#phonon noise in w space
#def phnoisew(gamma,wl,T,phcut):
#electron noise in w space 
#def enoisew(wl,efric,exim,exip,bias,T,ecut):

def Drs(wl,dyn,se):
    """
    calculate the harmonic spectrum using green's functions
    
    wl:     list of frequencies
    dyn:    dynamical matrix
    se:     retarded self-energy
    """
    se=N.array(se)
    wl=N.array(wl)
    n=len(dyn)
    m=len(wl)
    D=N.zeros((m,n,n),N.complex)
    for i in range(m):
        D[i]=Dr(wl[i],dyn,se[i])
    return D
def Dr(w,dyn,se):
    eta=10**-10
    n=len(dyn)
    return LA.inv((w+1.j*eta)**2*N.identity(n)-dyn-se)

def DDos(wl,Ds):
    """
    calculate the harmonic spectrum using green's functions
    """
    nc=chkShape(Ds[0])
    nw=len(wl)
    #As=N.zeros((nw,nc,nc))
    As=N.zeros((nw))
    for i in range(len(wl)):
        As[i] = -2.*wl[i]*N.trace(Ds[i].imag)
    return As

def UU(wl,Ds,noise):
    """
    semi-classical displacement correlation function
    UU(w)=Dr(w)*Pi(w)*Da(w)

    Ds:     retarded green's function
    noise:  noise correlation function
    """
    n=len(Ds)
    Ds=N.array(Ds)
    #Us=0.0*Ds
    Us=N.zeros(n)
    for i in range(n):
        Us[i]=N.trace(mm(Ds[i],noise[i],dagger(Ds[i]))).real
    return Us

def PP(wl,Ds,noise):
    """
    semi-classical displacement correlation function
    PP(w)=w**2*Dr(w)*Pi(w)*Da(w)

    Ds:     retarded green's function
    noise:  noise correlation function
    """
    n=len(Ds)
    Ds=N.array(Ds)
    #Us=0.0*Ds
    Us=N.zeros(n)
    for i in range(n):
        Us[i]=wl[i]**2*N.trace(mm(Ds[i],noise[i],dagger(Ds[i]))).real
    return Us

def SigE(wl,efric,exim,zeta1,zeta2,bias):
    """
    effective retarded self-energy from electrons
    """
    nw=len(wl)
    nc=chkShape(efric)
    Sig=N.zeros((nw,nc,nc),N.complex)

    for i in range(nw):
        Sig[i]=-1.j*wl[i]*(efric+bias*zeta2)+bias*zeta1-bias*exim
    return Sig




#--------------------------------------------------------------------------------------
#testing
#
#--------------------------------------------------------------------------------------

if __name__=="__main__":
    from myio import *

    eph=ReadNewEPHNCFile("EPH.nc")

    #def enoisew(wl,efric,exim,exip,bias,T,ecut):
    bias=0.0
    T=0.0
    ecut=0.05
    phcut=0.05
    enoise=enoisew(eph.wl,eph.efric,eph.xim,eph.xip,bias,T,ecut)

    nw=len(eph.wl)
    nc=chkShape(eph.efric)
    nl=chkShape(eph.SigL[0])
    nr=chkShape(eph.SigR[0])
    print(nl,nr)

    phSig=N.zeros((nw,nc,nc),N.complex)
    phnoise=N.zeros((nw,nc,nc),N.complex)

    #padding
    phSig[:,0:nl,0:nl]=phSig[:,0:nl,0:nl]+eph.SigL
    phSig[:,-nr:,-nr:]=phSig[:,-nr:,-nr:]+eph.SigR

    #---------------
    #Be careful:
    #gamma=-Im[Sig]/w
    #---------------
    gammaL=N.zeros((nw,nl,nl))
    gammaR=N.zeros((nw,nr,nr))
    for i in range(nw):
        gammaL[i]=-eph.SigL[i].imag/eph.wl[i]
        gammaR[i]=-eph.SigR[i].imag/eph.wl[i]
    phnoiseL=phnoisew(gammaL,eph.wl,T,phcut)
    phnoiseR=phnoisew(gammaR,eph.wl,T,phcut)
    phnoise[:,0:nl,0:nl]=phnoise[:,0:nl,0:nl]+phnoiseL
    phnoise[:,-nr:,-nr:]=phnoise[:,-nr:,-nr:]+phnoiseR

    #for i in range(nw):
    #    phSig[i,0:nl,0:nl]=phSig[i,0:nl,0:nl]+eph.SigL[i]
    #    phSig[i,-nr:,-nr:]=phSig[i,-nr:,-nr:]+eph.SigR[i]
    #    phnoise[i,0:nl,0:nl]=phnoise[i,0:nl,0:nl]+phnoiseL[i]
    #    phnoise[i,-nr:,-nr:]=phnoise[i,-nr:,-nr:]+phnoiseR[i]

    
    eSig=SigE(eph.wl,eph.efric,eph.xim,eph.zeta1,eph.zeta2,bias)
    Ds=Drs(eph.wl,eph.DynMat,phSig+eSig)

    DOS=DDos(eph.wl,Ds)
    DD=PP(eph.wl,Ds,enoise+phnoise)

    f = open("pp.dat","w")
    f.write("#PP col. gives the w**2*|x(w)|^2*dw/(2pi). The direct sum of this\
            column gives 2*kinetic energy of the system. Note that the factor\
            dw/(2pi)~~1/delta(0) has already been included in the noise\
            calculation automatically.\n")
    f.write("#w, DoS, PP\n")
    for i in range(len(eph.wl)):
        if eph.wl[i]<=max(ecut,phcut):
            f.write("%f  %f  %f\n"%(eph.wl[i],DOS[i],DD[i]))
        else:
            break
    f.close()

