#!/applications/mbrsoft/bin/python

import sys
import numpy as N
from numpy import linalg as LA

import units as U
from matrix import *
from functions import *
from myfft import *


#--------------------------------------------------------------------------------------

#phonon noise in w space
def phnoisew(gamma,wl,T,phcut,classical=False,zpmotion=True):
    """
    gamma       gamma[w]=-Im[\Pi^r[w]]/w
    wl          energy points of gamma (only the positive part is needed)
    T           temperature of phonon bath
    phcut       cutoff energy of phonon bath

    The noise spectrum in w space is: 2 w gamma(w) (bose(w,T)+0.5)
    20180816: equ() already gives 2*w*(bose(w,T)+zpmotion)
    """
    n=len(wl)
    gamma=N.array(gamma)
    noiw=0.0*gamma

    for i in range(n):
        #noiw[i] = equ(wl[i],phcut,T,classical,zpmotion)*wl[i]*gamma[i]
        #20180816 - equ() already has w included.
        noiw[i] = equ(wl[i],phcut,T,classical,zpmotion)*gamma[i]
    return noiw


#phonon noise generator
def phnoise(gamma,wl,T,phcut,dt,nmd,classical=False,zpmotion=True):
    """
    gamma       gamma[w]=-Im[\Pi^r[w]]/w
    wl          energy points of gamma (only the positive part is needed)
    T           temperature of phonon bath
    phcut       cutoff energy of phonon bath
    dt          time step of md 
    nmd         number of md steps


    The noise spectrum in w space is: 2 w gamma(w) (bose(w,T)+0.5)
    """
    print "phnoise: generate phonon noise from the correlation function"

    hlen = int(nmd/2)
    dw = 2.0*N.pi/dt/nmd
    delta = dt*nmd #Dirac delta delta

    nx,ny = N.shape(N.array(gamma[0]))


    #positive frequency
    phnoi1 = []
    for i in range(hlen+1):
        w = dw*i
        #flinterp: linear interpolation of gamma
        #equ(w) = 2.0*(bose(w,T)+0.5)
        amat = delta*equ(w,phcut,T,classical,zpmotion)*flinterp(w,wl,gamma)
        amat = hermitianize(amat)
        #eigh only takes half of the matrix, 
        #always return real eigenvalues and vectors
        av,au = LA.eigh(amat)
        #generate multivariance gaussian random numbers
        phnoi1.append(vargau(av,au))
    phnoi1 = N.array(phnoi1)

    #negative frequency
    phnoi2 = []
    for i in range(hlen):
        phnoi2.append(N.conjugate(phnoi1[hlen-i,:]))
    phnoi2 = N.array(phnoi2)

    phnoi1 = N.delete(phnoi1,-1,0)    #delete last row
    phnoi = N.concatenate((phnoi1,phnoi2),axis=0)


    phnoit = []
    fti = myfft(dt,nmd)
    for i in range(nx):
        phnoit.append(fti.iFourier1D(phnoi[:,i])) #w->t
    return N.transpose(N.array(phnoit))

#electron noise in w space 
def enoisew(wl,efric,exim,exip,bias,T,ecut,classical=False,zpmotion=True):
    """
    NO SHAPE CHECKING HERE

    shape of efric,exim,exip (md.nph,md.nph)

    wl          list of frequencies(energies) where the noise is calculated
    efric       electronic friction coefficient
    exim,exip   Im,Re[MA_LMA_R]
    bias        muL-muR
    T           temperature
    ecut        cutoff energy
    """
    nw=len(wl)
    #shape checking 
    nc=chkShape(efric)
    nm=chkShape(exim)
    np=chkShape(exip)

    if nc != nm or nc!=np: 
        print "enoisew: efric shape error!"
        stoppp

    enoi1=N.zeros((nw,nc,nc),N.complex)
    for i in range(nw):
        w = wl[i]
        #equilibrium part
        aw = equ(w,ecut,T,classical,zpmotion)
        amate = aw*N.array(efric)
        #nonequilibrium minus part
        awm = equ(U.hbar*w-bias,ecut,T,classical,zpmotion)
        amatm = -0.5*aw*N.array(exip)+0.5*awm*(N.array(exip)+1j*N.array(exim)) #yes,+
        #nonequilibrium plus part
        awp = equ(U.hbar*w+bias,ecut,T,classical,zpmotion)
        amatp = -0.5*aw*N.array(exip)+0.5*awp*(N.array(exip)-1j*N.array(exim)) #yes,-

        amat = amate+amatm+amatp
        enoi1[i]= hermitianize(amat)
    return enoi1




#electron noise generator
def enoise(efric,exim,exip,bias,T,ecut,dt,nmd,classical=False,zpmotion=True):
    """
    NO SHAPE CHECKING HERE

    shape of efric,exim,exip (md.nph,md.nph)

    efric       electronic friction coefficient
    exim,exip
    bias        
    T           temperature
    ecut        cutoff energy
    dt          md timestep
    nmd         total md steps
    """
    nx,ny = N.shape(N.array(efric))
    
    hlen = int(nmd/2)
    dw = 2.0*N.pi/dt/nmd
    delta = dt*nmd #Dirac delta delta

    enoi1 = []
    for i in range(hlen+1):
        w = dw*i
        #equilibrium part
        aw = delta*equ(w,ecut,T,classical,zpmotion)
        amate = aw*N.array(efric)
        #nonequilibrium minus part
        awm = delta*equ(U.hbar*w-bias,ecut,T,classical,zpmotion)
        amatm = -0.5*aw*N.array(exip)+0.5*awm*(N.array(exip)+1j*N.array(exim))
        #nonequilibrium minus part
        awp = delta*equ(U.hbar*w+bias,ecut,T,classical,zpmotion)
        amatp = -0.5*aw*N.array(exip)+0.5*awp*(N.array(exip)-1j*N.array(exim))

        amat = amate+amatm+amatp
        amath = hermitianize(amat)
        #eigh only takes half of the matrix, always return real eigenvalues and
        #vectors
        av,au = LA.eigh(amath)
        #generate multivariance gaussian random numbers
        enoi1.append(vargau(av,au))
    enoi1 = N.array(enoi1)

    enoi2 = []
    for i in range(hlen):
        enoi2.append(N.conjugate(enoi1[hlen-i,:]))
    enoi2 = N.array(enoi2)

    enoi1 = N.delete(enoi1,-1,0)    #delete last row
    enoi = N.concatenate((enoi1,enoi2),axis=0)


    enoit = []
    fti = myfft(dt,nmd)
    for i in range(nx):
        enoit.append(fti.iFourier1D(enoi[:,i])) #w->t
    return N.transpose(N.array(enoit))




#--------------------------------------------------------------------------------------
#driver routine for generating noise
def nonequm(w,bias,T,classical=False):
    """
    w       frequency
    bias    applied bias
    T       equilibrium electronic temperature
    """
    hw1 = U.hbar*w-bias
    hw2 = U.hbar*w
    small = 10e-20
    if classical:
        if hw1 == 0.: hw1=small
        if hw2 == 0.: hw2=small
        return 2.0*hw1*(U.kb*T/hw1-U.kb*T/hw2)
    else:
        return 2.0*hw1*(bose(hw1,T)-bose(hw2,T))

def nonequp(w,bias,T,classical=False):
    """
    w       frequency
    bias    applied bias
    T       equilibrium electronic temperature
    """
    hw1 = U.hbar*w+bias
    hw2 = U.hbar*w
    small = 10e-20
    if classical:
        if hw1 == 0.: hw1=small
        if hw2 == 0.: hw2=small
        return 2.0*hw1*(U.kb*T/hw1-U.kb*T/hw2)
    else:
        return 2.0*hw1*(bose(hw1,T)-bose(hw2,T))

def equ(w,cut,T,classical=False,zpmotion=True):
    """
    w       frequency
    cut     electron band cutoff energy
    T       equilibrium electronic temperature
    """
    small=10e-20
    hw = U.hbar*w
    if zpmotion is True:
        zp=0.5
    else:
        zp=0.0
    if(hw < cut):
        if classical:
            return 2.0*U.kb*T
        else:
            if hw==0:
                return 2.0*U.kb*T
            else:
                return 2.0*hw*(zp+bose(hw,T))
    else:
        return 0.0

def vargau(eval,evec,cof=1.0):
    """
    generate multivariant gaussian:
    eval        a real vector made from the covariance of each dimension
    evec        each COLUMN of evec is one eigen vector 
                NOTE THE DIFFERENCE WITH THE MATHEMATICA CODE
    cof         real coefficient

    All of them should be real.
    Not checked here.
    """
    #print "vargau: checking eigen values and vectors"
    lval = len(eval)
    nval = N.array(eval)
    nvec = N.array(evec)
    svec = N.shape(nvec)
    if(lval != svec[0] or svec[0] != svec[1]):
        print "vargau: shape error"
        sys.exit(0)
    ##if(min(nval) < -10**-10):
    ##    print "vargau: WARNING, found negative eigenvalue"
    ##    print "vargau: minimum negative value: ", min(nval)
    ##    #print nval
    ##    #stop
    sval = [cof*v for v in eval]
    rval = []
    for i in range(len(sval)):
        if(sval[i] <= 0):
            rval.append(0.0)
        else:
            rval.append(N.random.normal(0.0,N.sqrt(sval[i])))
    tmm = mdot(nvec,rval)
    return tmm
