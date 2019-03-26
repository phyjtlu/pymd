#!/applications/mbrsoft/bin/python

import sys
import numpy as N
import units as U
from myfft import *

def powerspec(qs,dt,nmd):
    """
    qs      list of trajectories, shape(nmd,nph)
    dt      time step of MD simulation
    nmd     number of MD steps
    """
    qst=N.transpose(N.array(qs))
    nph,nmd2 = qst.shape
    if nmd != nmd2:
        print "power: qs shape error!"
        sys.exit()
    dw = 2.*N.pi/dt/nmd

    fti = myfft(dt,nmd)
    qsw = N.array([fti.Fourier1D(a) for a in qst])
    qsw = N.real(N.transpose(qsw*N.conjugate(qsw)))
    dos = N.array([[i*dw,(dw*i)**2*N.sum(qsw[i])/dt/nmd] for i in range(nmd)])
    return dos


def powerspec2(ps,dt,nmd):
    """
    ps      list of trajectories, shape(nmd,nph)
    dt      time step of MD simulation
    nmd     number of MD steps
    """
    pst=N.transpose(N.array(ps))
    nph,nmd2 = pst.shape
    if nmd != nmd2:
        print "power: ps shape error!"
        sys.exit()
    dw = 2.*N.pi/dt/nmd

    fti = myfft(dt,nmd)
    psw = N.array([fti.Fourier1D(a) for a in pst])
    psw = N.real(N.transpose(psw*N.conjugate(psw)))
    dos2 = N.array([[i*dw,N.sum(psw[i])/dt/nmd] for i in range(nmd)])
    return dos2




