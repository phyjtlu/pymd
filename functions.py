import glob
import sys

import numpy as N

import units as U
from myfft import myfft

N.seterr(over="ignore")


def coth(x):
    """
    hyperbolic cotangent function
    """
    if x == 0.0:
        print("coth:coth(0) is infinity")
        sys.exit(0)
    else:
        return N.cosh(x)/N.sinh(x)


def xcoth(x):
    """
    x*coth(x)
    """
    if x == 0.0:
        return 1.0
    else:
        return x*N.cosh(x)/N.sinh(x)


def bose(w, T):
    """
    bose distribution
    """
    #small = 10e-20
    if T == 0.0:
        if w == 0.0:
            return 1/(N.exp(1.0/U.kb)-1)
        elif w < 0.0:
            return -1.0
        else:
            return 0.0
    else:
        if w == 0.0:
            # return 1/small
            # have problems for finite temperature for bias calculation
            # return 0 seems solves it
            return 0.0
        else:
            return 1.0/(N.exp(w/U.kb/T)-1.0)


def fermi(ep, mu, T):
    """
    fermi distribution
    """
    if T == 0.0:
        if(ep < mu):
            return 1.0
        elif(ep > mu):
            return 0.0
        else:
            return 0.5
    else:
        return 1/(N.exp((ep-mu)/U.kb/T)+1)


def flinterp(x, xs, ys):
    """
    do a linear interpolation of (xs,ys), return the interpolated value at x.
    """
    id = nearest(x, xs)

    # boundaries
    if id == len(xs)-1:
        return ys[-1]
    if id == 0:
        return ys[0]

    # linear interpolation
    dd = xs[id]-x
    if dd < 0:
        return ys[id]+dd/(xs[id]-xs[id-1])*(ys[id]-ys[id-1])
    else:
        return ys[id]+dd/(xs[id]-xs[id+1])*(ys[id]-ys[id+1])


def nearest(b, bs):
    """
    return the index of the element in bs wich is the nearest to b.
    """
    bsn = N.array(bs)
    bst = abs(bsn-b)
    return(list(bst).index(min(bst)))


def rpadleft(bs, b):
    if len(bs) > 1:
        return N.concatenate((N.array([b]), N.array(bs)[:-1]), axis=0)
    elif len(bs) == 1:
        return N.array([b])
    else:
        print("len(bs) is less than 1")
        sys.exit()


def mm(* args):
    tmp = args[0].copy()
    for ii in range(1, len(args)):
        tmp = N.dot(tmp, args[ii])
    return tmp


def chkShape(a):
    """
    check if a is a n by n matrix, if yes return n
    """
    aa = N.array(a)
    ash = N.shape(aa)
    if(ash[0] == ash[1]):
        return ash[0]
    else:
        print("The matrix should be a n by n matrix")
        sys.exit(0)


def symmetrize(a):
    aa = N.array(a)
    return 0.5*(aa+N.transpose(aa))


def antisymmetrize(a):
    aa = N.array(a)
    return 0.5*(aa-N.transpose(aa))


def dagger(a):
    aa = N.array(a)
    ash = N.shape(aa)
    if(ash[0] != ash[1]):
        print("Not sqaure matrix")
        sys.exit(0)
    return N.transpose(N.conjugate(aa))


def hermitianize(a):
    aa = N.array(a)
    return 0.5*(aa+dagger(aa))


def mdot(* args):
    # dot product with arbitrary number of arguments
    tmp = N.identity(len(args[0]))
    for ii in range(len(args)):
        tmp = N.dot(tmp, args[ii])
    return tmp


def powerspec(qs, dt, nmd):
    """
    qs      list of trajectories, shape(nmd,nph)
    dt      time step of MD simulation
    nmd     number of MD steps
    """
    qst = N.transpose(N.array(qs))
    nph, nmd2 = qst.shape
    if nmd != nmd2:
        print("power: qs shape error!")
        sys.exit()
    dw = 2.*N.pi/dt/nmd

    fti = myfft(dt, nmd)
    qsw = N.array([fti.Fourier1D(a) for a in qst])
    qsw = N.real(N.transpose(qsw*N.conjugate(qsw)))
    dos = N.array([[i*dw, (dw*i)**2*N.sum(qsw[i])/dt/nmd] for i in range(nmd)])
    return dos


def powerspec2(ps, dt, nmd):
    """
    ps      list of trajectories, shape(nmd,nph)
    dt      time step of MD simulation
    nmd     number of MD steps
    """
    pst = N.transpose(N.array(ps))
    nph, nmd2 = pst.shape
    if nmd != nmd2:
        print("power: ps shape error!")
        sys.exit()
    dw = 2.*N.pi/dt/nmd

    fti = myfft(dt, nmd)
    psw = N.array([fti.Fourier1D(a) for a in pst])
    psw = N.real(N.transpose(psw*N.conjugate(psw)))
    dos2 = N.array([[i*dw, N.sum(psw[i])/dt/nmd] for i in range(nmd)])
    return dos2


def calHF(dlist=1):
    # calculate average heat flux
    print("Calculate heat flux.")
    # temperture=temp
    for filename in glob.glob('./kappa.*.bath0.run0.dat'):
        with open(filename, 'r') as f:
            for line in f:
                temperture = float(line.split()[1])

    dlist = list(range(dlist))
    times = int(len(glob.glob('./kappa.*.bath*.run*.dat'))/2)
    kb = N.empty([2, times])

    for i in range(2):
        for j in range(times):
            kappafile = "./kappa." + \
                str(int(temperture))+".bath"+str(i)+".run"+str(j)+".dat"
            for files in glob.glob(kappafile):
                with open(files, 'r') as f:
                    for line in f:
                        kb[i][j] = line.split()[2]
#                        temperture=float(line.split()[1])
    oldkb = N.delete(kb, dlist, axis=1)
    balancekb = N.delete(kb, dlist, axis=1)
    for i in range(balancekb.shape[0]):
        for j in range(balancekb.shape[1]):
            balancekb[i][j] = N.mean(oldkb[i][0:j+1])

    heatflux = (balancekb[0]-balancekb[1])/2

    with open('heatflux.'+str(int(temperture))+'.dat', 'w') as f:
        f.write("Temperture\t"+str(temperture)+"\n"+"Bath0\t" +
                str(balancekb[0])+"\n"+"Bath1\t"+str(balancekb[1])+"\n"+"HeatFlux\t"+str(heatflux)+"\n"+"\n")


def calTC(delta, dlist=0):
    # calculate thermal conductance
    print("Calculate thermal conductance.")
    delta = delta
    # temperture=temp
    for filename in glob.glob('./kappa.*.bath0.run0.dat'):
        with open(filename, 'r') as f:
            for line in f:
                temperture = float(line.split()[1])
    dlist = list(range(dlist))
    times = int(len(glob.glob('./kappa.*.bath*.run*.dat'))/2)
    kb = N.empty([2, times])

    for i in range(2):
        for j in range(times):
            kappafile = "./kappa." + \
                str(int(temperture))+".bath"+str(i)+".run"+str(j)+".dat"
            for files in glob.glob(kappafile):
                with open(files, 'r') as f:
                    for line in f:
                        kb[i][j] = line.split()[2]
#                        temperture=float(line.split()[1])
    kappa = (kb[0]-kb[1])/2/(delta*temperture)
    kappa = N.delete(kappa, dlist)
    # for i in range(len(kappa)):
    #    kappa[i]=N.mean(kappa[0:i+1])

    with open('thermalconductance.'+str(int(temperture))+'.dat', 'w') as f:
        f.write("Temperture\t"+str(temperture)+"\n"+"ThermalConductance\t"+str(kappa)+"\n" +
                "Mean\t"+str(N.mean(kappa))+"\n"+"StandardDeviation\t"+str(N.std(kappa))+"\n")
