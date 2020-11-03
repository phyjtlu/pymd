import sys

import numpy as N

import units as U


class myfft:
    def __init__(self, dt, n):
        self.dt = dt
        self.N = n
        self.dw = 2*N.pi/dt/n

    def Fourier1D(self, a):
        """ 
        1D Fourier transform from t to w using our definition:
        f(w) = int f(t) e^Iwt dt
        In discrete form, it's
        f(j) = dt sum_{i=0}^{N-1} f(i) e^(I 2pi i j/N).
        N is length of the data. 
        domega/2/pi = 1/N/dt
        This corresponds to the inverse FFT in numpy:
            numpy.fft.ifft(a)*2*pi/domega
        """
        if(len(a) != self.N):
            print("MyFFT.Fourier1D: array length error!")
            sys.exit(0)
        else:
            nor = 2.*N.pi/self.dw
            b = N.fft.ifft(a)
            return nor*b

    def iFourier1D(self, a):
        """ 
        1D Fourier transform from w to t using our definition:
        f(t) = int f(w) e^-Iwt dw/2/pi
        In discrete form, it's
        f(i) = dw/2/pi sum_{j=0}^{N-1} f(j) e^(-I 2pi i j/N).
        N is length of the data. 
        domega/2/pi = 1/N/dt
        This corresponds to the FFT in numpy:
            numpy.fft.fft(a)*domega/2/pi
        """
        if(len(a) != self.N):
            print("MyFFT.iFourier1D: array length error!")
            sys.exit(0)
        else:
            nor = self.dw/2/N.pi
            b = N.fft.fft(a)
            return nor*b


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


def mdot(* args):
    #tmp = args[0].copy()
    #for ii in range(1, len(args)):
    #    tmp = N.dot(tmp, args[ii])
    #return tmp
    return N.linalg.multi_dot([im for im in args])

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

def powerspec(qs, dt, nmd):
    """
    qs      list of trajectories, shape(nmd,nph)
    dt      time step of MD simulation
    nmd     number of MD steps
    """
    qst = N.transpose(N.array(qs))
    nmd2 = qst.shape[1]
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
    nmd2 = pst.shape[1]
    if nmd != nmd2:
        print("power: ps shape error!")
        sys.exit()
    dw = 2.*N.pi/dt/nmd

    fti = myfft(dt, nmd)
    psw = N.array([fti.Fourier1D(a) for a in pst])
    psw = N.real(N.transpose(psw*N.conjugate(psw)))
    dos2 = N.array([[i*dw, N.sum(psw[i])/dt/nmd] for i in range(nmd)])
    return dos2
