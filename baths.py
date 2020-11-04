import sys

import numpy as N

from functions import antisymmetrize, chkShape, flinterp, mdot, symmetrize
from noise import enoise, mf, phnoise


def exlist(a, indices):
    return a[indices]
    #return N.array([a[i] for i in indices])

# fourier transform of gamma (calculated directly, no fft)


def gamt(tl, wl, gwl, gam, eta_ad=0):
    '''
    calculate gam(t)

    tl      equal-spaced time points starting from 0 
    wl      equal-spaced omega points starting from 0
    gam     corresponding friction in omega space


    return the friction kernel in time domain
    '''
    gt = []
    # print "TG test was here!"
    # ea=0.01
    # print "Weird: ", wl[-1]
    # print wl[1]-wl[0],wl[2]-wl[1]
    if eta_ad == 0:
        print("eta=0")
        for t in tl:
            tm = []
            for i in range(len(wl)):
                tm.append(N.array(flinterp(wl[i], gwl, gam))*N.cos(wl[i]*t))
            # 2.0 account for the negative frequency part
            gt.append(2.0*N.mean(N.array(tm), axis=0)*wl[-1]/N.pi)
    else:
        print("eta!=0")
        for t in tl:
            tm = []
            for i in range(len(wl)):
                tm.append(N.array(flinterp(wl[i], gwl, gam))*wl[i]/(wl[i]-1j*eta_ad)*N.exp(-1j*wl[i]*t-eta_ad*t)+N.array(
                    flinterp(wl[i], gwl, gam))*wl[i]/(wl[i]+1j*eta_ad)*N.exp(+1j*wl[i]*t-eta_ad*t))
            gt.append(N.mean(N.array(tm), axis=0)*wl[-1]/N.pi)
            # print "TG imag:", max(N.imag(N.mean(N.array(tm[-1]),axis=0)*wl[-1]/N.pi))
    return N.array(N.real(gt))


class ebath:
    '''
    cats    center degree's of freedom connecting to electrons
    T      equlibrium temperature of electrons
    bias    applied bias
    wmax    cutoff energy of electron bath
    nw      number of energy points
    efric   friction matrix
    exim,exip
    zeta1,zeta2

    dt      md time step
    nmd     md steps

    noise   noise series
    kernel  friction kernel 
    ml      length of kernel(=1)

    ebath   include ebath or not
    '''

    def __init__(self, cats, T, dt, nmd, wmax=None, nw=None, bias=0.,
                 efric=None, exim=None, exip=None, zeta1=None, zeta2=None, classical=False, zpmotion=True):
        self.cats = N.array(cats, dtype='int')
        #self.cids = N.array([[3*c+0,3*c+1,3*c+2] for c in cats]).flatten()
        self.cids = N.array(cats, dtype='int')
        self.nc = len(self.cids)
        self.T, self.wmax = T, wmax
        self.nw, self.bias = nw, bias
        self.dt, self.nmd = dt, nmd
        self.cur = N.zeros(nmd)
        self.classical = classical
        self.zpmotion = zpmotion

        if nw is None or wmax is None:
            self.wl = None
        else:
            self.wl = [self.wmax*i/nw for i in range(nw)]

        self.CheckEmat(efric, exim, exip, zeta1, zeta2)

        # local friction
        self.ml = 1
        self.noise = None

    def CheckEmat(self, efric=None, exim=None, exip=None, zeta1=None, zeta2=None):
        '''
        check the matrix, and set the following variables:
            efric
            exim
            exip
            zeta1
            zeta2
            ebath
        '''
        # friction
        if efric is not None:
            print("ebath.CheckEmat: got efric, checking")
            n = chkShape(efric)
            if n != self.nc:
                print("ebath.CheckEmat: efric shape error!")
                sys.exit()
            print("ebath.setEmat: symmetrizing efric")
            self.efric = symmetrize(efric)
            self.kernel = N.array([self.efric])
            self.exip = N.zeros(shape=(n, n))
            self.exim = N.zeros(shape=(n, n))
            self.zeta1 = N.zeros(shape=(n, n))
            self.zeta2 = N.zeros(shape=(n, n))
            self.ebath = True
        else:
            print("ebath.CheckEmat: no efric provided, setting ebath to False")
            self.efric = None
            self.kernel = None
            self.exim = None
            self.exip = None
            self.zeta1 = None
            self.zeta2 = None
            self.ebath = False
            return

        # nc
        if exim is not None:
            print("ebath.CheckEmat: got exim, checking")
            n = chkShape(exim)
            if self.nc != n:
                print("ebath.CheckEmat: the dimension of exim is wrong!")
                sys.exit(0)
            print("ebath.setEmat: antisymmetrizing exim")
            self.exim = antisymmetrize(exim)

        if exip is not None:
            print("ebath.CheckEmat: got exip, checking")
            n = chkShape(exip)
            if self.nc != n:
                print("ebath.setEmat: the dimension of exip is wrong!")
                sys.exit(0)
            print("ebath.setEmat: symmetrizing exip")
            self.exip = symmetrize(exip)

        # renormalization
        if zeta1 is not None:
            print("ebath.CheckEmat: got zeta1, checking")
            n = chkShape(zeta1)
            if self.nc != n:
                print("ebath.setEmat: the dimension of zeta1 is wrong!")
                sys.exit(0)
            print("ebath.setEmat: symmetrizing zeta1")
            self.zeta1 = symmetrize(zeta1)

        # berry
        if zeta2 is not None:
            print("ebath.CheckEmat: got zeta2, checking")
            n = chkShape(zeta2)
            if self.nc != n:
                print("ebath.setEmat: the dimension of zeta2 is wrong!")
                sys.exit(0)
            print("ebath.setEmat: antisymmetrizing zeta2")
            self.zeta2 = antisymmetrize(zeta2)
        return

    def gnoi(self):
        '''
        electronic noise 
        '''
        if self.nmd is None:
            print("ebath.gnoi: nmd not set!")
            sys.exit()
        if self.dt is None:
            print("ebath.gnoi: dt not set!")
            sys.exit()
        if self.ebath is False:
            print("ebath.gnoi: ebath is False!")
            sys.exit()
        print("ebath.gnoi:classical: %r" % self.classical)
        print("ebath.gnoi:including zero point motion: %r" % self.zpmotion)
        self.noise = N.real(enoise(self.efric, self.exim, self.exip,
                                   self.bias, self.T, self.wmax, self.dt, self.nmd, self.classical, self.zpmotion))

    def GetSig(self):
        '''
        effective retarded self-energy in the wideband limit
        '''
        if self.wl is None:
            print("ebath.GetSig:wl is not set")
            sys.exit()
        else:
            wl = self.wl
        nw = len(wl)
        nc = chkShape(self.efric)
        self.sig = N.zeros((nw, nc, nc), N.complex)

        for i in range(nw):
            self.sig[i] = -1.j*wl[i]*(self.efric+self.bias*self.zeta2)\
                + self.bias*self.zeta1-self.bias*self.exim

    def SetMDsteps(self, dt, nmd):
        self.dt, self.nmd = dt, nmd
        print("ebath.SetMDsteps: memory len reset, you need to \
                regenerate the noise")

    def setbias(self, bias=0.0):
        '''
        set bias applied to the system
        '''
        self.bias = bias
        print("ebath.setbias: bias set to:%f" % self.bias)
        print("ebath.setbias: WARNING--BIAS CHANGED! YOU NEED TO REGENERATE THE NOISE!")

    def bforce(self, t, phis, qhis):
        '''
        return the force from baths, including noise and friction

        note that some of them are bias dependent

        BIAS IS DEFINED AS MUL-MUR
        '''
        f = self.noise[t % self.nmd]  # noise
        if not (self.exim.any() and self.zeta1.any() and self.zeta2.any()):
            for i in range(self.ml):  # friction,nc,rn,berry
                if self.ml == 1:
                    f = f-mdot(self.kernel[i], exlist(phis[i], self.cids))
                else:
                    print("WARNING: nonlocal electronic force not implemented!")
                    # stophere
                    f = f-mdot(self.kernel[i], exlist(phis[i], self.cids))*self.dt
        else:
            for i in range(self.ml):  # friction,nc,rn,berry
                if self.ml == 1:
                    f = f-mdot(self.kernel[i], exlist(phis[i], self.cids))\
                        + mdot(self.bias*self.exim, exlist(qhis[0], self.cids))\
                        - mdot(self.bias*self.zeta1, exlist(qhis[0], self.cids))\
                        - mdot(self.bias*self.zeta2, exlist(phis[0], self.cids))
                else:
                    print("WARNING: nonlocal electronic force not implemented!")
                    # stophere
                    f = f-mdot(self.kernel[i], exlist(phis[i], self.cids))*self.dt
        return mf(f, self.cids, len(phis[0]))


class phbath:
    '''
    class for the phonon baths

    ##########################################################
    INPUT
    ##########################################################
    T           temperature of the bath
    cats        central system atoms connecting to the bath (python index)
    debye       Debye frequency of the lead material
    nw          number of sampling points in positive w space
    dt,nmd      time step and number of steps of md simulation
    ml          length of memory kernel
    mcof        the maximum frequency in wl is mcof*debye

    sig         self-energy matrix
    gamma       friction kernel in w space
    gwl         energy points of self-energy matrix

    K00,K01     the dynamical matrix of the semi-infinite leads
    V01         interaction matrix of lead with central system
    ##########################################################


    ##########################################################
    GENERATED VARIABLES
    ##########################################################
    nc          how many degrees of freedom are connected to central system
    local       if the friction is time-local or not
    wmax        cutoff frequency
    wl          list of frequencies, starting from 0
    cur         save the energy current to the center at each md time step
    noise       noise series
    kernel      friction kernel with length ml
    '''

    def __init__(self, T, cats, debye, nw, dt, nmd, ml=None, mcof=2.0, sig=None, gamma=None, gwl=None, K00=None, K01=None, V01=None, eta_ad=0, classical=False, zpmotion=True):
        self.classical = classical
        self.zpmotion = zpmotion
        self.T, self.debye, self.cats = T, debye, N.array(cats, dtype='int')
        #self.K00,self.K01,self.V01 = N.array(K00),N.array(K01),N.array(V01)
        self.K00, self.K01, self.V01 = K00, K01, V01
        self.dt, self.nmd, self.ml = dt, nmd, ml
        self.kernel = None
        self.cids = N.array([[3*c+0, 3*c+1, 3*c+2] for c in cats]).flatten()
        self.nc = len(self.cids)
        self.wmax = mcof*debye
        self.local = False
        self.nw = nw
        self.wl = [self.wmax*i/nw for i in range(nw)]
        self.gamma = gamma
        self.sig = sig
        self.gwl = gwl
        self.cur = N.zeros(nmd)
        self.eta_ad = eta_ad
        #
        # initialise self.gamma
        if self.UseK():
            # calculate the real self-energy
            self.local = False
            print("phbath: Calculating self-energy is not implemented yet.")
            sys.exit(0)

        elif self.UseG() or self.UsePi():
            if self.UsePi():
                if len(self.sig[0]) != self.nc:
                    print("phbath: inconsist between cids and sig!")
                    sys.exit()
                self.ggamma()
            if self.UseG():
                if len(self.gamma[0]) != self.nc:
                    print("phbath: inconsist between cids and gamma!")
                    sys.exit()
        else:
            # use Debye model
            # friction coefficient from Debye frequency
            # see, Adelman&Doll, JCP, Vol.64,2375 (1976)
            phfric = debye*N.pi/6.0
            self.gamma = N.array([N.diag(phfric+N.zeros(int(self.nc)))])
            self.gwl = N.array([0])
            self.local = True
            self.ml = 1

    def SetMDsteps(self, dt, nmd):
        self.dt, self.nmd = dt, nmd
        print("phbath.SetMDsteps: memory len reset, you need to \
                regenerate the memory kernel and noise")

    def SetMemlen(self, len):
        self.ml = len
        print("phbath.SetMemlen: memory len reset, you need to \
                regenerate the memory kernel")

    def SetT(self, T):
        print("phbath.SetT: Set temperature  to %s\n" % T)
        self.T = T

    def UseG(self):
        if self.gamma is not None and self.gwl is not None:
            return True
        else:
            return False

    def UsePi(self):
        if self.sig is not None and self.gwl is not None:
            return True
        else:
            return False

    def UseK(self):
        if self.K00 is not None and self.K01 is not None \
           and self.V01 is not None:
            return True
        else:
            return False

    def ggamma(self):
        '''
        calculate gamma from self-energy
        Gamma(w) = -Im(Sig(w))/w
        '''
        if self.sig is not None:
            print("2TG test was here2!")
            Sig = self.sig
            wl = self.gwl
            a = []
            for i in range(len(wl)):
                if wl[i] == 0:
                    # a.append(-N.imag(Sig[0]))#
                    a.append(-N.imag(Sig[i+1])/wl[i+1])
                else:
                    # a.append(-N.imag(Sig[0]))#
                    a.append(-N.imag(Sig[i])/wl[i])
            self.gamma = N.array(a)
        else:
            print("phbath.Gamma: self.sig is not set, need it to calculate gamma")
            sys.exit()

    def gnoi(self):
        '''
        generate phonon noise using the friction kernel
        '''
        if self.dt is None or self.nmd is None:
            print("phbath.gnoi: the md information dt and nmd are not set!")
            sys.exit()
        print("regenerating noise")
        print("WARNING: remember to reset t=0 to use the new noise!")
        print("phbath.gnoi:classical: %r" % self.classical)
        print("phbath.gnoi:including zero point motion:%r" % self.zpmotion)
        self.noise = N.real(phnoise(self.gamma, self.gwl, self.T,
                                    self.wmax, self.dt, self.nmd, self.classical, self.zpmotion))
        # N.linspace(self.gwl[0],self.gwl[-1])

    def gmem(self):
        '''
        generate the memory kernel in time domain
        '''
        if self.ml is None or self.dt is None:
            print("phbath.gmem: length of memory kernel not set!")
            sys.exit()
        if self.local:
            self.ml = 1
            self.kernel = self.gamma
        else:
            tl = [self.dt*i for i in range(self.ml)]
            self.kernel = N.real(
                gamt(tl, self.wl, self.gwl, self.gamma, self.eta_ad))
            # update gamma to include artificial damping:
            if self.eta_ad != 0:
                # print "TG test FFT was here"
                #gamnew = N.zeros(self.kernel.shape)
                #fti = myfft(self.dt,self.kernel.shape[0])
                # for i in range(self.kernel.shape[1]):
                #    for j in range(self.kernel.shape[2]):
                #        gamnew[:,i,j]=(fti.Fourier1D(self.kernel[:,i,j])) #t->w, but only positive t!
                # self.gamma=N.array(2*N.real(gamnew)) # Real since its really a cos-transform (kernel(-t)=kernel(t)).
                gamnew = N.zeros(self.gamma.shape)
                for i in range(len(self.gwl)):
                    for it in range(self.kernel.shape[0]):
                        gamnew[i, :, :] = gamnew[i, :, :]+self.dt * \
                            self.kernel[it, :, :]*N.cos(self.gwl[i]*tl[it])
                # self.gamma=N.array(len(self.gwl)/self.kernel.shape[0]*N.real(gamnew)) # Real since its really a cos-transform (kernel(-t)=kernel(t)).
                self.gammaOld = self.gamma
                # Real since its really a cos-transform (kernel(-t)=kernel(t)).2*self.gwl[i]*
                self.gamma = N.array(N.real(gamnew))
                # print "Test igen:",self.gamma.shape

    def bforce(self, t, phis, qhis):
        '''
        return the force from baths, including noise and friction
        '''
        f = self.noise[t % self.nmd]  # noise
        for i in range(self.ml):  # friction
            if self.ml == 1:
                f = f-mdot(self.kernel[i], exlist(phis[i], self.cids))
            else:
                f = f-mdot(self.kernel[i], exlist(phis[i], self.cids))*self.dt
        return mf(f, self.cids, len(phis[0]))
