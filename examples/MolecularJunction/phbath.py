#!/applications/mbrsoft/bin/python

import sys
import numpy as N
from noise import *
from functions import *


class phbath:
    """
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
    """
    def __init__(self,T,cats,debye,nw,dt,nmd,ml=None,mcof=2.0,sig=None,gamma=None,gwl=None,K00=None,K01=None,V01=None,eta_ad=0,classical=False,zpmotion=True):
        self.classical=classical
        self.zpmotion=zpmotion
        self.T,self.debye,self.cats = T,debye,N.array(cats,dtype='int')
        #self.K00,self.K01,self.V01 = N.array(K00),N.array(K01),N.array(V01)
        self.K00,self.K01,self.V01 = K00,K01,V01
        self.dt,self.nmd,self.ml = dt,nmd,ml
        self.kernel = None
        self.cids = N.array([[3*c+0,3*c+1,3*c+2] for c in cats]).flatten()
        self.nc = len(self.cids)
        self.wmax = mcof*debye
        self.local = False
        self.nw = nw
        self.wl = [self.wmax*i/nw for i in range(nw)]
        self.gamma = gamma
        self.sig = sig
        self.gwl = gwl
        self.cur = N.zeros(nmd)
        self.eta_ad=eta_ad
        #
        #initialise self.gamma
        if self.UseK():
            #calculate the real self-energy
            self.local = False
            print "phbath: Calculating self-energy is not implemented yet."
            sys.exit(0)

        elif self.UseG() or self.UsePi():
            if self.UsePi():
                if len(self.sig[0]) != self.nc:
                    print "phbath: inconsist between cids and sig!"
                    sys.exit()
                self.ggamma()
            if self.UseG():
                if len(self.gamma[0]) != self.nc:
                    print "phbath: inconsist between cids and gamma!"
                    sys.exit()
        else:
            #use Debye model
            #friction coefficient from Debye frequency
            #see, Adelman&Doll, JCP, Vol.64,2375 (1976)
            phfric = debye*N.pi/6.0
            self.gamma = N.array([N.diag(phfric+N.zeros(int(self.nc)))])
            self.gwl = N.array([0])
            self.local = True
            self.ml = 1


    def SetMDsteps(self,dt,nmd):
        self.dt,self.nmd=dt,nmd
        print "phbath.SetMDsteps: memory len reset, you need to \
                regenerate the memory kernel and noise"

    def SetMemlen(self,len):
        self.ml = len
        print "phbath.SetMemlen: memory len reset, you need to \
                regenerate the memory kernel"

    def SetT(self,T):
        print "phbath.SetT: Set temperature  to %s\n"%T
        self.T=T


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
        """
        calculate gamma from self-energy
        Gamma(w) = -Im(Sig(w))/w
        """
        if self.sig is not None:
            print "2TG test was here2!"
            Sig=self.sig
            wl=self.gwl
            a=[]
            for i in range(len(wl)):
                if wl[i] == 0:
                    a.append(-N.imag(Sig[i+1])/wl[i+1]) #a.append(-N.imag(Sig[0]))#
                else:
                    a.append(-N.imag(Sig[i])/wl[i]) #a.append(-N.imag(Sig[0]))#
            self.gamma=N.array(a)
        else:
            print "phbath.Gamma: self.sig is not set, need it to calculate gamma"
            sys.exit()

    def gnoi(self):
        """
        generate phonon noise using the friction kernel
        """
        if self.dt is None or self.nmd is None:
            print "phbath.gnoi: the md information dt and nmd are not set!"
            sys.exit()
        print "regenerating noise"
        print "WARNING: remember to reset t=0 to use the new noise!"
        print "phbath.gnoi:classical:", self.classical
        print "phbath.gnoi:including zero point motion:", self.zpmotion
        self.noise=N.real(phnoise(self.gamma,self.gwl,self.T,self.wmax,self.dt,self.nmd,self.classical,self.zpmotion))
        #N.linspace(self.gwl[0],self.gwl[-1])


    def gmem(self):
        """
        generate the memory kernel in time domain
        """
        if self.ml is None or self.dt is None:
            print "phbath.gmem: length of memory kernel not set!"
            sys.exit()
        if self.local:
            self.ml = 1
            self.kernel = self.gamma
        else:
            tl = [self.dt*i for i in range(self.ml)]
            self.kernel = N.real(gamt(tl,self.wl,self.gwl,self.gamma,self.eta_ad))
            # update gamma to include artificial damping:
            if self.eta_ad != 0:
               #print "TG test FFT was here"
               #gamnew = N.zeros(self.kernel.shape)
               #fti = myfft(self.dt,self.kernel.shape[0])
               #for i in range(self.kernel.shape[1]):
               #    for j in range(self.kernel.shape[2]):
               #        gamnew[:,i,j]=(fti.Fourier1D(self.kernel[:,i,j])) #t->w, but only positive t!
               #self.gamma=N.array(2*N.real(gamnew)) # Real since its really a cos-transform (kernel(-t)=kernel(t)).
               gamnew = N.zeros(self.gamma.shape)
               for i in range(len(self.gwl)):
                    for it in range(self.kernel.shape[0]):
                        gamnew[i,:,:]=gamnew[i,:,:]+self.dt*self.kernel[it,:,:]*N.cos(self.gwl[i]*tl[it])
               #self.gamma=N.array(len(self.gwl)/self.kernel.shape[0]*N.real(gamnew)) # Real since its really a cos-transform (kernel(-t)=kernel(t)).
               self.gammaOld=self.gamma
               self.gamma=N.array(N.real(gamnew)) # Real since its really a cos-transform (kernel(-t)=kernel(t)).2*self.gwl[i]*
               #print "Test igen:",self.gamma.shape

    def bforce(self,t,phis,qhis):
        """
        return the force from baths, including noise and friction
        """
        f = self.noise[t%self.nmd] #noise
        for i in range(self.ml): #friction
            if self.ml == 1:
                f=f-mm(self.kernel[i],exlist(phis[i],self.cids))
            else:
                f=f-mm(self.kernel[i],exlist(phis[i],self.cids))*self.dt
        return mf(f,self.cids,len(phis[0]))



#--------------------------------------------------------------------------------------
def mf(f,cats,lens):
    """
    padding f to dimension len
    """
    t = N.zeros(lens)
    for i in range(len(cats)):
        t[cats[i]]=f[i]
    return t


def exlist(a,indices):
    return N.array([a[i] for i in indices])

#fourier transform of gamma (calculated directly, no fft)
def gamt(tl,wl,gwl,gam,eta_ad=0):
    """
    calculate gam(t) = \theta(t) \int_{-\infty}^{+\infty} -Im{\Pi^r[w]}/w exp(-iwt)dw/pi 
                     =2\theta(t) \int_{0}^{+\infty} -Im{\Pi^r[w]}/w exp(-iwt)dw/pi 

    tl      equal-spaced time points starting from 0 
    wl      equal-spaced omega points starting from 0
    gam     corresponding friction in omega space


    return the friction kernel in time domain
    """
    gt=[]
    #print "TG test was here!"
    #ea=0.01
    #print "Weird: ", wl[-1]
    #print wl[1]-wl[0],wl[2]-wl[1]
    if eta_ad==0:
       print "eta=0"
       for t in tl:
           tm=[]
           for i in range(len(wl)):
               tm.append(N.array(flinterp(wl[i],gwl,gam))*N.cos(wl[i]*t))
           #2.0 account for the negative frequency part
           gt.append(2.0*N.mean(N.array(tm),axis=0)*wl[-1]/N.pi)
    else:
       print "eta!=0"
       for t in tl:
           tm=[]
           for i in range(len(wl)):
               tm.append(N.array(flinterp(wl[i],gwl,gam))*wl[i]/(wl[i]-1j*eta_ad)*N.exp(-1j*wl[i]*t-eta_ad*t)+N.array(flinterp(wl[i],gwl,gam))*wl[i]/(wl[i]+1j*eta_ad)*N.exp(+1j*wl[i]*t-eta_ad*t))
           gt.append(N.mean(N.array(tm),axis=0)*wl[-1]/N.pi)
           #print "TG imag:", max(N.imag(N.mean(N.array(tm[-1]),axis=0)*wl[-1]/N.pi))
    return N.array(N.real(gt))


#--------------------------------------------------------------------------------------
#testing
#
#--------------------------------------------------------------------------------------

if __name__=="__main__":
    #import matplotlib.pyplot as PP
#--------------------------------------------------------------------------------------
    print "testing ebath - begin"
    #set baths
    a = N.diag(1.+N.zeros(3))
    test = phbath(0.,[1],1.0,100,dt=None,nmd=None,ml=None,mcof=2.0,gamma=None,gwl=None,K00=None,K01=None,V01=None)
    print phbath.__doc__


    test.SetMDsteps(0.5,2**14)
    print "generating phonon noise"
    test.gnoi()
    test.gmem()
    print test.noise.shape
    print N.mean(test.noise,0)
    print test.kernel
    #PP.plot(test.noise[:,0])
    #PP.show()
    #PP.plot(test.kernel[:,0])
    #PP.show()
