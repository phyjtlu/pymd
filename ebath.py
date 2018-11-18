#!/applications/mbrsoft/bin/python

import sys
import numpy as N
from noise import *
from functions import *
from netCDF4 import Dataset


class ebath:
    """
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
    """
    def __init__(self,cats,T,dt,nmd,wmax=None,nw=None,bias=0.,\
            efric=None,exim=None,exip=None,zeta1=None,zeta2=None,classical=False,zpmotion=True):
        self.cats=N.array(cats,dtype='int')
        self.cids = N.array([[3*c+0,3*c+1,3*c+2] for c in cats]).flatten()
        self.nc = len(self.cids)
        self.T,self.wmax = T,wmax
        self.nw,self.bias = nw,bias
        self.dt,self.nmd = dt,nmd
        self.cur = N.zeros(nmd)
        self.classical=classical
        self.zpmotion=zpmotion

        if nw is None or wmax is None:
            self.wl=None
        else:
            self.wl = [self.wmax*i/nw for i in range(nw)]

        self.CheckEmat(efric,exim,exip,zeta1,zeta2)

        #local friction
        self.ml=1
        self.noise=None




    def CheckEmat(self,efric=None,exim=None,exip=None,zeta1=None,zeta2=None):
        """
        check the matrix, and set the following variables:
            efric
            exim
            exip
            zeta1
            zeta2
            ebath
        """
        #friction
        if efric is not None:
            print "ebath.CheckEmat: got efric, checking"
            n = chkShape(efric)
            if n != self.nc:
                print "ebath.CheckEmat: efric shape error!"
                sys.exit()
            print "ebath.setEmat: symmetrizing efric"
            self.efric = symmetrize(efric)
            self.kernel = N.array([self.efric])
            self.exip = N.zeros(shape=(n,n))
            self.exim = N.zeros(shape=(n,n))
            self.zeta1=N.zeros(shape=(n,n))
            self.zeta2=N.zeros(shape=(n,n))
            self.ebath = True
        else:
            print "ebath.CheckEmat: no efric provided, setting ebath to False"
            self.efric = None
            self.kernel = None
            self.exim = None
            self.exip = None
            self.zeta1= None
            self.zeta2= None
            self.ebath = False
            return

        #nc
        if exim is not None:
            print "ebath.CheckEmat: got exim, checking"
            n = chkShape(exim)
            if self.nc != n:
                print "ebath.CheckEmat: the dimension of exim is wrong!"
                sys.exit(0)
            print "ebath.setEmat: antisymmetrizing exim"
            self.exim = antisymmetrize(exim)

        if exip is not None:
            print "ebath.CheckEmat: got exip, checking"
            n = chkShape(exip)
            if self.nc != n:
                print "ebath.setEmat: the dimension of exip is wrong!"
                sys.exit(0)
            print "ebath.setEmat: symmetrizing exip"
            self.exip = symmetrize(exip)

        #renormalization
        if zeta1 is not None:
            print "ebath.CheckEmat: got zeta1, checking"
            n = chkShape(zeta1)
            if self.nc != n:
                print "ebath.setEmat: the dimension of zeta1 is wrong!"
                sys.exit(0)
            print "ebath.setEmat: symmetrizing zeta1"
            self.zeta1 = symmetrize(zeta1)

        #berry
        if zeta2 is not None:
            print "ebath.CheckEmat: got zeta2, checking"
            n = chkShape(zeta2)
            if self.nc != n:
                print "ebath.setEmat: the dimension of zeta2 is wrong!"
                sys.exit(0)
            print "ebath.setEmat: antisymmetrizing zeta2"
            self.zeta2 = antisymmetrize(zeta2)
        return

    def gnoi(self):
        """
        electronic noise 
        """
        if self.nmd is None:
            print "ebath.gnoi: nmd not set!"
            sys.exit()
        if self.dt is None:
            print "ebath.gnoi: dt not set!"
            sys.exit()
        if self.ebath is False:
            print "ebath.gnoi: ebath is False!"
            sys.exit()
        print "ebath.gnoi:classical:", self.classical
        print "ebath.gnoi:including zero point motion:", self.zpmotion
        self.noise=N.real(enoise(self.efric,self.exim,self.exip,\
                  self.bias,self.T,self.wmax,self.dt,self.nmd,self.classical,self.zpmotion))

    def GetSig(self):
        """
        effective retarded self-energy in the wideband limit
        """
        if self.wl is None:
            print "ebath.GetSig:wl is not set"
            sys.exit()
        else:
            wl=self.wl
        nw=len(wl)
        nc=chkShape(self.efric)
        self.sig=N.zeros((nw,nc,nc),N.complex)
    
        for i in range(nw):
            self.sig[i]=-1.j*wl[i]*(self.efric+self.bias*self.zeta2)\
                    +self.bias*self.zeta1-self.bias*self.exim

    def SetMDsteps(self,dt,nmd):
        self.dt,self.nmd=dt,nmd
        print "ebath.SetMDsteps: memory len reset, you need to \
                regenerate the noise"

    def setbias(self,bias=0.0):
        """
        set bias applied to the system
        """
        self.bias = bias
        print "ebath.setbias: bias set to:", self.bias
        print "ebath.setbias: WARNING--BIAS CHANGED! YOU NEED TO REGENERATE THE NOISE!"

    def bforce(self,t,phis,qhis):
        """
        return the force from baths, including noise and friction
        
        note that some of them are bias dependent

        BIAS IS DEFINED AS MUL-MUR
        """
        f = self.noise[t%self.nmd] #noise
        for i in range(self.ml): #friction,nc,rn,berry
            if self.ml == 1:
                f=f-mm(self.kernel[i],exlist(phis[i],self.cids))\
                        +mm(self.bias*self.exim,exlist(qhis[0],self.cids))\
                        -mm(self.bias*self.zeta1,exlist(qhis[0],self.cids))\
                        -mm(self.bias*self.zeta2,exlist(phis[0],self.cids))
            else:
                print "WARNING: nonlocal electronic force not implemented!"
                stophere
                f=f-mm(self.kernel[i],exlist(phis[i],self.cids))*self.dt
        return mf(f,self.cids,len(phis[0]))


#--------------------------------------------------------------------------------------
#driver
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
#--------------------------------------------------------------------------------------
#testing
#
#--------------------------------------------------------------------------------------

if __name__=="__main__":
    import matplotlib.pyplot as PP
#--------------------------------------------------------------------------------------
    def ReadEPHNCFile(filename):
        """
        Reads a NetCDF file that describes dynamical matrix, self-energies
        """
        class eph:
            pass
    
        #file = nc.NetCDFFile(filename,'r')
        file = Dataset(filename,'r')
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
    
    
    
    filename="EPH.nc"
    
    eph=ReadEPHNCFile(filename)
    
    #--------------------------------------------------------------------------------------
    #----------------------------------------------------
    # electron bath
    #    def __init__(self,cats,T=0.,wmax=None,nw=None,bias=0.,efric=None,exim=None,exip=None,dt=None,nmd=None):
    #----------------------------------------------------
    ecats=range(13)
    eb = ebath(ecats,T=300.0,wmax=1.0,nw=500,bias=1.,efric=eph.efric,exim=eph.xim,exip=eph.xip)
    eb.SetMDsteps(8,2**12)
    eb.genoi()
    PP.plot(eb.noise[:,0])
    PP.savefig("enoise.pdf")
    PP.close()
    sys.exit()
