#import matplotlib.pyplot as PP
import Scientific.IO.NetCDF as nc
from md import *
from phbath import *

def ReadEPHNCFile(filename):
    """
    Reads a NetCDF file that describes dynamical matrix, self-energies
    """
    class eph:
        pass

    file = nc.NetCDFFile(filename,'r')
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
nrep = 2
dt = 8
nmd = 2**14
T=300
dT=0.0*T
hwcut=0.03
#--------------------------------------------------------------------------------------
#    def __init__(self,xyz,syslist,dt,nmd,harmonic=False,dyn=None):
print "initialise md"
test = md()
print md.__doc__

test.setDyn(eph.DynMat)
syslist=range(len(eph.DynMat)/3)
test.SetSyslist(syslist)
test.SetHarm(True)
test.SetMD(dt,nmd)
test.SetT(T)
print "number of degrees of freedome: ", test.nph
print "number of atoms: ", test.na
print test.nmd

#----------------------------------------------------
#debye bath
#----------------------------------------------------
#    def __init__(self,T,cats,debye,nw,dt=None,nmd=None,ml=None,mcof=2.0,gamma=None,gwl=None,K00=None,K01=None,V01=None):
    #        self.T,self.debye,self.cats = T,debye,N.array(cats,dtype='int')
# initialise a phonon bath with debye damping
# phl = phbath(300,[0,1,2,3,4,5,6,7,8,9,10,11,12],0.022,100,test.dt,test.nmd)
# phl.gnoi()
# phl.gmem()
# PP.plot(phl.noise[:,0])
# PP.show()
# ##add the bath 
# test.AddBath(phl)
# print "memory length: %s\n"%test.ml
# 
#----------------------------------------------------
#----------------------------------------------------
# electron bath
#    def __init__(self,nph,cats,Te=0.,wmax=None,nw=None,bias=0.,efric=None,exim=None,exip=None,dt=None,nmd=None):
#----------------------------------------------------
###nph = len(eph.SigL[0])
###ecats=range(13)
###eb = ebath(nph,ecats,Te=300.0,wmax=1.0,nw=500,bias=0.,efric=eph.efric,exim=eph.xim,exip=eph.xip)
###eb.genoi()
###PP.plot(eb.noise[:,0])
###PP.show()
###sys.exit()
#----------------------------------------------------
# left bath
#----------------------------------------------------
#    def __init__(self,T,cats,debye,nw,dt=None,nmd=None,ml=None,mcof=2.0,gamma=None,gwl=None,K00=None,K01=None,V01=None):
    #        self.T,self.debye,self.cats = T,debye,N.array(cats,dtype='int')
# initialise a phonon bath with debye damping
phl = phbath(T+dT,[0,1,2,3,4],hwcut,200,test.dt,test.nmd,ml=100,sig=eph.SigL,gwl=eph.wl)

phl.gnoi()
phl.gmem()

print phl.kernel.shape

#PP.plot(phl.noise[:,0])
#PP.savefig("noise_left.pdf",format="pdf")
#PP.close()
#PP.plot([N.trace(a) for a in phl.kernel])
#PP.savefig("friction_left.pdf",format="pdf")
#PP.close()
#add the bath 
test.AddBath(phl)
print "memory length: %s\n"%test.ml
print phl.ml

#----------------------------------------------------
# right bath
#----------------------------------------------------
#    def __init__(self,T,cats,debye,nw,dt=None,nmd=None,ml=None,mcof=2.0,gamma=None,gwl=None,K00=None,K01=None,V01=None):
    #        self.T,self.debye,self.cats = T,debye,N.array(cats,dtype='int')
# initialise a phonon bath with debye damping
#phl = phbath(300,[0,1,2,3,4],0.022,200,test.dt,test.nmd,ml=100,gamma=gaml,gwl=eph.wl)
phr = phbath(T-dT,[8,9,10,11,12],hwcut,200,test.dt,test.nmd,ml=100,sig=eph.SigR,gwl=eph.wl)

phr.gnoi()
phr.gmem()

print phr.kernel.shape

#PP.plot(phr.noise[:,0])
#PP.savefig("noise_right.pdf",format="pdf")
#PP.close()
#PP.plot([N.trace(a) for a in phr.kernel])
#PP.savefig("friction_right.pdf",format="pdf")
#PP.close()
#add the bath 
test.AddBath(phr)
print "memory length: %s\n"%test.ml
print phr.ml
#----------------------------------------------------

test.initialise()
test.ResetHis()
print test.ml
print len(test.baths)


power=N.zeros((test.nmd,2))
for j in range(nrep):
    for i in range(len(test.baths)):
        test.baths[i].gnoi()

    test.ResetSavepq()
    for i in range(test.nmd):
        test.vv()
    print test.qs.shape
    test.GetPower()
    #print test.power
    #PP.plot(test.power)
    #PP.show()
    power=power+test.power
    
#    
#    print "average current:"+str(N.mean(test.baths[0].cur))
#    PP.plot(test.baths[0].cur)
#    PP.savefig("curl"+str(j)+".pdf",format="pdf")
#    PP.close()
#    
#    print "average current:"+ str(N.mean(test.baths[1].cur))
#    PP.plot(test.baths[1].cur)
#    PP.savefig("curr"+str(j)+".pdf",format="pdf")
#    PP.close()
    
power=power/float(nrep)

f = open("power.dat","w")
#f.write("#k-point averaged transmission and DoS from MAMA.py\n")
#f.write("#energy    transmission    DoSL    DoSR\n")
for i in range(len(power)):
    f.write("%f     %f \n"%(power[i,0],power[i,1]))
f.close()
