#!/applications/mbrsoft/bin/python

import sys,time
import os.path
import numpy as N
from numpy import linalg as LA
#import Scientific.IO.NetCDF as nc
from netCDF4 import Dataset

import units as U
from matrix import *
from functions import *
from myfft import *
from noise import *
from spectrum import *

"""
#TODO:
1. Now we always save fhis for each bath, this may do not work for large
    structures.
"""

class md:
    """
    do langevin molecular dynamics using a modified velocity verlet method.
    The friction from both electrons and phonons are time-local now.

    The following variables are initialied when called:
        xyz         position array of atoms [x1,y1,z1,x2,y2,z2,...]
        els         list of elements ["Au","Au",...]
        nta         total number of atoms, including system and buffering atoms
        syslist     list made from indices of system atoms (python indices,
                    starting from 0)
        dt          time step 
        nmd         len of md simulation
        na          number of system atoms
        nph         number of system degrees of freedom
        dyn         dynamical matrix of the system
        harmonic    how to calculate the potential force: harmonic approximation
                    or siesta
        T           Average temperature
        savepq      whether to save md trajectories
        ml          length of memory kernel (max(all baths))

        hw          eigen frequencies
        U           eigen vectors, each column corresponds to one eigen vector,
                    for example, hw, U = np.linalg.eigh(DynMat)

        t,p,q       md time step, velocity and position vector at t

        q0,f0       q and f of previous potforce run
        
        baths       List of baths connecting to the system
        fbaths      Force from all the baths at time t
        etot        total energy at each time step
        nstep       Output atomic position after nstep MD 


    """
    def __init__(self,dt,nmd,T,syslist=None,axyz=None,harmonic=False,\
                 dyn=None,savepq=True,nrep=1,npie=8,constr=None,nstep=100,md2ang=0.06466):
        #drivers
        self.sint = None #siesta instance
        self.brennerrun=None
        self.lammpsrun=None

        self.constraint=constr
        self.nrep=nrep
        self.dt,self.nmd = dt,nmd
        self.harmonic = harmonic
        self.T=T
        self.npie=npie
        self.power = N.zeros((self.nmd,2))
        self.power2 = N.zeros((self.nmd,2))


        #var: xyz,nta,els
        self.SetXyz(axyz)
        #var: syslist,na,nph
        if syslist is not None:
            if(len(syslist) > self.nta or min(syslist) < 0 or \
               max(syslist) > self.nta-1):
                print "syslist out of range"
                sys.exit(0)
            else:
                self.syslist = N.array(syslist,dtype='int')
            #number of system atoms
            self.na = len(syslist)
            #number of system degrees of freedom
            self.nph = 3*len(syslist)
        elif axyz is not None:
            #set using axyz
            #all are system atoms
            self.syslist = N.array(range(len(axyz)),dtype='int')
            self.na = len(self.syslist)
            self.nph = 3*len(self.syslist)
        else:
            self.syslist = None
            self.na = None
            self.nph = None



        self.ml = 1
        self.t=0
        self.p=[]
        self.q=[]
        self.pinit=[]
        self.qinit=[]
        self.q0=[]
        self.f0=[]

        #qhis and phis

        #list of baths
        self.baths=[]
        self.fhis=[]
        self.fbaths=[]

        self.etot=N.zeros(nmd)

        #vars: dyn,hw,U,nph
        self.setDyn(dyn)

        #var: ps,qs,power,savepq
        self.ResetSavepq(savepq)
#--------------------------------------------------------------
#Add tracking of atomic trajectories by Li Gen.
        self.md2ang = md2ang
        self.mass = []
        self.get_atommass()
        self.conv = self.md2ang*N.array([3*[1.0/N.sqrt(mass)]
                                         for mass in self.mass]).flatten()
        self.nstep = nstep

    def get_atommass(self):
        for atomsname in self.els:
            for key, value in U.AtomicMassTable.items():
                if atomsname==key:
                    self.mass.append(value)
#--------------------------------------------------------------

    def info(self):
        print "--------------------------------------------\n"
        print "Basis information of the MD simulation:\n\n"
        print "--------------------------------------------\n"
        print "Harmonic force: "+str(self.harmonic)+"\n"
        print "System atom number:"+str(self.na)+"\n"
        print "MD time step:"+str(self.dt)+"\n"
        print "MD number of steps:"+str(self.nmd)+"\n"
        print "MD memory kernel length:"+str(self.ml)+"\n"
        print "Number of baths attached:"+str(len(self.baths))+"\n\n\n"

        if self.dyn is None:
            print "md.info: No dynamical matrix input!"
            #sys.exit()


    def ResetSavepq(self,savepq=True):
        self.savepq = savepq
        if self.nmd is not None and self.nph is not None:
            self.ps = N.zeros((self.nmd,self.nph))
            self.qs = N.zeros((self.nmd,self.nph))
            #self.power = N.zeros((self.nmd,2))
        else:
            print "md.ResetSavepq: nmd or nph is not set!"

    def energy(self):
        return 0.5*mm(self.p,self.p)+0.5*mm(self.q,self.dyn,self.q)
    def energy(self):
        """
        kinetic energy
        """
        return 0.5*mm(self.p,self.p)

    def AddBath(self,bath):
        """
        Adding a bath
        """
        if self.dt != bath.dt:
            print "md.AddBath: md time step dt not consistent!"
            sys.exit()
        if self.nmd != bath.nmd:
            print "md.AddBath: number of md steps nmd not consistent!"
            sys.exit()
        self.baths.append(bath)
        #make sure we save enought memory for all the baths
        if(bath.ml > self.ml):
            self.ml = bath.ml
        self.fbaths.append(N.zeros(self.nph))
        #force history
        self.fhis.append(N.zeros((self.nmd,self.nph)))

    def SetT(self,T):
        self.T=T

    def SetMD(self,dt,nmd):
        self.dt,self.nmd=dt,nmd
        self.etot=N.zeros(nmd)

    def SetHarm(self,harmonic):
        self.harmonic = harmonic

    def SetXyz(self,axyz):
        if axyz is not None:
            print "md.SetXyz:Seting xyz and nta"
            self.xyz = N.array([a[1:] for a in axyz],dtype='d').flatten()
            self.els = [a[0] for a in axyz]
            self.nta = len(axyz)
        else:
            self.xyz = None
            self.els = None
            self.nta = None

    def SetSyslist(self,syslist):
        print "md.SetXyz:Seting syslist"
        self.syslist= N.array(syslist)
        #number of system atoms
        self.na = len(syslist)
        #number of system degrees of freedom
        self.nph = 3*len(syslist)
        if self.xyz is not None:
            if len(self.syslist) > self.nta:
                print "md.SetSyslist:system atom number larger than total atom number"
                sys.exit()

    def setDyn(self,dyn=None):
        """
        set up the dynamical matrix of the system
        """
        if dyn is not None:
            print "md.setDyn: getting dynamical matrix"
            ndyn = N.array(dyn)
            print "md.setDyn: checking dynamical matrix"
            n = chkShape(ndyn)
            if self.nph is not None and self.nph != n:
                print "md.setDyn: the dimension of dynamical matrix is wrong!"
                sys.exit(0)
            self.nph = n
            print "md.setDyn: symmetrizing dynamical matrix"
            self.dyn = symmetrize(ndyn)

            av,au = LA.eigh(self.dyn)
            if min(av)<0:
                print "md.setDyn: there are negative frequencies"
                print "md.setDyn: I will remove them"
                avn=0.*av
                for i in range(len(av)):
                    if av[i]<0:
                        avn[i]=0.
                    else:
                        avn[i]=av[i]
                av=avn
            self.hw = N.array(map(N.real,map(N.sqrt,av)))
            self.U = N.array(au)
            self.dyn=mm(self.U,N.diag(N.array(av)),N.transpose(self.U))
            #if min(av)>=0:
            #    print "the dynmat should not change much"
            #    print "max diff. of dynmatrix:", abs(self.dyn-ndyn).sum()
            print "md.setDyn: Done"
        else:
            self.dyn = None
            self.hw = None
            self.U = None
            print "md.setDyn: no dynamical matrix provided!"
            #sys.exit()
            

    def initialise(self):
        """
        initial displacement and velocity from the dynamical matrix
        """
        self.t = 0
        if self.dyn is None:
            print "md.initial: no dynamical matrix!!"
            #sys.exit()
            print "p,q set to 0"
            self.p=N.zeros(self.nph)
            self.q=N.zeros(self.nph)
            self.pinit=N.zeros(self.nph)
            self.qinit=N.zeros(self.nph)
        else:
            av=self.hw
            au=self.U

            dis=N.zeros(len(av))
            vel=N.zeros(len(av))
            for i in range(len(av)):
                #cutoff energy 0.005 eV
                #do not initialise motion due to slow modes
                #because it may gives large displacement
                if av[i] < 0.01:
                    am=0.0
                else:
                    am=((bose(av[i],self.T)+0.5)*2.0/av[i])**0.5
                r=N.random.rand()
                dis = dis + au[:,i]*am*N.cos(2.*N.pi*r)
                vel = vel - av[i]*au[:,i]*am*N.sin(2.*N.pi*r)

                dis = ApplyConstraint(dis,self.constraint)
                vel = ApplyConstraint(vel,self.constraint)

            self.p=vel
            self.q=dis
            self.pinit = vel
            self.qinit = dis

    def ResetHis(self):
        """
        set history list of the friction kernel to zeros
        """
        if self.nph is not None and self.ml is not None:
            self.qhis=N.zeros((self.ml,self.nph))
            self.phis=N.zeros((self.ml,self.nph))
        else:
            print "self.nph and self.ml are not set!"
            sys.exit()

    def GetPower(self):
        """
        calculate the power spectrum from the MD trajectories.
        """
        if not self.savepq:
            print "md.GetPower: trajectories not saved!"
            print "md.GetPower: you need to set savepq to True!"
            sys.exit()
        print "md.GetPower: generate power spectrum from trajectories!"
        self.power = powerspec(self.qs,self.dt,self.nmd)
        self.power2 = powerspec2(self.ps,self.dt,self.nmd)

    def vv(self,id):
        """
        velocity-verlet method integrator
        """
        #print "velocity-verlet integrator"
        t,p,q = self.t,self.p,self.q
        t = int(t)
        if self.savepq:
            self.ps[t%self.nmd] = p
            self.qs[t%self.nmd] = q

        #total energy
        #self.etot = N.append(self.etot,self.energy())
        self.etot[t%self.nmd] = self.energy()

        #update history here
        self.qhis=rpadleft(self.qhis,q)
        self.phis=rpadleft(self.phis,p)

        #calculate displacement at next time
        f = self.force(t,p,q,0)
        pthalf = p + f*self.dt/2.0
        qtt = q + p*self.dt + f*self.dt**2/2.0
        
        #evaluate current 
        for i in range(len(self.baths)):
            #self.baths[i].cur=N.append(self.baths[i].cur,mm(self.fbaths[i],p))
            self.baths[i].cur[t%self.nmd] = mm(self.fbaths[i],p)
            self.fhis[i][t%self.nmd] = self.fbaths[i]


        #calculate velocity at next time
        f = self.force(t,pthalf,qtt,1)
        ptt1=pthalf+self.dt*f/2.0
        f=self.force(t,ptt1,qtt,1)
        ptt2=pthalf+self.dt*f/2.0

        
        #constraint
        ptt2=ApplyConstraint(ptt2,self.constraint)
        qtt=ApplyConstraint(qtt,self.constraint)

        t=t+1
        self.t,self.p,self.q = t,ptt2,qtt
#-------------------------------------------------------------------------------------
#Add tracking of atomic trajectories by Li Gen.
        if self.t == 1 or self.t % self.nstep == 0:
            with open('trajectories'+str(self.T)+"."+str(id)+'.ani', 'a') as fileobject:
            #with open('OptimizationMJ'+str(self.t)+'.ang', 'w') as fileobject:
                fileobject.write(str(len(self.els)))
                fileobject.write('\n')
                fileobject.write('Timestep'+'   '+str(self.t))
                fileobject.write('\n')
                for ip in range(len(self.els)):
                    fileobject.write(str(self.els[ip])+'    ')
                    fileobject.write(str(self.xyz[ip*3]+self.conv[ip*3]*self.q[ip*3])+'   ')
                    fileobject.write(str(self.xyz[ip*3+1]+self.conv[ip*3+1]*self.q[ip*3+1])+'   ')
                    fileobject.write(str(self.xyz[ip*3+2]+self.conv[ip*3+2]*self.q[ip*3+2])+'   ')
                    fileobject.write('\n')
                fileobject.write('\n')
#-------------------------------------------------------------------------------------

    def force(self,t,p,q,id=0):
        """
        force due to the dynamical matrix
        """
        it = t+id

        #print "potential force"
        pf=self.potforce(q)

        #apply constraint
        #pf=ApplyConstraint(pf,self.constraint)

        #print "friction and random force "
        if id == 0:
            tphis = self.phis
            tqhis = self.qhis
        else:
            tphis = rpadleft(self.phis,p)
            tqhis = rpadleft(self.qhis,q)
        for i in range(len(self.baths)):
            self.fbaths[i] = self.baths[i].bforce(it,tphis,tqhis)
            pf = pf + self.fbaths[i]
        
        return pf

    def potforce(self,q):
        """
        potential force
        """
        if sameq(q,self.q0):
            return self.f0
        else:
            #use dynamical matrix
            if self.sint is None:
                f=-mdot(self.dyn,q)
            #use siesta force 
            else:
                slist=self.syslist
                if len(q)/3 != len(slist):
                    print "md:potforce: length error!"
                    sys.exit()
                extq = N.zeros(len(self.xyz))
                for i in range(len(slist)):
                    extq[3*slist[i]:3*(slist[i]+1)] = q[3*i:3*(i+1)]
                fa = self.sint.force(extq)
                f = N.zeros(len(q))
                for i in range(len(f)/3):
                    f[i*3:(i+1)*3] = fa[slist[i]*3:(slist[i]+1)*3]
            #save 
            self.q0=q
            self.f0=f
            return f

    def potforce(self,q):
        """
        Tue's version including brenner
        """
        if sameq(q,self.q0):
            return self.f0
        else:
            #use dynamical matrix
            if self.sint is None:
                if self.brennerrun is None:
                    f=-mdot(self.dyn,q)
                else:
                    f=self.brennerrun.force(q)
            #use siesta force 
            else:
                slist=self.syslist
                if len(q)/3 != len(slist):
                    print "md:potforce: length error!"
                    sys.exit()
                extq = N.zeros(len(self.xyz))
                for i in range(len(slist)):
                    extq[3*slist[i]:3*(slist[i]+1)] = q[3*i:3*(i+1)]
                fa = self.sint.force(extq)
                f = N.zeros(len(q))
                for i in range(len(f)/3):
                    f[i*3:(i+1)*3] = fa[slist[i]*3:(slist[i]+1)*3]
            #save 
            self.q0=q
            self.f0=f
            return f

    def potforce(self,q):
        """
        Potential force from drivers
        8Nov2018:
        Now we have the following 4 drivers:
            Siesta, Brenner, Lammps, harmonic
        
        q is an array of displacements of the atoms in self.syslist
        """
        #self.q0 and f0 are the displacment and
        #force of last call. If the displacement
        #did not change, use the old force.
        if sameq(q,self.q0):
            return self.f0

        if len(q)/3 != len(self.syslist):
            print "md:potforce: length error!"
            sys.exit()
        extq = N.zeros(len(self.xyz))
        for i in range(len(self.syslist)):
            extq[3*self.syslist[i]:3*(self.syslist[i]+1)] = q[3*i:3*(i+1)]

        #search for possible drivers
        #use siesta force 
        if self.sint is not None:
            fa = self.sint.force(extq)
            f = N.zeros(len(q))
            for i in range(len(f)/3):
                f[i*3:(i+1)*3] = fa[self.syslist[i]*3:(self.syslist[i]+1)*3]
        #use brenner force 
        elif self.brennerrun is not None:
            f=self.brennerrun.force(q)
        #use lammps force 
        elif  self.lammpsrun is not None:
            fa=self.lammpsrun.force(extq)
            f = N.zeros(len(q))
            for i in range(len(f)/3):
                f[i*3:(i+1)*3] = fa[self.syslist[i]*3:(self.syslist[i]+1)*3]
        #use dynamical matrix 
        elif self.dyn is not None:
            f=-mdot(self.dyn,q)
        else:
            print "no driver, no md"
            sys.exit()

        #save 
        self.q0=q
        self.f0=f
        return f


    def AddSint(self,sint):
        """
        add siesta instance
        """
        #if sint.xyz != self.xyz:
        #    print "md.AddSint: xyz not consistent!"
        #    sys.exit()
        #if sint.els != self.els:
        #    print "md.AddSint: els not consistent!"
        #    sys.exit()
        self.sint = sint
        #print "Starting Siesta Server..."
        #self.sint.start()

    def AddBint(self,bint):
        """
        add Brenner instance
        """
        self.brennerrun = bint


    def AddLMPint(self,lint):
        """
        add Lammps instance
        """
        self.lammpsrun = lint


    def Run(self):
        """
        define nrep,npie
        """
        #initialise t,p,q
        #if find old unfinished data,
        #they are over-written
        self.initialise()
        #reset qhis,phis
        self.ResetHis()
        self.info()
        
        #loop over independent md runs
        for j in range(self.nrep):   
            print "MD run: "+str(j)+"\n"
            fn="MD"+str(j)+".nc"
            fnm="MD"+str(j-1)+".nc"

            if os.path.isfile(fn):
                print "find file: "+fn+"\n"
                ipie = int(ReadNetCDFVar(fn,'ipie')[0])
                if(ipie+1 < self.npie):
                    print "unfinished run"
                    print "reading resume information"
                    self.p = ReadNetCDFVar(fn,'p')
                    self.q = ReadNetCDFVar(fn,'q')
                    self.t = ReadNetCDFVar(fn,'t')[0]
                    self.qhis = ReadNetCDFVar(fn,'qhis')
                    self.phis = ReadNetCDFVar(fn,'phis')
                    self.power =ReadNetCDFVar(fn,'power')
                    self.power2 =ReadNetCDFVar(fn,'power2')
                    if self.savepq:
                        self.qs = ReadNetCDFVar(fn,'qs')
                        self.ps = ReadNetCDFVar(fn,'ps')
                    for i in range(len(self.baths)):
                        self.baths[i].noise=ReadNetCDFVar(fn,'noise'+str(i))

                elif(ipie+1 == self.npie):
                    print "finished run"
                    self.power =ReadNetCDFVar(fn,'power')
                    self.power2 =ReadNetCDFVar(fn,'power2')
                    self.t = ReadNetCDFVar(fn,'t')[0]
                    continue
                else:
                    print "ipie error"
                    sys.exit()
            else:
                print "new run"
                ipie=-1  #yes,-1
                if os.path.isfile(fnm):
                    print "reading history from previous run"
                    self.p = ReadNetCDFVar(fnm,'p')
                    self.q = ReadNetCDFVar(fnm,'q')
                    qhis0 = ReadNetCDFVar(fnm,'qhis')
                    phis0 = ReadNetCDFVar(fnm,'phis')
                    self.t = ReadNetCDFVar(fnm,'t')[0]
                    if qhis0.shape == self.qhis.shape and\
                            phis0.shape == self.phis.shape:
                        self.qhis = qhis0
                        self.phis = phis0
                #noise generation
                for i in range(len(self.baths)): #loop over baths
                    self.baths[i].gnoi()
                    #print N.shape(self.baths[i].noise)
                    #stppp
        
                #reset qs and ps to zero
                self.ResetSavepq()
        
            #loop over md steps
            ipie1=ipie+1
            iss=ipie1+N.array(range(self.npie-ipie1))
            for i in iss:
                for jj in range(self.nmd/self.npie):
                    self.vv(j)
                self.dump(i,j)

            #power spectrum
            power=N.copy(self.power)
            power2=N.copy(self.power2)
            self.GetPower()
            power=(power*j+self.power)/float(j+1)
            power2=(power2*j+self.power2)/float(j+1)
            self.power=N.copy(power)
            self.power2=N.copy(power2)

            #dump again, to make sure power is all right
            self.dump(i,j)
            

            #----------------------------------------------------------------
            #heat current
            #----------------------------------------------------------------
            for ii in range(len(self.baths)):
                fk = open("kappa."+str(self.T)+"."+"bath"+str(ii)+".run"+str(j)+".dat","w")
                #write average current
                fk.write("%i %f    %f \n"%(j,self.T,N.mean(self.baths[ii].cur)*U.curcof))
                fk.close()

            #----------------------------------------------------------------
            #power spectrum
            #----------------------------------------------------------------
            f = open("power."+str(self.T)+"."+"run"+str(j)+".dat","w")
            #f.write("#k-point averaged transmission and DoS from MAMA.py\n")
            #f.write("#energy    transmission    DoSL    DoSR\n")
            for i in range(len(self.power)):
                #only write out power spectrum upto 1.5max(hw)
                if self.hw is not None:
                    if(self.power[i,0] < 1.5*max(self.hw)):
                        f.write("%f     %f \n"%(self.power[i,0],self.power[i,1]))
                    else:
                        break
                else:
                    f.write("%f     %f \n"%(self.power[i,0],self.power[i,1]))
            f.close()

            #----------------------------------------------------------------
            #power spectrum from velocity
            #----------------------------------------------------------------
            f = open("power2."+str(self.T)+"."+"run"+str(j)+".dat","w")
            for i in range(len(self.power2)):
                #only write out power spectrum upto 1.5max(hw)
                if self.hw is not None:
                    if(self.power2[i,0] < 1.5*max(self.hw)):
                        f.write("%f     %f \n"%(self.power2[i,0],self.power2[i,1]))
                    else:
                        break
                else:
                    f.write("%f     %f \n"%(self.power2[i,0],self.power2[i,1]))
            f.close()

    def dump(self,ipie,id):
        """
        dump md information
        """
        outfile="MD"+str(id)+".nc"
        NCfile = Dataset(outfile,'w','Created '+time.ctime(time.time()))
        NCfile.title = 'Output from md.py'
        NCfile.createDimension('nph',self.nph)
        #NCfile.createDimension('na',self.na)
        NCfile.createDimension('one',1)
        NCfile.createDimension('two',2)
        NCfile.createDimension('mem',self.ml)
        NCfile.createDimension('nmd',self.nmd)
        NCfile.createDimension('nnmd',None)
        for i in range(len(self.baths)):
            NCfile.createDimension('n'+str(i),self.baths[i].nc)

        #els
        #Write2NetCDFFile(NCfile,self.els,'elements',('na',),units='')

        #noise series
        for i in range(len(self.baths)):
            #NCfile.createDimension('n'+str(i),self.baths[i].nc)
            Write2NetCDFFile(NCfile,self.baths[i].noise,'noise'+str(i),\
                            ('nnmd','n'+str(i),),units='')
            Write2NetCDFFile(NCfile,self.fhis[i],'fhis'+str(i),('nnmd','nph',),units='')


        #save all the histories of p,q or not
        if self.savepq:
            Write2NetCDFFile(NCfile,self.ps,'ps',('nnmd','nph',),units='')
            Write2NetCDFFile(NCfile,self.qs,'qs',('nnmd','nph',),units='')

        #power spectrum
        Write2NetCDFFile(NCfile,self.power,'power',('nnmd','two',),units='')
        Write2NetCDFFile(NCfile,self.power2,'power2',('nnmd','two',),units='')
        #energy
        Write2NetCDFFile(NCfile,self.etot,'energy',('nnmd',),units='')
        #velocity
        Write2NetCDFFile(NCfile,self.p,'p',('nph',),units='')
        #displacement
        Write2NetCDFFile(NCfile,self.q,'q',('nph',),units='')
        #current time
        Write2NetCDFFile(NCfile,self.t,'t',('one',),units='')
        #current segment
        Write2NetCDFFile(NCfile,ipie,'ipie',('one',),units='')
        #memory kernel for p
        Write2NetCDFFile(NCfile,self.phis,'phis',('mem','nph',),units='')
        #memory kernel for q
        Write2NetCDFFile(NCfile,self.qhis,'qhis',('mem','nph',),units='')

        NCfile.close()


#--------------------------------------------------------------------------------------
#misc driver routines
def Write2NetCDFFile(file,var,varLabel,dimensions,units=None,description=None):
    #print 'Write2NetCDFFile:', varLabel, dimensions
    tmp = file.createVariable(varLabel,'d',dimensions)
    tmp[:] = var
    if units: tmp.units = units
    if description: tmp.description = description

def ReadNetCDFVar(file,var):
    print "ReadNetCDFFile: reading "+ var
    f = Dataset(file,'r')
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
    #if the difference is less than 10e-10
    #use the old force
    if(max(abs(dif)) < 10e-10):
        return True
    else:
        return False


def ApplyConstraint(f,constr=None):
    """
    apply constraint to the force

    constr is an array of vectors
    This subroutine zerofy the force along each vector
    """
    if constr is None:
        return f
    #check dimensions
    f=N.array(f)
    dd=len(f)
    constr=N.array(constr)
    n,d=N.shape(constr)
    if d!=dd:
        print "ApplyConstraint:shape error"

    nf=1.0*f
    for i in range(n):
        nn=LA.norm(constr[i])
        if nn != 0:
            constr[i]=constr[i]/nn
        nf=nf-N.dot(f,constr[i])*constr[i]
    return nf


#--------------------------------------------------------------------------------------
#testing
#
#--------------------------------------------------------------------------------------

if __name__=="__main__":
    #import matplotlib.pyplot as PP
#--------------------------------------------------------------------------------------
    print "testing md - begin"
    test = md([[0,0,0],[0,0,1],[0,0,-1]],[0],0.5,2**16,1.,1.,1.)
    print md.__doc__
    print "number of degrees of freedome: ", test.nph
    print "number of atoms: ", test.na

    #set dynamical matrix
    a = N.diag(1.+N.zeros(3))
    test.setDyn(a)
    print test.dyn

    #set baths
    test.setEmat(efric=a,exim=a,exip=a)
    print test.ebath
    test.setPhbath(Tl=300.,gammal=[[[0.1]],[[0.1]],[[0.1]]],wll=[0.0,0.1,0.2],ll=[1],Tr=300.,gammar=[[[0.1]],[[0.1]],[[0.1]]],wlr=[0.0,0.1,0.2],lr=[1])
    print test.lbath
    print test.rbath

    if test.ebath:
        print "generating electron noise"
        test.eno =test.genoi()
        print test.eno.shape
        print N.mean(test.eno,0)
        PP.plot(test.eno[:,0])
        PP.show()
    if test.lbath:
        print "generating left phonon noise"
        test.phnol =test.glnoi()
        print test.phnol.shape
        print N.mean(test.phnol,0)
        PP.plot(test.phnol[:,0])
        PP.show()
    if test.rbath:
        print "generating right phonon noise"
        test.phnor=test.grnoi()
        print test.phnor.shape
        print N.mean(test.phnor,0)
        PP.plot(test.phnor[:,0])
        PP.show()

    print "testing md - end"

#--------------------------------------------------------------------------------------
    print "initial conditions - begin"
    dis,vel=test.initialise()
    print dis
    print vel
    print "initial conditions - end"
#--------------------------------------------------------------------------------------
    print "testing vargau - begin"
    a = N.array([[1.,0.1,0.3],[0.1,2.,3.],[1.,2.,3.]])
    aa = a+dagger(a)
    print aa
    eval,evec = LA.eigh(aa)  #
    #column evec[:,i] corresponds to eval[i]
    print eval
    print evec

    aa = vargau(eval,evec)
    print aa
    print "testing vargau - end"

#--------------------------------------------------------------------------------------
    print "testing nonequ - begin"
    x0 = -0.01*N.arange(100)
    x00 = x0[::-1]  #reverse of an array
    x = N.concatenate((x00[:-1],-x0[1:]),axis=0)
    print "x is:", x
    T = 0.0
    y = [nonequp(xi,1.0,T)+nonequm(xi,1.0,T) for xi in x]
    PP.plot(x,y)
    print "y=nonequp()+nonequm()"
    print "At T=0, it should be a symmetric linear curve"
    PP.show()
    print "testing nonequ - end"
#--------------------------------------------------------------------------------------
