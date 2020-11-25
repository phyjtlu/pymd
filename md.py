#!/usr/bin/env python
# -*- coding: utf-8 -*
import os
import sys
import time

import numpy as N
from netCDF4 import Dataset
from numpy import linalg as LA
from tqdm import tqdm

import units as U
from functions import bose, chkShape, mdot, powerspecp, rpadleft, symmetrize


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
        saveq       whether to save md trajectories
        savep       whether to save md trajectories
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
        saveall     Whether to save all information to .nc file
        rmnc        Remove NC files after calculation
    """

    def __init__(self, dt, nmd, T, syslist=None, axyz=None, harmonic=False, dyn=None, nstart=0, nstop=1, npie=1, md2ang=0.06466):
        self.nstart = nstart
        self.nstop = nstop
        self.dt, self.nmd = dt, nmd
        self.harmonic = harmonic
        self.T = T
        self.npie = npie
        self.saveall = False
        self.savep = False
        self.saveq = False
        self.rmnc = False
        self.nstep = None
        self.pforce = None
        self.constraint = None
        self.atomlist = None
        # var: xyz,nta,els
        self.SetXyz(axyz)
        # var: syslist,na,nph
        if syslist is not None:
            if(len(syslist) > self.nta or min(syslist) < 0 or
               max(syslist) > self.nta-1):
                print("syslist out of range")
                sys.exit(0)
            else:
                self.syslist = N.array(syslist, dtype='int')
            # number of system atoms
            self.na = len(syslist)
            # number of system degrees of freedom
            self.nph = 3*len(syslist)
        elif axyz is not None:
            # set using axyz
            # all are system atoms
            self.syslist = N.array(range(len(axyz)), dtype='int')
            self.na = len(self.syslist)
            self.nph = 3*len(self.syslist)
        else:
            self.syslist = None
            self.na = None
            self.nph = None

        self.ml = 1
        self.t = 0
        self.p = []
        self.q = []
        self.pinit = []
        self.qinit = []
        self.q0 = []
        self.f0 = []

        #qhis and phis

        # list of baths
        self.baths = []
        self.fhis = []
        self.fbaths = []

        self.etot = N.zeros(nmd)

        # vars: dyn,hw,U,nph
        self.setDyn(dyn)

        # var: ps,qs,power
        self.ResetSavepq()
        self.md2ang = md2ang
        self.mass = []
        self.get_atommass()
        self.conv = self.md2ang*N.array([3*[1.0/N.sqrt(mass)]
                                         for mass in self.mass]).flatten()

    def get_atommass(self):
        for atomsname in self.els:
            for key, value in list(U.AtomicMassTable.items()):
                if atomsname == key:
                    self.mass.append(value)

    def info(self):
        print("--------------------------------------------")
        print("Basis information of the MD simulation:")
        # print("Harmonic force: "+str(self.harmonic))
        print("System atom number:"+str(self.na))
        print("MD time step:"+str(self.dt))
        print("MD number of steps:"+str(self.nmd))
        print("MD memory kernel length:"+str(self.ml))
        print("Number of baths attached:"+str(len(self.baths))+"\n")

        # if self.dyn is None:
        #    print("md.info: No dynamical matrix input")
        # sys.exit()

    def ResetSavepq(self):
        if self.savep and self.nmd is not None and self.nph is not None:
            self.ps = N.zeros((self.nmd, self.nph))
            self.power = N.zeros((self.nmd, 2))
        if self.saveq and self.nmd is not None and self.nph is not None:
            self.qs = N.zeros((self.nmd, self.nph))
        print("md.save all momentum: %r" % self.savep)
        print("md.save all postions: %r" % self.saveq)

    # def energy(self):
    #    return 0.5*mdot(self.p,self.p)+0.5*mdot(self.q,self.dyn,self.q)
    def energy(self):
        """
        kinetic energy
        """
        return 0.5*mdot(self.p, self.p)

    def AddBath(self, bath):
        """
        Adding a bath
        """
        if self.dt != bath.dt:
            print("md.AddBath: md time step dt not consistent")
            sys.exit()
        if self.nmd != bath.nmd:
            print("md.AddBath: number of md steps nmd not consistent")
            sys.exit()
        self.baths.append(bath)
        # make sure we save enought memory for all the baths
        if(bath.ml > self.ml):
            self.ml = bath.ml
        self.fbaths.append(N.zeros(self.nph))
        # force history
        self.fhis.append(N.zeros((self.nmd, self.nph)))

    def AddPowerSection(self, atomlist):
        self.atomlist = atomlist
        self.poweratomlist = N.empty((len(self.atomlist), self.nmd, 2))

    def AddConstr(self, constr):
        self.constraint = constr

    def CalPowerSpec(self, cal=True):
        self.savep = cal

    def CalAveStruct(self, cal=True):
        self.saveq = cal

    def SaveAll(self, save=True):
        self.saveall = save

    def SaveTraj(self, nstep=100):
        self.nstep = nstep

    def RemoveNC(self, rmnc=True):
        self.rmnc = rmnc

    def SetT(self, T):
        self.T = T

    def SetMD(self, dt, nmd):
        self.dt, self.nmd = dt, nmd
        self.etot = N.zeros(nmd)

    def SetHarm(self, harmonic):
        self.harmonic = harmonic

    def SetXyz(self, axyz):
        if axyz is not None:
            print("md.SetXyz:Seting xyz and nta")
            self.xyz = N.array([a[1:] for a in axyz], dtype='d').flatten()
            self.els = [a[0] for a in axyz]
            self.nta = len(axyz)
        else:
            self.xyz = None
            self.els = None
            self.nta = None

    def SetSyslist(self, syslist):
        print("md.SetXyz:Seting syslist")
        self.syslist = N.array(syslist)
        # number of system atoms
        self.na = len(syslist)
        # number of system degrees of freedom
        self.nph = 3*len(syslist)
        if self.xyz is not None:
            if len(self.syslist) > self.nta:
                print("md.SetSyslist:system atom number larger than total atom number")
                sys.exit()

    def setDyn(self, dyn=None):
        """
        set up the dynamical matrix of the system
        """
        if dyn is not None:
            print("md.setDyn: getting dynamical matrix")
            ndyn = N.array(dyn)
            print("md.setDyn: checking dynamical matrix")
            n = chkShape(ndyn)
            if self.nph is not None and self.nph != n:
                print("md.setDyn: the dimension of dynamical matrix is wrong")
                sys.exit(0)
            self.nph = n
            print("md.setDyn: symmetrizing dynamical matrix")
            self.dyn = symmetrize(ndyn)

            av, au = LA.eigh(self.dyn)
            if min(av) < 0:
                print("md.setDyn: there are negative frequencies")
                print("md.setDyn: I will remove them")
                avn = 0.*av
                for i in range(len(av)):
                    if av[i] < 0:
                        avn[i] = 0.
                    else:
                        avn[i] = av[i]
                av = avn
            self.hw = N.array(list(map(N.real, list(map(N.sqrt, av)))))
            self.U = N.array(au)
            self.dyn = mdot(self.U, N.diag(N.array(av)), N.transpose(self.U))
            # if min(av)>=0:
            #    print "the dynmat should not change much"
            #    print "max diff. of dynmatrix:", abs(self.dyn-ndyn).sum()
            print("md.setDyn: Done")
        else:
            self.dyn = None
            self.hw = [1.0]
            self.U = None
            print("md.setDyn: no dynamical matrix provided")
            print("Set max eigen frequencies 1.0eV")
            # sys.exit()

    def initialise(self):
        """
        initial displacement and velocity from the dynamical matrix
        """
        self.t = 0
        if self.dyn is None:
            print("md.initial: no dynamical matrix")
            # sys.exit()
            print("p,q set to 0")
            self.p = N.zeros(self.nph)
            self.q = N.zeros(self.nph)
            self.pinit = N.zeros(self.nph)
            self.qinit = N.zeros(self.nph)
        else:
            av = self.hw
            au = self.U

            dis = N.zeros(len(av))
            vel = N.zeros(len(av))
            for i in range(len(av)):
                # cutoff energy 0.005 eV
                # do not initialise motion due to slow modes
                # because it may gives large displacement
                if av[i] < 0.01:
                    am = 0.0
                else:
                    am = ((bose(av[i], self.T)+0.5)*2.0/av[i])**0.5
                r = N.random.rand()
                dis = dis + au[:, i]*am*N.cos(2.*N.pi*r)
                vel = vel - av[i]*au[:, i]*am*N.sin(2.*N.pi*r)

                dis = ApplyConstraint(dis, self.constraint)
                vel = ApplyConstraint(vel, self.constraint)

            self.p = vel
            self.q = dis
            self.pinit = vel
            self.qinit = dis

    def ResetHis(self):
        """
        set history list of the friction kernel to zeros
        """
        if self.nph is not None and self.ml is not None:
            self.qhis = N.zeros((self.ml, self.nph))
            self.phis = N.zeros((self.ml, self.nph))
        else:
            print("self.nph and self.ml are not set")
            sys.exit()

    def GetPower(self):
        """
        calculate the power spectrum from the MD trajectories.
        """
        print("md.GetPower: generate power spectrum from trajectories.")
        self.power = powerspecp(self.ps, self.dt, self.nmd)
        if self.atomlist is not None:
            for layers in range(len(self.atomlist)):
                self.poweratomlist[layers] = powerspecp(
                    self.ps[:, self.atomlist[layers]], self.dt, self.nmd)

    def vv(self, id):
        """
        velocity-verlet method integrator
        """
        # print "velocity-verlet integrator"
        t, p, q = self.t, self.p, self.q
        t = int(t)
        if self.savep:
            self.ps[t % self.nmd] = p
        if self.saveq:
            self.qs[t % self.nmd] = q

        # total energy
        #self.etot = N.append(self.etot,self.energy())
        self.etot[t % self.nmd] = self.energy()

        # update history here
        self.qhis = rpadleft(self.qhis, q)
        self.phis = rpadleft(self.phis, p)

        # calculate displacement at next time
        f = self.force(t, p, q, 0)
        pthalf = p + f*self.dt/2.0
        qtt = q + p*self.dt + f*self.dt**2/2.0

        # evaluate current
        for i in range(len(self.baths)):
            # self.baths[i].cur=N.append(self.baths[i].cur,mdot(self.fbaths[i],p))
            self.baths[i].cur[t % self.nmd] = mdot(self.fbaths[i], p)
            self.fhis[i][t % self.nmd] = self.fbaths[i]

        # calculate velocity at next time
        f = self.force(t, pthalf, qtt, 1)
        ptt1 = pthalf+self.dt*f/2.0
        f = self.force(t, ptt1, qtt, 1)
        ptt2 = pthalf+self.dt*f/2.0

        # constraint
        ptt2 = ApplyConstraint(ptt2, self.constraint)
        qtt = ApplyConstraint(qtt, self.constraint)

        t = t+1
        self.t, self.p, self.q, self.f = t, ptt2, qtt, f

    def force(self, t, p, q, id=0):
        """
        force due to the dynamical matrix
        """
        it = t+id

        # print "potential force"
        pf = self.potforce(q)

        # apply constraint
        # pf=ApplyConstraint(pf,self.constraint)

        # print "friction and random force "
        if id == 0:
            tphis = self.phis
            tqhis = self.qhis
        else:
            tphis = rpadleft(self.phis, p)
            tqhis = rpadleft(self.qhis, q)
        for i in range(len(self.baths)):
            self.fbaths[i] = self.baths[i].bforce(it, tphis, tqhis)
            pf = pf + self.fbaths[i]
        return pf

    def potforce(self, q):
        """
        Potential force from drivers
        Now we have the following 4 drivers:
            Siesta, Brenner, Lammps, harmonic

        q is an array of displacements of the atoms in self.syslist
        """
        # TODO Simplified code
        # self.q0 and f0 are the displacment and
        # force of last call. If the displacement
        # did not change, use the old force.
        if sameq(q, self.q0):
            return self.f0

        # if len(q)/3 != len(self.syslist):
        #     print("md:potforce: length error")
        #     sys.exit()

        # search for possible drivers
        if self.potforce is not None:
            f = self.pforce.force(q)
        # use dynamical matrix
        elif self.dyn is not None:
            f = -1*mdot(self.dyn, q)
        else:
            print("no driver, no md")
            sys.exit()
        # save
        self.q0 = q
        self.f0 = f
        return f

    def AddPotential(self, pint):
        """
        add siesta, lammps or Brenner instance
        """
        self.pforce = pint

    def Run(self):
        """
        define nrep,npie
        """
        # initialise t,p,q
        # if find old unfinished data,
        # they are over-written
        self.initialise()
        # reset qhis,phis
        self.ResetHis()
        self.info()

        # loop over independent md runs
        for j in range(self.nstart, self.nstop):
            print("MD run: "+str(j)+"\n")
            fn = "MD"+str(j)+".nc"
            fnm = "MD"+str(j-1)+".nc"

            if os.path.isfile(fn):
                print("find file: "+fn+"\n")
                ipie = int(ReadNetCDFVar(fn, 'ipie')[0])
                if(ipie+1 < self.npie):
                    print("unfinished run")
                    print("reading resume information")
                    self.p = ReadNetCDFVar(fn, 'p')
                    self.q = ReadNetCDFVar(fn, 'q')
                    self.t = ReadNetCDFVar(fn, 't')[0]
                    self.qhis = ReadNetCDFVar(fn, 'qhis')
                    self.phis = ReadNetCDFVar(fn, 'phis')
                    if self.savep:
                        self.power = ReadNetCDFVar(fn, 'power')
                        if self.atomlist is not None:
                            self.poweratomlist = ReadNetCDFVar(
                                fn, 'poweratomlist')
                    if self.saveall and self.q and self.p:
                        self.qs = ReadNetCDFVar(fn, 'qs')
                        self.ps = ReadNetCDFVar(fn, 'ps')
                    else:
                        print("saveall need to be set true to continue")
                        sys.exit(0)
                    for i in range(len(self.baths)):
                        self.baths[i].noise = ReadNetCDFVar(fn, 'noise'+str(i))

                elif(ipie+1 == self.npie):
                    print("finished run")
                    if self.savep:
                        self.power = ReadNetCDFVar(fn, 'power')
                        if self.atomlist is not None:
                            self.poweratomlist = ReadNetCDFVar(
                                fn, 'poweratomlist')
                    self.t = ReadNetCDFVar(fn, 't')[0]
                    continue
                else:
                    print("ipie error")
                    sys.exit()
            else:
                print("new run")
                ipie = -1  # yes,-1
                if os.path.isfile(fnm):
                    print("reading history from previous run")
                    self.p = ReadNetCDFVar(fnm, 'p')
                    self.q = ReadNetCDFVar(fnm, 'q')
                    qhis0 = ReadNetCDFVar(fnm, 'qhis')
                    phis0 = ReadNetCDFVar(fnm, 'phis')
                    self.t = ReadNetCDFVar(fnm, 't')[0]
                    if qhis0.shape == self.qhis.shape and\
                            phis0.shape == self.phis.shape:
                        self.qhis = qhis0
                        self.phis = phis0
                elif j == 0:
                    print("initialize a new simulation")
                else:
                    print("no previous nc file exists")
                    sys.exit()
                # noise generation
                for i in range(len(self.baths)):  # loop over baths
                    self.baths[i].gnoi()
                    # print N.shape(self.baths[i].noise)
                    # stppp

                # reset qs and ps to zero
                self.ResetSavepq()

            # loop over md steps
            ipie1 = ipie+1
            iss = ipie1+N.array(range(self.npie-ipie1))
            trajfile = open('trajectories'+"."+str(self.T) +
                            "."+"run"+str(j)+'.ani', 'w')
            for i in iss:
                print("Progress of MD")
                for _ in tqdm(range(int(self.nmd/self.npie)), unit="steps", mininterval=1):
                    self.vv(j)
                    if self.nstep is not None and ((self.t-1) == 0 or (self.t-1) % self.nstep == 0):
                        #head = str(len(self.els))+'\n'+str(self.t-1)+'\n'
                        # N.savetxt(trajfile, N.column_stack((
                        #    self.els, N.transpose(N.reshape(self.xyz+self.conv*self.q,(3,len(self.els))))[:], N.transpose(N.reshape(self.f,(3,len(self.els))))[:])), header=head)
                        trajfile.write(str(len(
                            self.els))+'\n'+str(self.t-1)+'\n')
                        structure = self.xyz+self.conv*self.q
                        for ip in range(len(self.els)):
                            trajfile.write(str(self.els[ip])+'    '+str(structure[ip*3])+'   '+str(structure[ip*3+1])+'   '+str(
                                structure[ip*3+2])+'   '+str(self.f[ip*3])+'   '+str(self.f[ip*3+1])+'   '+str(self.f[ip*3+2])+'\n')
                self.dump(i, j)
            trajfile.close()

            if self.savep:
                # power spectrum
                power = N.copy(self.power)
                if self.atomlist is not None:
                    poweratomlist = [None]*len(self.atomlist)
                    for layers in range(len(self.atomlist)):
                        poweratomlist[layers] = N.copy(
                            self.poweratomlist[layers])
                self.GetPower()
                self.power = (power*(j-self.nstart)+self.power) / \
                    float(j-self.nstart+1)
                if self.atomlist is not None:
                    for layers in range(len(self.atomlist)):
                        self.poweratomlist[layers] = (
                            poweratomlist[layers]*(j-self.nstart)+self.poweratomlist[layers])/float(j-self.nstart+1)
                # power spectrum
                f = open("power."+str(self.T)+"."+"run"+str(j)+".dat", "w")
                # f.write("#k-point averaged transmission and DoS from MAMA.py\n")
                # f.write("#energy    transmission    DoSL    DoSR\n")
                for i in range(len(self.power)):
                    # only write out power spectrum upto 1.5max(hw)
                    if self.hw is not None:
                        if(self.power[i, 0] < 1.5*max(self.hw)):
                            f.write("%f     %f \n" %
                                    (self.power[i, 0], self.power[i, 1]))
                        else:
                            break
                    else:
                        f.write("%f     %f \n" %
                                (self.power[i, 0], self.power[i, 1]))
                f.close()
                if self.atomlist is not None:
                    for layers in range(len(self.atomlist)):
                        # power spectrum atomlist
                        f = open("poweratomlist."+str(layers)+"." +
                                 str(self.T)+"."+"run"+str(j)+".dat", "w")
                        # f.write("#k-point averaged transmission and DoS from MAMA.py\n")
                        # f.write("#energy    transmission    DoSL    DoSR\n")
                        for i in range(len(self.poweratomlist[layers])):
                            # only write out power spectrum upto 1.5max(hw)
                            if self.hw is not None:
                                if(self.poweratomlist[layers][i, 0] < 1.5*max(self.hw)):
                                    f.write("%f     %f \n" % (
                                        self.poweratomlist[layers][i, 0], self.poweratomlist[layers][i, 1]))
                                else:
                                    break
                            else:
                                f.write("%f     %f \n" % (
                                    self.poweratomlist[layers][i, 0], self.poweratomlist[layers][i, 1]))
                        f.close()
            # dump again, to make sure power is all right
            self.dump(i, j)

            # heat current
            for ii in range(len(self.baths)):
                fk = open("kappa."+str(self.T)+"."+"bath" +
                          str(ii)+".run"+str(j)+".dat", "w")
                # write average current
                fk.write("%i %f    %f \n" %
                         (j, self.T, N.mean(self.baths[ii].cur)*U.curcof))
                fk.close()
            if self.saveq:
                # save average structure
                f = open("avestructure."+str(self.T) +
                         "."+"run"+str(j)+".dat", "w")
                avestructure = self.conv*(self.qs.mean(axis=0)) + self.xyz
                f.write(str(len(
                    self.els))+'\n'+"average structure"+'\n')
                for ip in range(len(self.els)):
                    f.write(str(self.els[ip])+'    '+str(avestructure[ip*3])+'   '+str(
                        avestructure[ip*3+1])+'   '+str(avestructure[ip*3+2])+'\n')
                f.close()

            if self.rmnc:
                if os.path.exists("MD"+str(j-1)+".nc"):
                    print("Remove MD"+str(j-1)+".nc")
                    os.remove("MD"+str(j-1)+".nc")
                # else:
                #    print("No NC file exists.")

    def dump(self, ipie, id):
        """
        dump md information
        """
        outfile = "MD"+str(id)+".nc"
        NCfile = Dataset(outfile, 'w', 'Created '+time.ctime(time.time()))
        NCfile.title = 'Output from md.py'
        NCfile.createDimension('nph', self.nph)
        # NCfile.createDimension('na',self.na)
        NCfile.createDimension('one', 1)
        NCfile.createDimension('two', 2)
        NCfile.createDimension('mem', self.ml)
        NCfile.createDimension('nmd', self.nmd)
        NCfile.createDimension('nnmd', None)
        if self.atomlist is not None:
            NCfile.createDimension('atomlist', len(self.atomlist))

        for i in range(len(self.baths)):
            NCfile.createDimension('n'+str(i), self.baths[i].nc)

        # els
        # Write2NetCDFFile(NCfile,self.els,'elements',('na',),units='')

        if self.saveall and self.p and self.q:
            # noise series
            for i in range(len(self.baths)):
                # NCfile.createDimension('n'+str(i),self.baths[i].nc)
                Write2NetCDFFile(NCfile, self.baths[i].noise, 'noise'+str(i),
                                 ('nnmd', 'n'+str(i),), units='')
                Write2NetCDFFile(
                    NCfile, self.fhis[i], 'fhis'+str(i), ('nnmd', 'nph',), units='')
            # save all the histories of p,q or not
            Write2NetCDFFile(NCfile, self.ps, 'ps', ('nnmd', 'nph',), units='')
            Write2NetCDFFile(NCfile, self.qs, 'qs', ('nnmd', 'nph',), units='')
        if self.savep:
            # power spectrum
            Write2NetCDFFile(NCfile, self.power, 'power',
                             ('nnmd', 'two',), units='')
            Write2NetCDFFile(NCfile, self.poweratomlist,
                             'poweratomlist', ('atomlist', 'nnmd', 'two'), units='')

        # energy
        Write2NetCDFFile(NCfile, self.etot, 'energy', ('nnmd',), units='')
        # velocity
        Write2NetCDFFile(NCfile, self.p, 'p', ('nph',), units='')
        # displacement
        Write2NetCDFFile(NCfile, self.q, 'q', ('nph',), units='')
        # current time
        Write2NetCDFFile(NCfile, self.t, 't', ('one',), units='')
        # current segment
        Write2NetCDFFile(NCfile, ipie, 'ipie', ('one',), units='')
        # memory kernel for p
        Write2NetCDFFile(NCfile, self.phis, 'phis', ('mem', 'nph',), units='')
        # memory kernel for q
        Write2NetCDFFile(NCfile, self.qhis, 'qhis', ('mem', 'nph',), units='')

        NCfile.close()


# misc driver routines
def Write2NetCDFFile(file, var, varLabel, dimensions, units=None, description=None):
    # print 'Write2NetCDFFile:', varLabel, dimensions
    tmp = file.createVariable(varLabel, 'd', dimensions, zlib=True)
    tmp[:] = var
    if units:
        tmp.units = units
    if description:
        tmp.description = description


def ReadNetCDFVar(file, var):
    print("ReadNetCDFFile: reading " + var)
    f = Dataset(file, 'r')
    vv = N.array(f.variables[var])
    f.close()
    return vv


def sameq(q1, q2):
    if(len(q1) != len(q2)):
        # print "sameq: qq1 and qq2 not same length"
        return False
    qq1 = N.array(q1)
    qq2 = N.array(q2)
    dif = qq1-qq2
    # if the difference is less than 10e-10
    # use the old force
    if(max(abs(dif)) < 10e-10):
        return True
    else:
        return False


def ApplyConstraint(f, constr=None):
    """
    apply constraint to the force

    constr is an array of vectors
    This subroutine zerofy the force along each vector
    """
    if constr is None:
        return f
    nf = N.array(f)*1.0
    for i in range(len(constr)):
        nf[constr[i]] = 0
    return nf


if __name__ == "__main__":
    import time
    import numpy as N
    from baths import ebath
    from tools import calHF, calTC
    from lammpsdriver import lammpsdriver
    from md import md

    lammpsinfile = [
        #"log none",
        "units metal ",
        "dimension 3 ",
        "boundary f p p",
        "atom_style full",
        "read_data structure.data",
        "pair_style rebo ",
        "pair_coeff * * CH.rebo C H",
    ]
    # temperature
    T = 300
    delta = 0.1
    nstart = 0
    nstop = 2
    # time = 0.658fs #time unit
    dt = 0.25/0.658
    # number of md steps
    nmd = 2**10
    # initialise lammps run
    lmp = lammpsdriver(infile=lammpsinfile)
    time_start = time.time()

    print("initialise md")
    fixatoms = [range(0*3, (19+1)*3), range(181*3, (200+1)*3)]

    # print(("constraint:",constraint))
    # Molecular Junction atom indices
    slist = range(70*3, (130+1)*3)
    cutslist = [range(70*3, (89+1)*3),
                range(90*3, (109+1)*3), range(110*3, (130+1)*3)]
    # atom indices that are connecting to debyge bath
    ecatsl = range(20*3, (69+1)*3)
    ecatsr = range(131*3, (180+1)*3)

    # if slist is not given, md will initialize it using xyz
    mdrun = md(dt, nmd, T, axyz=lmp.axyz,
               nstart=nstart, nstop=nstop)
    # attache lammps driver to md
    mdrun.AddPotential(lmp)
    # unit in 0.658211814201041 fs
    damp = 100/0.658211814201041

    etal = (1.0/damp)*N.identity(len(ecatsl), N.float)
    etar = (1.0/damp)*N.identity(len(ecatsr), N.float)
    # atom indices that are connecting to bath
    ebl = ebath(ecatsl, T*(1+delta/2), mdrun.dt, mdrun.nmd,
                wmax=1., nw=500, bias=0.0, efric=etal, classical=False, zpmotion=True)
    mdrun.AddBath(ebl)
    ebr = ebath(ecatsr, T*(1-delta/2), mdrun.dt, mdrun.nmd,
                wmax=1., nw=500, bias=0.0, efric=etar, classical=False, zpmotion=True)
    mdrun.AddBath(ebr)

    mdrun.AddConstr(fixatoms)
    # mdrun.CalPowerSpec()
    # mdrun.AddPowerSection([ecatsl, slist, ecatsr])
    # mdrun.CalAveStruct()
    mdrun.SaveTraj()
    # mdrun.RemoveNC()

    mdrun.Run()
    lmp.quit()
    calHF()
    calTC(delta=delta, dlist=0)
    time_end = time.time()
    print('time cost', time_end-time_start, 's.')
