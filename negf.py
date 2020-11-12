#!/usr/bin/env python
# -*- coding: utf-8 -*
import sys
import numpy as np
from lammps import lammps


class bpt:
    # Use NEGF to calculate ballistic phonon transport
    def __init__(self, infile, maxomega, damp, dofatomofbath, dofatomfixed=[[], []], num=1000, vector=False):
        print('Class init')
        # reduced Planck constant unit in: eV*ps
        self.rpc = 6.582119569e-4
        # Boltzmann constant unit in: eV/K
        self.bc = 8.617333262e-5
        self.infile = infile
        self.damp = damp
        self.maxomega = maxomega/self.rpc
        self.intnum = num
        self.dofatomfixed = dofatomfixed
        self.dofatomofbath = dofatomofbath
        self.getdynmat()
        self.gettm(vector)

    def getdynmat(self):
        #import os
        lmp = lammps()
        #lmp = lammps(cmdargs=['-screen', 'none', '-log', 'none'])
        # if os.path.exists('atoms.xyz'):
        #    print('Remove atoms.xyz')
        #    os.remove('atoms.xyz')
        print('LAMMPS init')
        lmp.commands_list(self.infile)
        print('Calculate dynamical matrix')
        lmp.command('dynamical_matrix all eskm 0.000001 file dynmat.dat')
        #lmp.command('dump 1 all xyz 1 atoms.xyz')
        #lmp.command('run 0')
        self.natoms = lmp.get_natoms()
        lmp.close()
        print('Calculate angular frequency')

        self.dynmat = []
        self.omegas = []
        self.doffreeatom = 0
        dynmatdat = np.loadtxt('dynmat.dat')
        dynlen = int(3*np.sqrt(len(dynmatdat)/3))
        self.dynmat = dynmatdat.reshape((dynlen, dynlen))
        self.dynmat = np.delete(self.dynmat, self.dofatomfixed[0], axis=0)
        self.dynmat = np.delete(self.dynmat, self.dofatomfixed[0], axis=1)
        self.dynmat = np.delete(self.dynmat, [
                                dof-len(self.dofatomfixed[0]) for dof in self.dofatomfixed[1]], axis=0)
        self.dynmat = np.delete(self.dynmat, [
                                dof-len(self.dofatomfixed[0]) for dof in self.dofatomfixed[1]], axis=1)
        #eigvals, eigvecs = np.linalg.eig(self.dynmat)
        eigvals = np.linalg.eigh(self.dynmat)[0]
        # frequencies in eV
        #self.omegas = np.sqrt(np.abs(eigvals))*self.rpc
        # Or remove false frequency
        reigvals = [item for item in eigvals if item >= 0]
        ieigvals = [item for item in eigvals if item < 0]
        self.omegas = np.sqrt(np.abs(reigvals))*self.rpc
        np.savetxt('omegas.dat', self.omegas)
        if len(ieigvals) != 0:
            print('False frequency exists in system')
            np.savetxt('iomegas.dat', np.sqrt(np.abs(ieigvals))*self.rpc,
                       header='Frequency obtained by taking the absolute value of negative eigenvalue')
        # print(len(reigvals),'>=0',len(ieigvals),'<0')

    def gettm(self, vector):
        print('Calculate transmission')
        x = np.linspace(0, self.maxomega, self.intnum+1)
        if vector:
            function = np.vectorize(self.tm)
            self.tmnumber = np.array(
                np.column_stack((x, np.array(function(x)))))
        else:
            from tqdm import tqdm
            tm = []
            for var in tqdm(x, unit="steps", mininterval=1):
                tm.append(self.tm(var))
            self.tmnumber = np.array(np.column_stack((x, np.array(tm))))
        np.savetxt('transmission.dat', np.column_stack(
            (self.tmnumber[:, 0]*self.rpc, self.tmnumber[:, 1])))
        print('Transmission saved')

    def selfenergy(self, omega, dofatoms):
        return -1j*omega*(1/self.damp)*self.atomofbath(dofatoms)

    def atomofbath(self, dofatoms):
        semat = np.zeros((self.natoms*3, self.natoms*3))
        for dofatom in dofatoms:
            semat[dofatom][dofatom] = 1
        semat = np.delete(semat, self.dofatomfixed[0], axis=0)
        semat = np.delete(semat, self.dofatomfixed[0], axis=1)
        semat = np.delete(semat, [dof-len(self.dofatomfixed[0])
                                  for dof in self.dofatomfixed[1]], axis=0)
        semat = np.delete(semat, [dof-len(self.dofatomfixed[0])
                                  for dof in self.dofatomfixed[1]], axis=1)
        if len(semat) != len(self.dynmat) or self.natoms*3 != len(self.dofatomfixed[0]) + len(self.dofatomfixed[1]) + len(semat):
            print('System DOF test failed, check again')
            sys.exit()
        return semat

    def retargf(self, omega):
        # retarded Green function
        return np.linalg.inv(omega*omega*np.identity(len(self.dynmat))-self.dynmat-self.selfenergy(omega, self.dofatomofbath[0])-self.selfenergy(omega, self.dofatomofbath[1]))

    def gamma(self, Pi):
        return -1j*(Pi-Pi.conjugate().transpose())

    def bosedist(self, omega, T):
        # Bose Einstein distribution
        if abs(T) < 1e-30:
            # print('T %e is too small. Set kBT min.' % T)
            return 1/(np.exp(self.rpc*omega*np.iinfo(np.int32).max)-1)
        elif abs(omega/T) < 1e-30:
            # print('omega %e is too small. Set bose einstein distribution max.' % omega)
            return np.iinfo(np.int32).max
        else:
            return 1/(np.exp(self.rpc*omega/self.bc/T)-1)

    def tm(self, omega):
        # Transmission
        return np.real(np.trace(np.dot(np.dot(np.dot(self.retargf(omega), self.gamma(self.selfenergy(
            omega, self.dofatomofbath[0]))), self.retargf(omega).conjugate().transpose()), self.gamma(self.selfenergy(omega, self.dofatomofbath[1])))))

    def thermalcurrent(self, T, delta):
        # def f(omega):
        #    return self.rpc*omega/2 / \
        #        np.pi*self.tm(omega)*(self.bosedist(omega, T*(1+0.5*delta)) -
        #                              self.bosedist(omega, T*(1-0.5*delta)))

        # def trape(function, maxnumber, n):
        #    function = np.vectorize(function)
        #    arr = function(np.linspace(0, maxnumber, n+1))
        #    return (float(maxnumber - 0)/n/2.)*(2*arr.sum() - arr[0] - arr[-1])

        def f(i):
            return self.rpc*self.tmnumber[i, 0]/2 / \
                np.pi*self.tmnumber[i, 1]*(self.bosedist(self.tmnumber[i, 0], T*(1+0.5*delta)) -
                                           self.bosedist(self.tmnumber[i, 0], T*(1-0.5*delta)))

        def trape(function):
            n = len(self.tmnumber[:, 0]) - 1
            if n != self.intnum:
                print('Error in number of omega')
                sys.exit()
            function = np.vectorize(function)
            arr = function(range(n+1))
            return (float(self.tmnumber[-1, 0] - self.tmnumber[0, 0])/n/2.)*(2*arr.sum() - arr[0] - arr[-1])
        # Unit in nW
        # return trape(f, self.maxomega, self.intnum)*1.60217662*1e2
        return trape(f)*1.60217662*1e2

    def thermalconductance(self, T, delta):
        return self.thermalcurrent(T, delta)/(T*delta)

    def plotresult(self, lines=180):
        from matplotlib import pyplot as plt
        plt.figure(0)
        plt.hist(self.omegas, bins=lines)
        plt.xlabel('Frequence(eV)')
        plt.ylabel('Number')
        #plt.xlim(0, self.maxomega*self.rpc)
        plt.savefig('omegas.png')
        plt.figure(1)
        plt.plot(self.tmnumber[:, 0]*self.rpc, self.tmnumber[:, 1])
        plt.xlabel('Frequence(eV)')
        plt.ylabel('Transmission')
        plt.savefig('transmission.png')


if __name__ == '__main__':
    '''
    Units
    Time: ps
    Frequence: eV
    Temperture: K
    Heat Current: nW
    '''
    import time
    import numpy as np
    from negf import bpt
    from matplotlib import pyplot as plt
    infile = ['atom_style full',
              'units metal',
              'boundary f p p',
              'read_data structure.data',
              'pair_style rebo',
              'pair_coeff * * CH.rebo C H',
              ]
    time_start = time.time()
    atomfixed = [range(0*3, (19+1)*3), range(181*3, (200+1)*3)]
    atomofbath = [range(20*3, (69+1)*3), range(131*3, (180+1)*3)]
    mybpt = bpt(infile, 0.25, 0.1, atomofbath, atomfixed, 100)
    mybpt.plotresult()
    # T_H/C = T*(1±delta/2)
    T = [100, 200, 300, 400, 500, 600, 700,
         800, 900, 1000]
    delta = 0.1
    thermalconductance = []
    for temp in T:
        thermalconductance.append(
            [temp, mybpt.thermalconductance(temp, delta)])
    np.savetxt('thermalconductance.dat', thermalconductance)
    plt.figure(5)
    plt.plot(np.array(thermalconductance)[
        :, 0], np.array(thermalconductance)[:, 1])
    plt.xlabel('Temperature(K)')
    plt.ylabel('Thermal Conductance(nW/K)')
    plt.savefig('thermalconductance.png')
    time_end = time.time()
    print('time cost', time_end-time_start, 's')
