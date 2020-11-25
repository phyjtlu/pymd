# LAMMPS driver for Langevin molecular dynamics
# Adapted from LAMMPS python wraper

import ctypes
import sys

import numpy as N
from lammps import lammps

import units as U

args = "-screen none"


def get_atomname(mass):
    """
    get the element name from its atomic mass by checking the dictionary
    """
    for key, value in list(U.AtomicMassTable.items()):
        if abs(mass-value) < 0.01:
            return key


class lammpsdriver(lammps):
    # create instance of LAMMPS
    def __init__(self, infile, label="",
                 constraints=[], tdir="./", lunit="Ang", eunit="eV", md2ang=0.06466,
                 name="", cmdargs=args.split(), ptr=None, comm=None
                 ):
        lammps.__init__(self, name="", cmdargs=cmdargs, ptr=None, comm=None)
        self.infile = infile
        self.md2ang = md2ang
        self.constraints = constraints
        self.label = label
        self.lunit = lunit
        self.eunit = eunit
        if self.eunit == "eV":
            self.para = 1.0
        elif self.eunit == "Kcal/mole":
            self.para = 0.043
        else:
            print("Wrong vaule in eunit")
            sys.exit(0)
        # start lammps
        self.start()

    def start(self, np=1):
        print("LAMMPS launched")
        # todo:better to set the unit to metals here again
        #self.command("units metal")

        #lines = open(self.infile, 'r').readlines()
        #for line in lines: self.command(line)
        self.commands_list(self.infile)
        self.type = N.array(self.gather_atoms("type", 0, 1))
        #self.mass = N.array(self.gather_atoms("mass",1,1))
        self.mass = self.extract_atom("mass", 2)
        self.number = self.get_natoms()
        self.els = []
        for type in self.type:
            self.els.append(self.mass[type])
        self.xyz = self.gather_atoms("x", 1, 3)
        self.conv = self.md2ang*N.array([3*[1.0/N.sqrt(mass)]
                                         for mass in self.els]).flatten()
        self.type = N.array(self.gather_atoms("type", 0, 1))
        #self.mass = N.array(self.gather_atoms("mass",1,1))
        self.axyz = []
        for i, a in enumerate(self.els):
            self.axyz.append([get_atomname(a), self.xyz[i*3],
                              self.xyz[i*3+1], self.xyz[i*3+2]])
        # print(self.conv.shape)
        self.initforce()

    def quit(self):
        self.close()
        print("Mission completed.")

    def newx(self, q):
        newxyz = self.xyz + self.conv*q
        return newxyz.ctypes

    def absforce(self, q):
        self.scatter_atoms("x", 1, 3, self.newx(q))
        self.command("run 0")
        return self.para*self.conv*N.array(self.gather_atoms("f", 1, 3))

    def initforce(self):
        print("Calculate zero displacement force")
        self.f0 = self.absforce(N.zeros(3*self.number))

    def force(self, q):
        return self.absforce(q) - self.f0

    def energy(self, eargs="pe"):  # energy,eargs:"pe","ke" or "etotal".
        return self.get_thermo(eargs)
