#LAMMPS driver for Langevin molecular dynamics
#Adapted from LAMMPS python wraper
#By Li Gen
#
# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   http://lammps.sandia.gov, Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

# Python wrappers on LAMMPS library via ctypes

# for python3 compatibility

from __future__ import print_function

# imports for simple LAMMPS python wrapper module "lammps"

# import sys,traceback,types
from ctypes import *
from os.path import dirname, abspath, join
from inspect import getsourcefile

# imports for advanced LAMMPS python wrapper modules "PyLammps" and "IPyLammps"

# from collections import namedtuple
# import os
# import select
# import re
import sys
import numpy as N

#<<<<<<< master
#=======
#atomic table and periodic table from units.py
#from units import *
PeriodicTable = {'H':1,1:'H','D':1001,1001:'D','He':2,2:'He','Li':3,3:'Li','Be':4,4:'Be','B':5,5:'B','C':6,6:'C','N':7,7:'N','O':8,8:'O','F':9,9:'F','Ne':10,10:'Ne','Na':11,11:'Na','Mg':12,12:'Mg','Al':13,13:'Al','Si':14,14:'Si','P':15,15:'P','S':16,16:'S','Cl':17,17:'Cl','Ar':18,18:'Ar','K':19,19:'K','Ca':20,20:'Ca','Sc':21,21:'Sc','Ti':22,22:'Ti','V':23,23:'V','Cr':24,24:'Cr','Mn':25,25:'Mn','Fe':26,26:'Fe','Co':27,27:'Co','Ni':28,28:'Ni','Cu':29,29:'Cu','Zn':30,30:'Zn','Ga':31,31:'Ga','Ge':32,32:'Ge','As':33,33:'As','Se':34,34:'Se','Br':35,35:'Br','Kr':36,36:'Kr','Rb':37,37:'Rb','Sr':38,38:'Sr','Y':39,39:'Y','Zr':40,40:'Zr','Nb':41,41:'Nb','Mo':42,42:'Mo','Tc':43,43:'Tc','Ru':44,44:'Ru','Rh':45,45:'Rh','Pd':46,46:'Pd','Ag':47,47:'Ag','Cd':48,48:'Cd','In':49,49:'In','Sn':50,50:'Sn','Sb':51,51:'Sb','Te':52,52:'Te','I':53,53:'I','Xe':54,54:'Xe','Cs':55,55:'Cs','Ba':56,56:'Ba','La':57,57:'La','Ce':58,58:'Ce','Pr':59,59:'Pr','Nd':60,60:'Nd','Pm':61,61:'Pm','Sm':62,62:'Sm','Eu':63,63:'Eu','Gd':64,64:'Gd','Tb':65,65:'Tb','Dy':66,66:'Dy','Ho':67,67:'Ho','Er':68,68:'Er','Tm':69,69:'Tm','Yb':70,70:'Yb','Lu':71,71:'Lu','Hf':72,72:'Hf','Ta':73,73:'Ta','W':74,74:'W','Re':75,75:'Re','Os':76,76:'Os','Ir':77,77:'Ir','Pt':78,78:'Pt','Au':79,79:'Au','Hg':80,80:'Hg','Tl':81,81:'Tl','Pb':82,82:'Pb','Bi':83,83:'Bi','Po':84,84:'Po','At':85,85:'At','Rn':86,86:'Rn','Fr':87,87:'Fr','Ra':88,88:'Ra','Ac':89,89:'Ac','Th':90,90:'Th','Pa':91,91:'Pa','U':92,92:'U','Np':93,93:'Np','Pu':94,94:'Pu','Am':95,95:'Am','Cm':96,96:'Cm','Bk':97,97:'Bk','Cf':98,98:'Cf','Es':99,99:'Es','Fm':100,100:'Fm','Md':101,101:'Md','No':102,102:'No'}

AtomicMassTable={'H':1.00794, 'He':4.002602, 'Li':6.941, 'Be':9.012182, \
    'B':10.811, 'C':12.0107, 'N':14.0067, 'O':15.9994, \
    'F':18.9984032, 'Ne':20.1791, 'Na':22.98976928, 'Mg':24.3050, \
    'Al':26.9815386, 'Si':28.0855, 'P':30.973762, 'S':32.065, \
    'Cl':35.453, 'Ar':39.948, 'K':39.0983, 'Ca':40.078, \
    'Sc':44.955912, 'Ti':47.867, 'V':50.9415, 'Cr':51.9961, \
    'Mn':54.938045, 'Fe':55.845, 'Co':58.933195, 'Ni':58.6934, \
    'Cu':63.546, 'Zn':65.38, 'Ga':69.723, 'Ge':72.64, \
    'As':74.92160, 'Se':78.96, 'Br':79.904, 'Kr':83.798, \
    'Rb':85.4678, 'Sr':87.62, 'Y':88.90585, 'Zr':91.224, \
    'Nb':92.90638, 'Mo':95.96, 'Tc':98, 'Ru':101.07, \
    'Rh':102.90550, 'Pd':106.42, 'Ag':107.8682, 'Cd':112.411, \
    'In':114.818, 'Sn':118.710, 'Sb':121.760, 'Te':127.60, \
    'I':126.90447, 'Xe':131.293, 'Cs':132.9054519, 'Ba':137.327, \
    'La':138.90547, 'Ce':140.116, 'Pr':140.90765, 'Nd':144.242, \
    'Pm':145, 'Sm':150.36, 'Eu':151.964, 'Gd':157.25, \
    'Tb':158.92535, 'Dy':162.500, 'Ho':164.93032, 'Er':167.259, \
    'Tm':168.93421, 'Yb':173.054, 'Lu':174.9668, 'Hf':178.49, \
    'Ta':180.94788, 'W':183.84, 'Re':186.207, 'Os':190.23, \
    'Ir':192.217, 'Pt':195.084, 'Au':196.966569, 'Hg':200.59, \
    'Tl':204.3833, 'Pb':207.2, 'Bi':208.98040, 'Po':209, \
    'At':210, 'Rn':222, 'Fr':223, 'Ra':226, 'Ac':227, \
    'Th':232.03806, 'Pa':231.03586, 'U':238.02891, 'Np':237, \
    'Pu':244, 'Am':243, 'Cm':247, 'Bk':247, 'Cf':251, \
    'Es':252, 'Fm':257, 'Md':258, 'No':259, 'Lr':262, \
    'Rf':265, 'Db':268, 'Sg':271, 'Bh':272, 'Hs':270, \
    'Mt':276, 'Ds':281, 'Rg':280, 'Cn':285, 'Uut':284, \
    'Uuq':289, 'Uup':288, 'Uuh':293, 'Uus':294, 'Uuo':294}

def get_atomname(mass):
    """
    get the element name from its atomic mass by checking the dictionary
    """
    for key, value in AtomicMassTable.items():
        if abs(mass-value) < 0.01:
            return key

#>>>>>>> master
def get_ctypes_int(size):
    if size == 4:
        return c_int32
    elif size == 8:
        return c_int64
    return c_int


class MPIAbortException(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


class lammps(object):
    # detect if Python is using version of mpi4py that can pass a communicator

    has_mpi4py = False
    try:
        from mpi4py import MPI
        from mpi4py import __version__ as mpi4py_version
        if mpi4py_version.split('.')[0] in ['2', '3']: has_mpi4py = True
    except:
        pass

    # create instance of LAMMPS

#<<<<<<< master
#    def __init__(self, infile, label="", mesh=100., dmtol=0.001, \
#=======
    def __init__(self, infile, label="", \
                 constraints=[], tdir="./", lunit="Ang", eunit="eV", md2ang=0.06466, \
                 name="", cmdargs=None, ptr=None, comm=None
                 ):
        self.comm = comm
        self.opened = 0

        # determine module location

        modpath = dirname(abspath(getsourcefile(lambda: 0)))
        self.lib = None

        # if a pointer to a LAMMPS object is handed in,
        # all symbols should already be available

        try:
            if ptr: self.lib = CDLL("", RTLD_GLOBAL)
        except:
            self.lib = None

        # load liblammps.so unless name is given
        #   if name = "g++", load liblammps_g++.so
        # try loading the LAMMPS shared object from the location
        #   of lammps.py with an absolute path,
        #   so that LD_LIBRARY_PATH does not need to be set for regular install
        # fall back to loading with a relative path,
        #   typically requires LD_LIBRARY_PATH to be set appropriately

        if not self.lib:
            try:
                if not name:
                    self.lib = CDLL(join(modpath, "liblammps.so"), RTLD_GLOBAL)
                else:
                    self.lib = CDLL(join(modpath, "liblammps_%s.so" % name),
                                    RTLD_GLOBAL)
            except:
                if not name:
                    self.lib = CDLL("liblammps.so", RTLD_GLOBAL)
                else:
                    self.lib = CDLL("liblammps_%s.so" % name, RTLD_GLOBAL)

        # define ctypes API for each library method
        # NOTE: should add one of these for each lib function

        self.lib.lammps_extract_box.argtypes = \
            [c_void_p, POINTER(c_double), POINTER(c_double),
             POINTER(c_double), POINTER(c_double), POINTER(c_double),
             POINTER(c_int), POINTER(c_int)]
        self.lib.lammps_extract_box.restype = None

        self.lib.lammps_reset_box.argtypes = \
            [c_void_p, POINTER(c_double), POINTER(c_double), c_double, c_double, c_double]
        self.lib.lammps_reset_box.restype = None

        self.lib.lammps_gather_atoms.argtypes = \
            [c_void_p, c_char_p, c_int, c_int, c_void_p]
        self.lib.lammps_gather_atoms.restype = None

        self.lib.lammps_gather_atoms_concat.argtypes = \
            [c_void_p, c_char_p, c_int, c_int, c_void_p]
        self.lib.lammps_gather_atoms_concat.restype = None

        self.lib.lammps_gather_atoms_subset.argtypes = \
            [c_void_p, c_char_p, c_int, c_int, c_int, POINTER(c_int), c_void_p]
        self.lib.lammps_gather_atoms_subset.restype = None

        self.lib.lammps_scatter_atoms.argtypes = \
            [c_void_p, c_char_p, c_int, c_int, c_void_p]
        self.lib.lammps_scatter_atoms.restype = None

        self.lib.lammps_scatter_atoms_subset.argtypes = \
            [c_void_p, c_char_p, c_int, c_int, c_int, POINTER(c_int), c_void_p]
        self.lib.lammps_scatter_atoms_subset.restype = None

        # if no ptr provided, create an instance of LAMMPS
        #   don't know how to pass an MPI communicator from PyPar
        #   but we can pass an MPI communicator from mpi4py v2.0.0 and later
        #   no_mpi call lets LAMMPS use MPI_COMM_WORLD
        #   cargs = array of C strings from args
        # if ptr, then are embedding Python in LAMMPS input script
        #   ptr is the desired instance of LAMMPS
        #   just convert it to ctypes ptr and store in self.lmp

        if not ptr:

            # with mpi4py v2, can pass MPI communicator to LAMMPS
            # need to adjust for type of MPI communicator object
            # allow for int (like MPICH) or void* (like OpenMPI)

            if comm:
                if not lammps.has_mpi4py:
                    raise Exception('Python mpi4py version is not 2 or 3')
                if lammps.MPI._sizeof(lammps.MPI.Comm) == sizeof(c_int):
                    MPI_Comm = c_int
                else:
                    MPI_Comm = c_void_p

                narg = 0
                cargs = 0
                if cmdargs:
                    cmdargs.insert(0, "lammps.py")
                    narg = len(cmdargs)
                    for i in range(narg):
                        if type(cmdargs[i]) is str:
                            cmdargs[i] = cmdargs[i].encode()
                    cargs = (c_char_p * narg)(*cmdargs)
                    self.lib.lammps_open.argtypes = [c_int, c_char_p * narg, \
                                                     MPI_Comm, c_void_p()]
                else:
                    self.lib.lammps_open.argtypes = [c_int, c_int, \
                                                     MPI_Comm, c_void_p()]

                self.lib.lammps_open.restype = None
                self.opened = 1
                self.lmp = c_void_p()
                comm_ptr = lammps.MPI._addressof(comm)
                comm_val = MPI_Comm.from_address(comm_ptr)
                self.lib.lammps_open(narg, cargs, comm_val, byref(self.lmp))

            else:
                if lammps.has_mpi4py:
                    from mpi4py import MPI
                    self.comm = MPI.COMM_WORLD
                self.opened = 1
                if cmdargs:
                    cmdargs.insert(0, "lammps.py")
                    narg = len(cmdargs)
                    for i in range(narg):
                        if type(cmdargs[i]) is str:
                            cmdargs[i] = cmdargs[i].encode()
                    cargs = (c_char_p * narg)(*cmdargs)
                    self.lmp = c_void_p()
                    self.lib.lammps_open_no_mpi(narg, cargs, byref(self.lmp))
                else:
                    self.lmp = c_void_p()
                    self.lib.lammps_open_no_mpi(0, None, byref(self.lmp))
                    # could use just this if LAMMPS lib interface supported it
                    # self.lmp = self.lib.lammps_open_no_mpi(0,None)

        else:
            # magic to convert ptr to ctypes ptr
            if sys.version_info >= (3, 0):
                # Python 3 (uses PyCapsule API)
                pythonapi.PyCapsule_GetPointer.restype = c_void_p
                pythonapi.PyCapsule_GetPointer.argtypes = [py_object, c_char_p]
                self.lmp = c_void_p(pythonapi.PyCapsule_GetPointer(ptr, None))
            else:
                # Python 2 (uses PyCObject API)
                pythonapi.PyCObject_AsVoidPtr.restype = c_void_p
                pythonapi.PyCObject_AsVoidPtr.argtypes = [py_object]
                self.lmp = c_void_p(pythonapi.PyCObject_AsVoidPtr(ptr))

        # optional numpy support (lazy loading)
        self._numpy = None 
        # set default types
        self.c_bigint = get_ctypes_int(self.extract_setting("bigint"))
        self.c_tagint = get_ctypes_int(self.extract_setting("tagint"))
        self.c_imageint = get_ctypes_int(self.extract_setting("imageint"))
        self._installed_packages = None
        self.infile = infile
        self.md2ang = md2ang
        self.constraints = constraints
        self.label = label
        self.lunit = lunit
        self.eunit = eunit


        #start lammps
        self.start()

    def start(self, np=1):
        print("lammps launched")
        #todo:better to set the unit to metals here again
        self.command("units metal")

        lines = open(self.infile, 'r').readlines()
        for line in lines: self.command(line)
        self.type = N.array(self.gather_atoms("type", 0, 1))
        #self.mass = N.array(self.gather_atoms("mass",1,1))
        self.mass = self.extract_atom("mass",2)
        self.els = []
        for type in self.type:
            self.els.append(self.mass[type])
        self.xyz = self.gather_atoms("x", 1, 3)
        self.newxyz = self.gather_atoms("x", 1, 3)
        self.conv = self.md2ang*N.array([3*[1.0/N.sqrt(mass)]
                                         for mass in self.els]).flatten()
        self.number = self.get_natoms()
        #<<<<<<< Updated upstream
#<<<<<<< master
#=======

#        #conversion factor from eV/Ang (force from lammps)
#        #to the intermal unit of MD
#        #todo
#        #self.els = self.extract_atom("mass",2)
#        self.els = self.gather_atoms("mass",1,1)
#        print("self.els:",self.els[1])
#        #self.conv = self.md2ang*N.array([3*[1.0/N.sqrt(AtomicMassTable[a])]\
#                   #                        for a in self.els]).flatten()
#        self.conv = 1.
#>>>>>>> master
#=======

        #conversion factor from eV/Ang (force from lammps)
        #to the intermal unit of MD
        #todo
        #self.els = self.extract_atom("mass",2)
        #self.els = self.gather_atoms("mass",1,1)
        #print("self.els:",self.els[1])
        #self.conv = self.md2ang*N.array([3*[1.0/N.sqrt(AtomicMassTable[a])]\
                   #                        for a in self.els]).flatten()
        #self.conv = 1.


        self.type = N.array(self.gather_atoms("type", 0, 1))
        #self.mass = N.array(self.gather_atoms("mass",1,1))
        self.mass = self.extract_atom("mass",2)
        self.els = []
        for type in self.type:
            self.els.append(self.mass[type])
        self.xyz = self.gather_atoms("x", 1, 3)
        self.newxyz = self.gather_atoms("x", 1, 3)
        self.conv = self.md2ang*N.array([3*[1.0/N.sqrt(mass)]
                                         for mass in self.els]).flatten()
        self.axyz = []
        for i, a in enumerate(self.els):
            self.axyz.append([get_atomname(a), self.xyz[i*3],
                         self.xyz[i*3+1], self.xyz[i*3+2]])
        #print(self.conv.shape)
        self.initforce()

    def quit(self):
        self.close()
        print("Quit lammps!")

    def newx(self, q):
        for i in range(3*self.number):
            self.newxyz[i] = self.xyz[i] + self.conv[i] * q[i]
        return self.newxyz

    def absforce(self, q):
        self.scatter_atoms("x", 1, 3, self.newx(q))
        self.command("run 1")
        return self.conv*N.array(self.gather_atoms("f", 1, 3))
        
    def initforce(self):
        print("Calculate zero displacement force")
        extq = N.zeros(3*self.number)
        self.f0 = self.absforce(extq)

    def force(self, q):
        f = self.absforce(q) - self.f0
        return f



    #----------------------------------------------------------------------
    # send a single command
    def command(self, cmd):
        if cmd: cmd = cmd.encode()
        self.lib.lammps_command(self.lmp, cmd)

        if self.has_exceptions and self.lib.lammps_has_error(self.lmp):
            sb = create_string_buffer(100)
            error_type = self.lib.lammps_get_last_error_message(self.lmp, sb, 100)
            error_msg = sb.value.decode().strip()

            if error_type == 2:
                raise MPIAbortException(error_msg)
            raise Exception(error_msg)

    # send a list of commands
    def commands_list(self, cmdlist):
        cmds = [x.encode() for x in cmdlist if type(x) is str]
        args = (c_char_p * len(cmdlist))(*cmds)
        self.lib.lammps_commands_list(self.lmp, len(cmdlist), args)

    # send a string of commands
    def commands_string(self, multicmd):
        if type(multicmd) is str: multicmd = multicmd.encode()
        self.lib.lammps_commands_string(self.lmp, c_char_p(multicmd))

    # extract lammps type byte sizes
    def extract_setting(self, name):
        if name: name = name.encode()
        self.lib.lammps_extract_atom.restype = c_int
        return int(self.lib.lammps_extract_setting(self.lmp, name))

    # extract global info
    def extract_global(self, name, type):
        if name: name = name.encode()
        if type == 0:
            self.lib.lammps_extract_global.restype = POINTER(c_int)
        elif type == 1:
            self.lib.lammps_extract_global.restype = POINTER(c_double)
        else:
            return None
        ptr = self.lib.lammps_extract_global(self.lmp, name)
        return ptr[0]

    # extract global info
    def extract_box(self):
        boxlo = (3 * c_double)()
        boxhi = (3 * c_double)()
        xy = c_double()
        yz = c_double()
        xz = c_double()
        periodicity = (3 * c_int)()
        box_change = c_int()

        self.lib.lammps_extract_box(self.lmp, boxlo, boxhi,
                                    byref(xy), byref(yz), byref(xz),
                                    periodicity, byref(box_change))

        boxlo = boxlo[:3]
        boxhi = boxhi[:3]
        xy = xy.value
        yz = yz.value
        xz = xz.value
        periodicity = periodicity[:3]
        box_change = box_change.value

        return boxlo, boxhi, xy, yz, xz, periodicity, box_change

    # extract per-atom info
    # NOTE: need to insure are converting to/from correct Python type
    #   e.g. for Python list or NumPy or ctypes

    def extract_atom(self, name, type):
        if name: name = name.encode()
        if type == 0:
            self.lib.lammps_extract_atom.restype = POINTER(c_int)
        elif type == 1:
            self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_int))
        elif type == 2:
            self.lib.lammps_extract_atom.restype = POINTER(c_double)
        elif type == 3:
            self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_double))
        else:
            return None
        ptr = self.lib.lammps_extract_atom(self.lmp, name)
        return ptr

    # shut-down LAMMPS instance

    def __del__(self):
        if self.lmp and self.opened:
            self.lib.lammps_close(self.lmp)
            self.opened = 0

    def close(self):
        if self.opened: self.lib.lammps_close(self.lmp)
        self.lmp = None
        self.opened = 0

    def version(self):
        return self.lib.lammps_version(self.lmp)

    def file(self, file):
        if file: file = file.encode()
        self.lib.lammps_file(self.lmp, file)

    @property
    def numpy(self):
        if not self._numpy:
            import numpy as np
            class LammpsNumpyWrapper:
                def __init__(self, lmp):
                    self.lmp = lmp

                def _ctype_to_numpy_int(self, ctype_int):
                    if ctype_int == c_int32:
                        return np.int32
                    elif ctype_int == c_int64:
                        return np.int64
                    return np.intc

                def extract_atom_iarray(self, name, nelem, dim=1):
                    if name in ['id', 'molecule']:
                        c_int_type = self.lmp.c_tagint
                    elif name in ['image']:
                        c_int_type = self.lmp.c_imageint
                    else:
                        c_int_type = c_int

                    np_int_type = self._ctype_to_numpy_int(c_int_type)

                    if dim == 1:
                        tmp = self.lmp.extract_atom(name, 0)
                        ptr = cast(tmp, POINTER(c_int_type * nelem))
                    else:
                        tmp = self.lmp.extract_atom(name, 1)
                        ptr = cast(tmp[0], POINTER(c_int_type * nelem * dim))

                    a = np.frombuffer(ptr.contents, dtype=np_int_type)
                    a.shape = (nelem, dim)
                    return a

                def extract_atom_darray(self, name, nelem, dim=1):
                    if dim == 1:
                        tmp = self.lmp.extract_atom(name, 2)
                        ptr = cast(tmp, POINTER(c_double * nelem))
                    else:
                        tmp = self.lmp.extract_atom(name, 3)
                        ptr = cast(tmp[0], POINTER(c_double * nelem * dim))

                    a = np.frombuffer(ptr.contents)
                    a.shape = (nelem, dim)
                    return a

            self._numpy = LammpsNumpyWrapper(self)
        return self._numpy

    # extract compute info

    def extract_compute(self, id, style, type):
        if id: id = id.encode()
        if type == 0:
            if style > 0: return None
            self.lib.lammps_extract_compute.restype = POINTER(c_double)
            ptr = self.lib.lammps_extract_compute(self.lmp, id, style, type)
            return ptr[0]
        if type == 1:
            self.lib.lammps_extract_compute.restype = POINTER(c_double)
            ptr = self.lib.lammps_extract_compute(self.lmp, id, style, type)
            return ptr
        if type == 2:
            if style == 0:
                self.lib.lammps_extract_compute.restype = POINTER(c_int)
                ptr = self.lib.lammps_extract_compute(self.lmp, id, style, type)
                return ptr[0]
            else:
                self.lib.lammps_extract_compute.restype = POINTER(POINTER(c_double))
                ptr = self.lib.lammps_extract_compute(self.lmp, id, style, type)
                return ptr
        return None

    # extract fix info
    # in case of global datum, free memory for 1 double via lammps_free()
    # double was allocated by library interface function

    def extract_fix(self, id, style, type, i=0, j=0):
        if id: id = id.encode()
        if style == 0:
            self.lib.lammps_extract_fix.restype = POINTER(c_double)
            ptr = self.lib.lammps_extract_fix(self.lmp, id, style, type, i, j)
            result = ptr[0]
            self.lib.lammps_free(ptr)
            return result
        elif (style == 1) or (style == 2):
            if type == 1:
                self.lib.lammps_extract_fix.restype = POINTER(c_double)
            elif type == 2:
                self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
            else:
                return None
            ptr = self.lib.lammps_extract_fix(self.lmp, id, style, type, i, j)
            return ptr
        else:
            return None

    # extract variable info
    # free memory for 1 double or 1 vector of doubles via lammps_free()
    # for vector, must copy nlocal returned values to local c_double vector
    # memory was allocated by library interface function

    def extract_variable(self, name, group, type):
        if name: name = name.encode()
        if group: group = group.encode()
        if type == 0:
            self.lib.lammps_extract_variable.restype = POINTER(c_double)
            ptr = self.lib.lammps_extract_variable(self.lmp, name, group)
            result = ptr[0]
            self.lib.lammps_free(ptr)
            return result
        if type == 1:
            self.lib.lammps_extract_global.restype = POINTER(c_int)
            nlocalptr = self.lib.lammps_extract_global(self.lmp, "nlocal".encode())
            nlocal = nlocalptr[0]
            result = (c_double * nlocal)()
            self.lib.lammps_extract_variable.restype = POINTER(c_double)
            ptr = self.lib.lammps_extract_variable(self.lmp, name, group)
            for i in range(nlocal): result[i] = ptr[i]
            self.lib.lammps_free(ptr)
            return result
        return None

    # return current value of thermo keyword

    def get_thermo(self, name):
        if name: name = name.encode()
        self.lib.lammps_get_thermo.restype = c_double
        return self.lib.lammps_get_thermo(self.lmp, name)

    # return total number of atoms in system

    def get_natoms(self):
        return self.lib.lammps_get_natoms(self.lmp)

    # set variable value
    # value is converted to string
    # returns 0 for success, -1 if failed

    def set_variable(self, name, value):
        if name: name = name.encode()
        if value: value = str(value).encode()
        return self.lib.lammps_set_variable(self.lmp, name, value)

    # reset simulation box size

    def reset_box(self, boxlo, boxhi, xy, yz, xz):
        cboxlo = (3 * c_double)(*boxlo)
        cboxhi = (3 * c_double)(*boxhi)
        self.lib.lammps_reset_box(self.lmp, cboxlo, cboxhi, xy, yz, xz)

    # return vector of atom properties gathered across procs
    # 3 variants to match src/library.cpp
    # name = atom property recognized by LAMMPS in atom->extract()
    # type = 0 for integer values, 1 for double values
    # count = number of per-atom valus, 1 for type or charge, 3 for x or f
    # returned data is a 1d vector - doc how it is ordered?
    # NOTE: need to insure are converting to/from correct Python type
    #   e.g. for Python list or NumPy or ctypes

    def gather_atoms(self, name, type, count):
        if name: name = name.encode()
        natoms = self.lib.lammps_get_natoms(self.lmp)
        if type == 0:
            data = ((count * natoms) * c_int)()
            self.lib.lammps_gather_atoms(self.lmp, name, type, count, data)
        elif type == 1:
            data = ((count * natoms) * c_double)()
            self.lib.lammps_gather_atoms(self.lmp, name, type, count, data)
        else:
            return None
        return data

    def gather_atoms_concat(self, name, type, count):
        if name: name = name.encode()
        natoms = self.lib.lammps_get_natoms(self.lmp)
        if type == 0:
            data = ((count * natoms) * c_int)()
            self.lib.lammps_gather_atoms_concat(self.lmp, name, type, count, data)
        elif type == 1:
            data = ((count * natoms) * c_double)()
            self.lib.lammps_gather_atoms_concat(self.lmp, name, type, count, data)
        else:
            return None
        return data

    def gather_atoms_subset(self, name, type, count, ndata, ids):
        if name: name = name.encode()
        if type == 0:
            data = ((count * ndata) * c_int)()
            self.lib.lammps_gather_atoms_subset(self.lmp, name, type, count, ndata, ids, data)
        elif type == 1:
            data = ((count * ndata) * c_double)()
            self.lib.lammps_gather_atoms_subset(self.lmp, name, type, count, ndata, ids, data)
        else:
            return None
        return data

    # scatter vector of atom properties across procs
    # 2 variants to match src/library.cpp
    # name = atom property recognized by LAMMPS in atom->extract()
    # type = 0 for integer values, 1 for double values
    # count = number of per-atom valus, 1 for type or charge, 3 for x or f
    # assume data is of correct type and length, as created by gather_atoms()
    # NOTE: need to insure are converting to/from correct Python type
    #   e.g. for Python list or NumPy or ctypes

    def scatter_atoms(self, name, type, count, data):
        if name: name = name.encode()
        self.lib.lammps_scatter_atoms(self.lmp, name, type, count, data)

    def scatter_atoms_subset(self, name, type, count, ndata, ids, data):
        if name: name = name.encode()
        self.lib.lammps_scatter_atoms_subset(self.lmp, name, type, count, ndata, ids, data)

    # create N atoms on all procs
    # N = global number of atoms
    # id = ID of each atom (optional, can be None)
    # type = type of each atom (1 to Ntypes) (required)
    # x = coords of each atom as (N,3) array (required)
    # v = velocity of each atom as (N,3) array (optional, can be None)
    # NOTE: how could we insure are passing correct type to LAMMPS
    #   e.g. for Python list or NumPy, etc
    #   ditto for gather_atoms() above

    def create_atoms(self, n, id, type, x, v, image=None, shrinkexceed=False):
        if id:
            id_lmp = (c_int * n)()
            id_lmp[:] = id
        else:
            id_lmp = id

        if image:
            image_lmp = (c_int * n)()
            image_lmp[:] = image
        else:
            image_lmp = image

        type_lmp = (c_int * n)()
        type_lmp[:] = type
        self.lib.lammps_create_atoms(self.lmp, n, id_lmp, type_lmp, x, v, image_lmp,
                                     shrinkexceed)

    @property
    def has_exceptions(self):
        """ Return whether the LAMMPS shared library was compiled with C++ exceptions handling enabled """
        return self.lib.lammps_config_has_exceptions() != 0

    @property
    def has_gzip_support(self):
        return self.lib.lammps_config_has_gzip_support() != 0

    @property
    def has_png_support(self):
        return self.lib.lammps_config_has_png_support() != 0

    @property
    def has_jpeg_support(self):
        return self.lib.lammps_config_has_jpeg_support() != 0

    @property
    def has_ffmpeg_support(self):
        return self.lib.lammps_config_has_ffmpeg_support() != 0

    @property
    def installed_packages(self):  #
        if self._installed_packages is None:
            self._installed_packages = []
            npackages = self.lib.lammps_config_package_count()
            sb = create_string_buffer(100)
            for idx in range(npackages):
                self.lib.lammps_config_package_name(idx, sb, 100)
                self._installed_packages.append(sb.value.decode())
        return self._installed_packages
