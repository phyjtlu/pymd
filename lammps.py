#LAMMPS driver for Langevin molecular dynamics
#Adapted from LAMMPS python wraper
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

#atomic table and periodic table from units.py
from units import *

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

    def __init__(self, label, infile, mesh=100., dmtol=0.001, \
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
        self.meshcutoff = mesh
        self.dmtol = dmtol
        self.constraints = constraints
        self.label = label
        self.lunit = lunit
        self.eunit = eunit

    def start(self, np=1):
        print("lammps launched")
        lines = open(self.infile, 'r').readlines()
        for line in lines: self.command(line)
        self.newxyz = self.xyz = self.gather_atoms("x", 1, 3)
        self.number = self.get_natoms()
        #self.els = self.extract_atom("mass",2)[1]
        self.els = self.gather_atoms("mass",1,1)[1]
        self.conv = 1.
        self.initforce()

    def quit(self):
        self.close()
        print("Quit lammps!")

    def newx(self, q):
        for i in range(3*self.number):
            # self.newxyz[i] = self.xyz[i]+self.conv*q[i]
            self.newxyz[i] = self.xyz[i] + self.conv * q[i]
        return self.newxyz

    #def absforce(self, q):
    #    self.scatter_atoms("x", 1, 3, self.newx(q))
    #    self.command("run 1")
    #    #self.absf = N.zeros((self.get_natoms(), 3), dtype=N.float_)
    #    #self.lmpf = self.extract_atom("f", 3)
    #    #for m in range(self.get_natoms()):
    #    #    for n in range(3):
    #    #        self.absf[m][n] = self.lmpf[m][n]
    #    # return self.conv * N.array(self.gather_atoms("f", 1, 3)) wrong!
    #    self.absf = N.zeros(3*self.number, dtype=N.float_)
    #    self.lmpf = self.extract_atom("f", 3)
    #    for m in range(self.number):
    #        for n in range(3):
    #            self.absf[3*m+n] = self.lmpf[m][n]
    #    return self.conv * self.absf

    def absforce(self, q):
        self.scatter_atoms("x", 1, 3, self.newx(q))
        self.command("run 1")
        return self.conv * N.array(self.gather_atoms("f",1,3))
        
    def initforce(self):
        print("Calculate zero displacement force")
        extq = N.zeros(len(self.xyz))
        self.f0 = self.absforce(extq)

    def force(self, q):
        return self.absforce(q) - self.f0

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
