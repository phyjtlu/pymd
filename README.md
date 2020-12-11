# Semi-classical Langevin Molecular Dynamics

This is a set of scripts to run semi-classical Langevin molecular dynamics on junctions connecting to several electron or phonon baths. The details of the method are presented in Prog. Surf. Sci. [<https://doi.org/10.1016/j.progsurf.2018.07.002]>.

To do the molecular dynamics, we need a force driver. Currently, we implemented Siesta, Brenner and LAMMPS potentials.

The important feature of this script is that, the statistics of the bath degrees of freeom is quantum. For example, the zero point fluctuations are included; they fulfill the quantum-mechanical Bose-Einstein distribution. Moreover, the electron bath is allowed to be in a nonequilibrium steady state (non-thermal).

## Install

Install LAMMPS and required python packages, which can be done by following script.

```bash
#!/bin/sh
# install.sh Compile & install LAMMPS shared library with Intel parallel studio & required packages in python
# Usage: sh install.sh
# -------------------------------------------------------
if [ -d 'lammps' ]
then
    echo 'LAMMPS already exists' 
else
    echo 'Download LAMMPS source code'
    git clone https://github.com/lammps/lammps.git -b stable
fi
echo 'Activate Intel compilation environment'
source /opt/intel/parallel_studio_xe_2020/psxevars.sh
cd lammps/src
echo 'Install select packages'
make yes-body yes-class2 yes-manybody yes-molecule yes-kspace yes-user-reaxc yes-user-phonon
echo 'Compile LAMMPS shared library'
make -j 8 mode=shlib intel_cpu
echo 'Install LAMMPS shared library in python'
make install-python
echo 'Compile the LAMMPS executable file'
make -j 8 intel_cpu
echo 'Install required python packages'
python -m pip install -U numpy netCDF4 tqdm
```

## Input files

examples/runmd.py: calculate thermal conductance of single molecular junction with quantum thermal bath molecular dynmics.
examples/runnegf.py: calculate thermal conductance of single molecular junction with none equilibrium Green' s function.

## Output files
