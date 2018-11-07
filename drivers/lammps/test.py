#using LAMMPS as a driver
from lammps import lammps
import numpy as N

lmp = lammps(label="test", infile="in.test")
lmp.start()
print lmp.f0
lmp.quit()
