# using LAMMPS as a driver
from lammps import *
import numpy as N

infile = 'in.test'

lmp = lammps(label="test", infile=infile)
lmp.start()

deltaq0 = [0, 0, 0, 0, 0, 0, 0, 0, 0]
deltaq1 = [0, 0, 0, 0, 0, 0, 1.4, -1.4, 0]

# f0=lmp.f0
f0 = lmp.absforce(deltaq0)
deltaf0 = lmp.force(deltaq0)
f1 = lmp.absforce(deltaq1)
deltaf1 = lmp.force(deltaq1)

print f0, '\n', deltaf0, '\n', f1, '\n', deltaf1
lmp.quit()
