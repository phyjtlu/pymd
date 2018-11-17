# using LAMMPS as a driver
from lammps import *
import numpy as N

lmp = lammps(infile='in.atoms')
lmp.start()

deltaq0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ]
deltaq1 = [0, 0, 0, 0, 0, 0, 1.4, -1.4, 0, 1.4, -1.4, 0, -1.4, 1.4, 0,]

# f0=lmp.f0
f0 = lmp.absforce(deltaq0)
deltaf0 = lmp.force(deltaq0)
f1 = lmp.absforce(deltaq1)
deltaf1 = lmp.force(deltaq1)

print lmp.type, '\n', lmp.els, '\n', f0, '\n', deltaf0, '\n', f1, '\n', deltaf1
lmp.quit()
