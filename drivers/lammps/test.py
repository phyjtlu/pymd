#using LAMMPS as a driver
from lammps import lammps
import numpy as N

lmp = lammps(label="test", infile="in.test")
lmp.start()
#print lmp.f0
#print lmp.newxyz[:]

qt = N.zeros(9)
qt[0]=0.0

#print "new displacment:", qt

#print "new force:" 
print lmp.force(qt)

lmp.quit()
