# using LAMMPS as a driver
from lammps import *
import numpy as N
import re

infile = 'in.test'
datafile = 'test.dat'

lines = open(datafile, 'r').readlines()
data = []
for line in lines:
    data.append(re.split(r'[;,\s]\s*', line))

id1 = float(data[20][1])
#mass1 = float(data[20][2])
mass1 = 6
#id2 = float(data[21][1])
#mass2 = float(data[21])

anr = []
xyzonly = []
for i in range(24, len(lines)):
    if float(data[i][2]) == id1:
        anr.append(mass1)
    #if float(data[i][2])==id2:anr.append(mass2)
    xyzonly.append([float(data[i][5]), float(data[i][6]), float(data[i][7])])

xyz = [[PeriodicTable[a], b[0], b[1], b[2]] for (a, b) in zip(anr, xyzonly)]

lmp = lammps(label="test", infile=infile, xyz=xyz)
lmp.start()

#mass0 = lmp.extract_atom("mass", 2)
#mass1 = N.array(lmp.gather_atoms("mass", 1, 1))

#print lmp.f0,'\n'
# print lmp.extract_atom("mass", 2)[0], lmp.extract_atom("mass", 2)[
#    1], lmp.extract_atom("mass", 2)[2], lmp.extract_atom("mass", 2)[3], '\n'
# print lmp.gather_atoms("mass", 1, 1)[0], lmp.gather_atoms("mass", 1, 1)[
#    1], lmp.gather_atoms("mass", 1, 1)[2], lmp.gather_atoms("mass", 1, 1)[3], '\n'

deltaq0 = [0, 0, 0, 0, 0, 0, 0, 0, 0]
deltaq1 = [0, 0, 0, 0, 0, 0, 1.4, -1.4, 0]

# f0=lmp.f0
f0 = lmp.absforce(deltaq0)
deltaf0 = lmp.force(deltaq0)
f1 = lmp.absforce(deltaq1)
deltaf1 = lmp.force(deltaq1)

print f0, '\n', deltaf0, '\n', f1, '\n', deltaf1
lmp.quit()
