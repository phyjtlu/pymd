import time
import numpy as np
from negf import bpt
from matplotlib import pyplot as plt
lammpsinfile = [
    'atom_style full',
    'units metal',
    'boundary f p p',
    'read_data structure.data',
    'pair_style rebo',
    'pair_coeff * * CH.rebo C H',
]
time_start = time.time()
atomfixed = [range(0*3, (19+1)*3), range(181*3, (200+1)*3)]
atomofbath = [range(20*3, (69+1)*3), range(131*3, (180+1)*3)]
mybpt = bpt(lammpsinfile, 0.25, 0.1, atomofbath, atomfixed, 100)
mybpt.plotresult()
# T_H/C = T*(1±delta/2)
T = [100, 200, 300, 400, 500, 600, 700,
     800, 900, 1000]
delta = 0.1
thermalconductance = []
for temp in T:
    thermalconductance.append([temp, mybpt.thermalconductance(temp, delta)])
np.savetxt('thermalconductance.dat', thermalconductance)
plt.figure(5)
plt.plot(np.array(thermalconductance)[
    :, 0], np.array(thermalconductance)[:, 1])
plt.xlabel('Temperature(K)')
plt.ylabel('Thermal Conductance(nW/K)')
plt.savefig('thermalconductance.png')
time_end = time.time()
print('time cost', time_end-time_start, 's.')
