import glob
import numpy as np
#from matplotlib import pyplot as plt

#calculate thermal conductance
def CalTC(deltaT,temp,dlist=0):
    print("Calculate thermal conductance.")
    deltaT=deltaT
    temperture=temp
    dlist=list(range(dlist))
    times=int(len(glob.glob('./kappa.*.bath*.run*.dat'))/2)
    kb=np.empty([2,times])

    for i in range(2):
        for j in range(times):
            kappafile="./kappa.*.bath"+str(i)+".run"+str(j)+"*.dat"
            for files in glob.glob(kappafile): 
                with open(files, 'r') as f: 
                    for line in f: 
                        kb[i][j]=line.split()[2]
#                        temperture=float(line.split()[1])
    heatflux=(kb[0]-kb[1])/2
    kappa=(kb[0]-kb[1])/2/(deltaT)
    kappa=np.delete(kappa,dlist)
    
    with open('thermalconductance.'+str(int(temperture))+'.dat','w') as f:
        f.write("Temperture\t"+str(temperture)+"\n"+"HeatFlux\t"+str(heatflux)+"\n"+"ThermalConductance\t"+str(kappa)+"\n"+"Mean\t"+str(np.mean(kappa))+"\n"+"StandardDeviation\t"+str(np.std(kappa))+"\n")

#powerspectra = np.loadtxt('power2.300.run1.dat')
#
#plt.plot(powerspectra[:,0],powerspectra[:,1])
#plt.show()
