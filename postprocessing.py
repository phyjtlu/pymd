import glob
import numpy as np
#from matplotlib import pyplot as plt

#calculate thermal conductance
def CTC(delta,temp,dlist=list(range(1))):
    print("Calculate thermal conductance.")
    delta=delta
    temperture=temp
    dlist=dlist
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

    kappa=(kb[0]-kb[1])/2/(2*delta*temperture)
    kappa=np.delete(kappa,dlist)
    
    with open('thermalconductance.'+str(int(temperture))+'.dat','w') as f:
        f.write("Temperture\t"+str(temperture)+"\n"+"Thermal Conductance\t"+str(kappa)+"\n"+"Mean\t"+str(np.mean(kappa))+"\n"+"Standard Deviation\t"+str(np.std(kappa))+"\n")

#powerspectra = np.loadtxt('power2.300.run1.dat')
#
#plt.plot(powerspectra[:,0],powerspectra[:,1])
#plt.show()
