import glob
import numpy as np
#from matplotlib import pyplot as plt


def CalHF(dlist=1):
    # calculate average heat flux
    print("Calculate heat flux.")
    # temperture=temp
    for filename in glob.glob('./kappa.*.bath0.run0.dat'):
        with open(filename, 'r') as f:
            for line in f:
                temperture = float(line.split()[1])

    dlist = list(range(dlist))
    times = int(len(glob.glob('./kappa.*.bath*.run*.dat'))/2)
    kb = np.empty([2, times])

    for i in range(2):
        for j in range(times):
            kappafile = "./kappa." + \
                str(int(temperture))+".bath"+str(i)+".run"+str(j)+".dat"
            for files in glob.glob(kappafile):
                with open(files, 'r') as f:
                    for line in f:
                        kb[i][j] = line.split()[2]
#                        temperture=float(line.split()[1])
    oldkb = np.delete(kb, dlist, axis=1)
    balancekb = np.delete(kb, dlist, axis=1)
    for i in range(balancekb.shape[0]):
        for j in range(balancekb.shape[1]):
            balancekb[i][j] = np.mean(oldkb[i][0:j+1])

    heatflux = (balancekb[0]-balancekb[1])/2

    with open('heatflux.'+str(int(temperture))+'.dat', 'w') as f:
        f.write("Temperture\t"+str(temperture)+"\n"+"Bath0\t" +
                str(balancekb[0])+"\n"+"Bath1\t"+str(balancekb[1])+"\n"+"HeatFlux\t"+str(heatflux)+"\n"+"\n")


def CalTC(delta, dlist=0):
    # calculate thermal conductance
    print("Calculate thermal conductance.")
    delta = delta
    # temperture=temp
    for filename in glob.glob('./kappa.*.bath0.run0.dat'):
        with open(filename, 'r') as f:
            for line in f:
                temperture = float(line.split()[1])
    dlist = list(range(dlist))
    times = int(len(glob.glob('./kappa.*.bath*.run*.dat'))/2)
    kb = np.empty([2, times])

    for i in range(2):
        for j in range(times):
            kappafile = "./kappa." + \
                str(int(temperture))+".bath"+str(i)+".run"+str(j)+".dat"
            for files in glob.glob(kappafile):
                with open(files, 'r') as f:
                    for line in f:
                        kb[i][j] = line.split()[2]
#                        temperture=float(line.split()[1])
    kappa = (kb[0]-kb[1])/2/(delta*temperture)
    kappa = np.delete(kappa, dlist)
    # for i in range(len(kappa)):
    #    kappa[i]=np.mean(kappa[0:i+1])

    with open('thermalconductance.'+str(int(temperture))+'.dat', 'w') as f:
        f.write("Temperture\t"+str(temperture)+"\n"+"ThermalConductance\t"+str(kappa)+"\n" +
                "Mean\t"+str(np.mean(kappa))+"\n"+"StandardDeviation\t"+str(np.std(kappa))+"\n")


def CalPS():
    # calculate average powerspectra
    pass

#powerspectra = np.loadtxt('power2.300.run1.dat')
#
# plt.plot(powerspectra[:,0],powerspectra[:,1])
# plt.show()
