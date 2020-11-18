import numpy as np


def dumpavetraj(lammpsdata, trajectoriesfiles, outputname="avestructure.data"):
    from ovito.io import export_file, import_file

    # export average atom postions from MD trajectories files
    print("import LAMMPS data file %s" % lammpsdata)
    data = import_file(lammpsdata).source.compute(0)
    aveposition = np.zeros(
        [len(trajectoriesfiles), data.number_of_particles, 3])
    for i, trajfile in enumerate(trajectoriesfiles):
        print("import trajectorie file %s" % trajfile)
        traj = import_file(trajfile, columns=[
            "Particle Type", "Position.X", "Position.Y", "Position.Z", "Force.X", "Force.Y", "Force.Z"])
        for frame_index in range(traj.source.num_frames):
            position = np.array(traj.source.compute(
                frame_index).particles.positions)
            aveposition[i] = (aveposition[i]*frame_index +
                              position)/(frame_index+1)
    data.particles_.positions_[:] = np.mean(aveposition, axis=0)
    print("export LAMMPS data file %s" % outputname)
    export_file(data, outputname, "lammps/data", atom_style="full")


def calHF(dlist=1):
    import glob

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

    np.savetxt('heatflux.'+str(int(temperture))+'.dat',np.transpose((balancekb[0],balancekb[1],(balancekb[0]-balancekb[1])/2)),header="Bath0 Bath1 heatflux")



def calTC(delta, dlist=1):
    import glob

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

    np.savetxt('thermalconductance.'+str(int(temperture))+'.dat',(np.mean(kappa),np.std(kappa)),header="Mean Std")


if __name__ == "__main__":
    lammps = "structure.data"
    trajectories = ["trajectories.300.run0.ani", "trajectories.300.run1.ani",
                    "trajectories.300.run2.ani", "trajectories.300.run3.ani", "trajectories.300.run4.ani"]
    dumpavetraj(lammps, trajectories)
