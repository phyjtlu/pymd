#!/applications/mbrsoft/bin/python

import numpy as N
import pysiesta as PS

from units import *

#print PS.__doc__
#print PS.pysiestaforce.__doc__
#print PS.pysiestalaunch.__doc__
#print PS.pysiestaunits.__doc__

# Map atomnumbers into elemental labels
#PeriodicTable = {'H':1,1:'H','D':1001,1001:'D','He':2,2:'He','Li':3,3:'Li','Be':4,4:'Be','B':5,5:'B','C':6,6:'C','N':7,7:'N','O':8,8:'O','F':9,9:'F','Ne':10,10:'Ne','Na':11,11:'Na','Mg':12,12:'Mg','Al':13,13:'Al','Si':14,14:'Si','P':15,15:'P','S':16,16:'S','Cl':17,17:'Cl','Ar':18,18:'Ar','K':19,19:'K','Ca':20,20:'Ca','Sc':21,21:'Sc','Ti':22,22:'Ti','V':23,23:'V','Cr':24,24:'Cr','Mn':25,25:'Mn','Fe':26,26:'Fe','Co':27,27:'Co','Ni':28,28:'Ni','Cu':29,29:'Cu','Zn':30,30:'Zn','Ga':31,31:'Ga','Ge':32,32:'Ge','As':33,33:'As','Se':34,34:'Se','Br':35,35:'Br','Kr':36,36:'Kr','Rb':37,37:'Rb','Sr':38,38:'Sr','Y':39,39:'Y','Zr':40,40:'Zr','Nb':41,41:'Nb','Mo':42,42:'Mo','Tc':43,43:'Tc','Ru':44,44:'Ru','Rh':45,45:'Rh','Pd':46,46:'Pd','Ag':47,47:'Ag','Cd':48,48:'Cd','In':49,49:'In','Sn':50,50:'Sn','Sb':51,51:'Sb','Te':52,52:'Te','I':53,53:'I','Xe':54,54:'Xe','Cs':55,55:'Cs','Ba':56,56:'Ba','La':57,57:'La','Ce':58,58:'Ce','Pr':59,59:'Pr','Nd':60,60:'Nd','Pm':61,61:'Pm','Sm':62,62:'Sm','Eu':63,63:'Eu','Gd':64,64:'Gd','Tb':65,65:'Tb','Dy':66,66:'Dy','Ho':67,67:'Ho','Er':68,68:'Er','Tm':69,69:'Tm','Yb':70,70:'Yb','Lu':71,71:'Lu','Hf':72,72:'Hf','Ta':73,73:'Ta','W':74,74:'W','Re':75,75:'Re','Os':76,76:'Os','Ir':77,77:'Ir','Pt':78,78:'Pt','Au':79,79:'Au','Hg':80,80:'Hg','Tl':81,81:'Tl','Pb':82,82:'Pb','Bi':83,83:'Bi','Po':84,84:'Po','At':85,85:'At','Rn':86,86:'Rn','Fr':87,87:'Fr','Ra':88,88:'Ra','Ac':89,89:'Ac','Th':90,90:'Th','Pa':91,91:'Pa','U':92,92:'U','Np':93,93:'Np','Pu':94,94:'Pu','Am':95,95:'Am','Cm':96,96:'Cm','Bk':97,97:'Bk','Cf':98,98:'Cf','Es':99,99:'Es','Fm':100,100:'Fm','Md':101,101:'Md','No':102,102:'No'}

#AtomicMassTable={'H':1.00794, 'He':4.002602, 'Li':6.941, 'Be':9.012182, \
#    'B':10.811, 'C':12.0107, 'N':14.0067, 'O':15.9994, \
#    'F':18.9984032, 'Ne':20.1791, 'Na':22.98976928, 'Mg':24.3050, \
#    'Al':26.9815386, 'Si':28.0855, 'P':30.973762, 'S':32.065, \
#    'Cl':35.453, 'Ar':39.948, 'K':39.0983, 'Ca':40.078, \
#    'Sc':44.955912, 'Ti':47.867, 'V':50.9415, 'Cr':51.9961, \
#    'Mn':54.938045, 'Fe':55.845, 'Co':58.933195, 'Ni':58.6934, \
#    'Cu':63.546, 'Zn':65.38, 'Ga':69.723, 'Ge':72.64, \
#    'As':74.92160, 'Se':78.96, 'Br':79.904, 'Kr':83.798, \
#    'Rb':85.4678, 'Sr':87.62, 'Y':88.90585, 'Zr':91.224, \
#    'Nb':92.90638, 'Mo':95.96, 'Tc':98, 'Ru':101.07, \
#    'Rh':102.90550, 'Pd':106.42, 'Ag':107.8682, 'Cd':112.411, \
#    'In':114.818, 'Sn':118.710, 'Sb':121.760, 'Te':127.60, \
#    'I':126.90447, 'Xe':131.293, 'Cs':132.9054519, 'Ba':137.327, \
#    'La':138.90547, 'Ce':140.116, 'Pr':140.90765, 'Nd':144.242, \
#    'Pm':145, 'Sm':150.36, 'Eu':151.964, 'Gd':157.25, \
#    'Tb':158.92535, 'Dy':162.500, 'Ho':164.93032, 'Er':167.259, \
#    'Tm':168.93421, 'Yb':173.054, 'Lu':174.9668, 'Hf':178.49, \
#    'Ta':180.94788, 'W':183.84, 'Re':186.207, 'Os':190.23, \
#    'Ir':192.217, 'Pt':195.084, 'Au':196.966569, 'Hg':200.59, \
#    'Tl':204.3833, 'Pb':207.2, 'Bi':208.98040, 'Po':209, \
#    'At':210, 'Rn':222, 'Fr':223, 'Ra':226, 'Ac':227, \
#    'Th':232.03806, 'Pa':231.03586, 'U':238.02891, 'Np':237, \
#    'Pu':244, 'Am':243, 'Cm':247, 'Bk':247, 'Cf':251, \
#    'Es':252, 'Fm':257, 'Md':258, 'No':259, 'Lr':262, \
#    'Rf':265, 'Db':268, 'Sg':271, 'Bh':272, 'Hs':270, \
#    'Mt':276, 'Ds':281, 'Rg':280, 'Cn':285, 'Uut':284, \
#    'Uuq':289, 'Uup':288, 'Uuh':293, 'Uus':294, 'Uuo':294}

class siesta:
    """
    construct necessory information to do a siesta force run

    assuming all the fdf files are available
    """
    def __init__(self,label,xyz,cell,mesh=100.,dmtol=0.001,\
                 constraints=[],tdir="./",lunit="Ang",eunit="eV"):
        self.md2ang = 0.06466;
        self.meshcutoff = mesh
        self.dmtol = dmtol
        self.constraints = constraints
        self.label = label
        #1D array gives the coordinates of each atom
        #[xi,yi,zi],i=1,na
        self.xyz = N.array([a[1:] for a in xyz],dtype='d').flatten()
        self.els = [a[0] for a in xyz]
        self.conv = self.md2ang*N.array([3*[1.0/N.sqrt(AtomicMassTable[a])]\
                                         for a in self.els]).flatten()
        #1D array gives the unit cell
        #x1,y1,z1,x2,y2,z2,x3,y3,z3
        self.cell = N.array(cell,dtype='d').flatten()
        self.lunit=lunit
        self.eunit=eunit
        self.genfdf(tdir)


    def genfdf(self,tdir="./"):
        """
        generate the fdf files.
        It includes extra fdf files:
            STRUCT.fdf
            Default.fdf
        """
        fname = self.label+".fdf"
        fn = open(fname,"w")

        fn.write("#fdf generated by siesta:genfdf\n")
        fn.write("SystemName   "+self.label+"\n")
        fn.write("SystemLabel   "+self.label+"\n")
        fn.write("MD.TypeOfRUN   forces\n")
        fn.write("MeshCutoff    "+str(self.meshcutoff)+" Ry\n")
        fn.write("DM.Tolerance  "+str(self.dmtol)+"\n\n\n")
        for i in range(len(self.constraints)):
            if i == 0:
                fn.write("%block GeometryConstraints\n")
            fn.write("position from "+str(self.constraints[i][0])+" to\
                     "+str(self.constraints[i][1])+"\n")
            if i == len(self.constraints)-1:
                fn.write("%endblock GeometryConstraints\n")
        fn.write("%include STRUCT.fdf\n")
        #fn.write("%include "+tdir+"MD.fdf\n")
        fn.write("%include "+tdir+"Default.fdf\n")
        fn.close()
                     





    def start(self,np=1):
        """
        start siesta
        np - number of cores
        """
        #set units to : Angstrom, eV
        PS.pysiestaunits("Ang", "eV")

        #PS.pysiestalaunch(label, 8, 'mpirun -np')	#parallel siesta
        PS.pysiestalaunch(self.label,np,'mpirun -np')	#serial siesta
        #PS.pysiestalaunch(self.label)	#serial siesta
        print "siesta launched!"

        print "running test..."
        #energy,force = PS.pysiestaforce( self.label , self.xyz, self.cell )
        #equilibrium force
        self.initforce()

        #print "siesta equilibrium energy:", energy
        print "siesta equilibrium force:", self.f0 
        print "test finished!"

    def quit(self):
        """
        quit siesta
        """
        print "Quit siesta!"
        PS.pysiestaquit(self.label)
        print "Done!"

    def newx(self,q):
        """
        return the real coordinates from displacements got from MD

        performing unit conversion and remove the mass factor
        """
        return self.xyz + self.conv*q

    def absforce(self,q):
        """
        calculate the force from siesta
        q:  displacement list of all atoms, including those fixed
        """
        energy,force = PS.pysiestaforce(self.label , self.newx(q), self.cell)
        return self.conv*force

    def initforce(self):
        """
        """
        print "Calculate zero displacement force..."
        #equilibrium force
        extq = N.zeros(len(self.xyz))
        self.f0=self.absforce(extq)

    def force(self,q):
        """
        calculate the relative force 
        q:  displacement list of all atoms, including those fixed
        """
        return self.absforce(q)-self.f0
        
        

if __name__=="__main__":
    import Scientific.IO.NetCDF as nc

    def ReadMDNCFile(filename):
        """
        Reads a NetCDF file 
        """
        class mdmath:
            pass
    
        file = nc.NetCDFFile(filename,'r')
        print 'Reading from %s' % filename
    
        # General attributes
        mdmath.filename = filename
        mdmath.cell= N.array(file.variables['UnitCell'])
        mdmath.xyz = N.array(file.variables['XYZ'])
        mdmath.dynatom = N.array(file.variables['DynamicAtoms'])
        mdmath.atomlist= N.array(file.variables['AtomList'])
    
        file.close()
    
        return mdmath 

    mdmath = ReadMDNCFile("./MD1.nc")


    xyz = [["Au",a[0],a[1],a[2]] for a in mdmath.xyz]          

    tests = siesta("test",xyz,mdmath.cell)
    print tests.els

    #tests.start()
    #tests.quit()

