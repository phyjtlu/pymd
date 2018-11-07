import brenner 
import numpy as np

#---------------------------------------------------------------------------------#
#                          	    Test System					  #
#---------------------------------------------------------------------------------#
# [Ftest,EAtest, ETOTtest] = brennermex([3,3,100],[0,0,0,0,1,0,1,0,0], [0,0,0], [0,0,0,0,0,0,0,0,0]);
# call brennerf(NA,cell,an,R, eain,Fin,etotTest)
NA=np.array(3,dtype='i')				# number of atoms.
cell=np.array([3,3,100],dtype='d')		# Unitcell dimensions.
R=np.array([0,0,0,0,1,0,1,0,0],dtype='d')		# Position vector.
typevec=np.array([6,6,6],dtype='i')		# type of atoms (atomic numbers).
eain=np.array([0,0,0],dtype='d')			# Energies input.
Fin=np.array([0,0,0,0,0,0,0,0,0],dtype='d')	# Force vector.
etotTest=np.array([0,0,0],dtype='d')
#print 'was here',cell,cell.shape

#---------------------------------------------------------------------------------#
#                          	    Call Force					  #
#---------------------------------------------------------------------------------#
#npfromC,CellfromC,AnfromC,RfromC,Ea2C,F2C,Etot2C=py_brenner_module.brennerf(NA,cell,typevec,R,eain,Fin,etotTest)
#py_brenner_module.brennerf(NA,cell,typevec,R,eain,Fin,etotTest)
#NA,typevec,R,eain,Fin,etotTest=
#u1,u2,u3=py_brenner_module.brennerf(3,cell,typevec,R,eain,Fin,etotTest)
#py_brenner_module.brennerf(1,np.array([1,1,1]),[10,10,10],int(6),np.array([0,0,0]),[0,0,0],[0,0,0])
#u1,u2,u3=py_brenner_module.brennerf([10,10,10],int(6),[0,0,0],[0],[0,0,0],[0,0,0])
#u1,u2,u3=py_brenner_module.brennerf([10,10,10],typevec,R,[0,0,0])
#print u1,type(u1)




u1,u2,u3=brenner.brennerf(cell,typevec,R)
print "Energy of each atom (eV):", u1
print "Force vector (eV/Ang):", u2 
print "Total energy (eV):", u3
