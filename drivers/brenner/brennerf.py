import numpy as N
import brenner as PS

#1D array gives the coordinates of each atom
#[xi,yi,zi],i=1,na
x = [0.0, 0.0, 0.0,0.7, 0.7, 0.0,-0.7, 0.7, 0.0]
xa = N.array(x,dtype='d')
#1D array gives the unit cell
#x1,y1,z1,x2,y2,z2,x3,y3,z3
c = [10.0,10.0,10.0]
cell = N.array(c,dtype='d')

tt=N.array([6,6,6],dtype='i')
print tt

#set units to : Angstrom, eV
a,b,c=PS.brennerf(cell,tt,xa)
print a
print b
print c

