#!/applications/mbrsoft/bin/python

import sys
import numpy as N
import numpy.matlib

def chkShape(a):
    """
    check if a is a n by n matrix, if yes return n
    """
    aa = N.array(a)
    ash = N.shape(aa)
    if(ash[0] == ash[1]):
        return  ash[0]
    else:
        print "The matrix should be a n by n matrix"
        sys.exit(0)

def symmetrize(a):
    aa = N.array(a)
    return 0.5*(aa+N.transpose(aa))
    
def antisymmetrize(a):
    aa = N.array(a)
    return 0.5*(aa-N.transpose(aa))

def dagger(a):
    aa = N.array(a)
    ash = N.shape(aa)
    if(ash[0] != ash[1]):
        print "Not sqaure matrix"
        sys.exit(0)
    return N.transpose(N.conjugate(aa))

def hermitianize(a):
    aa = N.array(a)
    return 0.5*(aa+dagger(aa))

def mdot(* args):
    # dot product with arbitrary number of arguments
    tmp=N.identity(len(args[0]))
    for ii in range(len(args)):
        tmp=N.dot(tmp,args[ii])
    return tmp
    
def mm(* args):
    # mm with arbitrary number of arguments
    tmp=args[0].copy()
    for mat in args[1:]:
        tmp=N.dot(tmp,mat)
    return tmp

if __name__ == "__main__":
    import numpy as N
    test = [[1,-0.1],[-0.09,1]]
    print symmetrize(test)
    print antisymmetrize(test)
    test = [[1,-0.1],[-0.09,1,2]]
    print symmetrize(test)
    print antisymmetrize(test)

