#!/applications/mbrsoft/bin/python

import sys
import numpy as N
import units as U

N.seterr(over="ignore")

def coth(x):
    """
    hyperbolic cotangent function
    """
    if x == 0.0:
        print "coth:coth(0) is infinity"
        sys.exit(0)
    else:
        return N.cosh(x)/N.sinh(x)

def xcoth(x):
    """
    x*coth(x)
    """
    if x == 0.0:
        return 1.0
    else:
        return x*N.cosh(x)/N.sinh(x)

def bose(w,T):
    """
    bose distribution
    """
    small = 10e-20
    if T == 0.0:
        if w == 0.0:
            return 1/(N.exp(1.0/U.kb)-1)
        elif w < 0.0:
            return -1.0
        else:
            return 0.0
    else:
        if w == 0.0:
            #return 1/small
            #have problems for finite temperature for bias calculation
            #return 0 seems solves it
            return 0.0
        else:
            return 1.0/(N.exp(w/U.kb/T)-1.0)

def fermi(ep,mu,T):
    """
    fermi distribution
    """
    if T == 0.0:
        if(ep < mu):
            return 1.0
        elif(ep > mu):
            return 0.0
        else:
            return 0.5
    else:
        return 1/(N.exp((ep-mu)/U.kb/T)+1)

def flinterp(x,xs,ys):
    """
    do a linear interpolation of (xs,ys), return the interpolated value at x.
    """
    id=nearest(x,xs)

    #boundaries
    if id == len(xs)-1:
        return ys[-1]
    if id == 0:
        return ys[0]

    #linear interpolation
    dd=xs[id]-x
    if dd<0:
        return ys[id]+dd/(xs[id]-xs[id-1])*(ys[id]-ys[id-1])
    else:
        return ys[id]+dd/(xs[id]-xs[id+1])*(ys[id]-ys[id+1])
        
def nearest(b,bs):
    """
    return the index of the element in bs wich is the nearest to b.
    """
    bsn = N.array(bs)
    bst = abs(bsn-b)
    return(list(bst).index(min(bst)))

def rpadleft(bs,b):
    if len(bs) > 1:
        return N.concatenate((N.array([b]),N.array(bs)[:-1]),axis=0)
    elif len(bs) == 1:
        return N.array([b])
    else:
        print "len(bs) is less than 1"
        sys.exit()

def mm(* args):
    tmp=args[0].copy()
    for ii in range(1,len(args)):
        tmp=N.dot(tmp,args[ii])
    return tmp

