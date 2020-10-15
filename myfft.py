#!/usr/bin/env python
# -*- coding: utf-8 -*
import sys
import numpy as N

class myfft:
    '''
    '''

    def __init__(self, dt, n):
        self.dt = dt
        self.N = n
        self.dw = 2*N.pi/dt/n

    def Fourier1D(self, a):
        """ 1D Fourier transform from t to w using our definition:
            f(w) = \int f(t) e^Iwt dt
            In discrete form, it's
            f(j) = dt \sum_{i=0}^{N-1} f(i) e^(I 2\pi i j/N).
            N is length of the data. 
            d\omega/2/\pi = 1/N/dt

            This corresponds to the inverse FFT in numpy:
                numpy.fft.ifft(a)*2*\pi/d\omega
            """
        if(len(a) != self.N):
            print("MyFFT.Fourier1D: array length error!")
            sys.exit(0)
        else:
            nor = 2.*N.pi/self.dw
            b = N.fft.ifft(a)
            return nor*b

    def iFourier1D(self, a):
        """ 1D Fourier transform from w to t using our definition:
            f(t) = \int f(w) e^-Iwt dw/2/\pi
            In discrete form, it's
            f(i) = dw/2/\pi \sum_{j=0}^{N-1} f(j) e^(-I 2\pi i j/N).
            N is length of the data. 
            d\omega/2/\pi = 1/N/dt

            This corresponds to the FFT in numpy:
                numpy.fft.fft(a)*d\omega/2/\pi
            """
        if(len(a) != self.N):
            print("MyFFT.iFourier1D: array length error!")
            sys.exit(0)
        else:
            nor = self.dw/2/N.pi
            b = N.fft.fft(a)
            return nor*b


if __name__ == "__main__":
    import numpy as N

    a = N.zeros(1024)+1.
    a2 = N.zeros(1024)+2.
    al = N.array([a, a2])
    mfft = myfft(0.1, len(a))
    bl = list(map(mfft.Fourier1D, al))
    bll = N.transpose(list(map(mfft.iFourier1D, bl)))
    print(bll.shape)
    print((bll.T-al).max())
    a = N.arange(10)
    mfft = myfft(0.1, len(a))
    print(N.abs(N.array(mfft.iFourier1D(mfft.Fourier1D(a))-a)).max())
