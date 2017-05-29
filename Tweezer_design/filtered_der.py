# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 22:33:49 2016

@author: Antoine
"""
from numpy import append as npappend
from numpy import diff,interp
from numpy import pi as np_pi
#from scipy.fftpack import diff as fftpackdiff
from scipy.signal import butter as signalbutter
from scipy.signal import filtfilt as signalfiltfilt
import matplotlib.pyplot as plt
from copy import deepcopy

twopi = 2*np_pi




def periodic_derivative(x,y,max_periods):
    plot = False
    Ns =len(x)
    b,a = signalbutter(8,2.0*max_periods/Ns)
    ymid =interp(x+0.5*(x[1]-x[0]),x,y,period=2*np_pi)
    yder = diff(ymid)/diff(x) 
    #yder = Ns/(max(x)-min(x))*fftpackdiff(y,1,Ns)
    yder_filt = deepcopy(yder)
    x_filt = deepcopy(x)
    x_filt = npappend(x_filt,x_filt[-1]+x_filt[1]-x_filt[0])
    yder_filt = signalfiltfilt(b,a,npappend(yder_filt,yder_filt[0]))
    if plot:
        plt.figure(1)
        plt.subplot(311)
        plt.plot(x, y)
        
        plt.subplot(312)
        plt.plot(x[0:-1],yder)
        
        plt.subplot(313)
        plt.plot(x_filt[0:-1],yder_filt)
        plt.show()
    return yder_filt
        
        
#x = numpy.array(range(100))
#y = numpy.sin(twopi*x/100)+numpy.sin(twopi*x/10)
#periodic_derivative(x,y,4)