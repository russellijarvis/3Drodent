
import matplotlib
#matplotlib.use('Agg')# //Must be before importing matplotlib.pyplot or pylab!
#import Image
#import matplotlib.pyplot as plt

import pylab as plt
from pylab import *

import numpy as np
#import scipy does not work.
import scipy as scipy
import time; tic=time.clock()
#from scipy import loadtxt, size, shape, zeros, mod, floor, mean
from pylab import figure, plot, xlabel, ylabel, legend, xlim, ylim, show, hold, squeeze, sqrt
from bsmart import granger # Load the Granger calculation tool
import nrn
import hoc
h=hoc.HocObject()


#x=zip(time_courses[int(indexi)], time_courses[int(indexj)] )#x is the data for at least two channels, e.g. a 2x8000 array consisting of two LFP time series
nctrls=1	#ntrls is the number of trials (whatever that means -- just leave it at 1)
npts=40001	#npts is the number of points in the data (in this example, 8000)
p=20	#p is the order of the polynomial fit (e.g. 10 for a smooth fit, 20 for a less smooth fit)
fs=40000	#fs is the sampling rate (e.g. 200 Hz)
freq=fs/2	#freq is the maximum frequency to calculate (e.g. fs/2=100, which will return 0:100 Hz)

bc='between cells i= '+str(int(indexi)) +'between cells j= '+str(int(indexj))+ 'ff= '+str(int(h.ff))+'prune net= '+str(int(h.prunenet))

print "Plotting..."

#def granger(vec1,vec2,order=10,rate=200,maxfreq=0):

execfile('pyhoc.py')

plt.figure(figure_cnt)
###
plt.hold(True) #Matlab style hold

plt.plot(F,pp,label='spectral power of channel X to channel Y',linewidth=2,c=colors[2])
plt.plot(F,cohe,label='coherance between two signals Y to channel X',linewidth=2,c=colors[3])


#plt.xlim(0,40)
#plt.ylim(0,0.5)
plt.xlim(0,60)
plt.ylim(0,20)
plt.xlabel("Hz")
plt.title('GSC components'+bc)
plt.legend()
plt.hold(False) #Matlab style hold

figure_cnt=figure_cnt+1
###
plt.figure(figure_cnt)
F,pp,cohe,Fx2y,Fy2x,Fxy=granger(time_courses[int(indexi)],time_courses[int(indexj)],20)#,200,100)# nfs2[0], nfs2[1] ,
plt.title('Spectral power'+bc)
plt.plot(F,pp)
plt.xlim(0,60)
plt.ylim(0,30)
plt.xlabel("Hz")

