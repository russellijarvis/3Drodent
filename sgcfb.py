#os.chdir('/home/zaza3/scipy-master/scipy/fftpack')
#execfile('setup.py')
#os.chdir(workingdir)

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

outi=downsample(time_courses[int(indexi)],oldrate=40000,newrate=200)
outj=downsample(time_courses[int(indexj)],oldrate=40000,newrate=200)

F,pp,cohe,Fx2y,Fy2x,Fxy=granger(outi,outj,20)#,200,100)# nfs2[0], nfs2[1] ,
###
plt.hold(False) #Matlab style hold
plt.figure(figure_cnt)
###

# Plot Granger spectra
plt.hold(True) #Matlab style hold
labels=list()
labels=['Fx2y','Fy2x','Fxy','cohe','pp']
alldata=list()
colors=[[0,0,1],[1,0,0],[0,0.5,0.5],[0.5,1,0],[1,0.5,1]]

plt.plot(F,Fxy,label='causality of channel X to channel Y',linewidth=4,c=colors[0])
plt.plot(F,Fy2x,label='causality of channel Y to channel X',linewidth=4,c=colors[1])
#plt.xlim(0,2)
#plt.ylim(-0.4,0.4)
plt.xlabel("Hz")
plt.xlim(0,100)
plt.ylim(0,1.5)
plt.title('GSC components'+bc)
plt.legend()
plt.hold(False) #Matlab style hold
figure_cnt=figure_cnt+1
###

plt.figure(figure_cnt)
###
plt.hold(True) #Matlab style hold

plt.plot(F,pp,label='spectral power of channel X to channel Y',linewidth=4,c=colors[2])
plt.plot(F,cohe,label='coherance between two signals Y to channel X',linewidth=4,c=colors[3])


#plt.xlim(0,40)
#plt.ylim(0,0.5)
plt.xlim(0,100)
plt.ylim(0,1.5)
plt.xlabel("Hz")
plt.title('GSC components'+bc)
plt.legend()
plt.hold(False) #Matlab style hold

figure_cnt=figure_cnt+1
###
plt.figure(figure_cnt)
F,pp,cohe,Fx2y,Fy2x,Fxy=granger(time_courses[int(indexi)],time_courses[int(indexj)],20)#,200,100)# nfs2[0], nfs2[1] ,

plt.plot(F,pp,linewidth=4,c=colors[4])
plt.title('Spectral power'+bc)
plt.xlim(0,60)
plt.ylim(0,30)
plt.xlabel("Hz")
