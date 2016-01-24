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
lfp1=np.array(h.ecpv.to_python())
lfp2=np.array(h.ecpv2.to_python())

vsum1=np.array(h.vsum1.to_python())
vsum2=np.array(h.vsum2.to_python())

"""
#anasig.signal = scipy.signal.filtfilt(B, A, anasig.signal)
#Chebyshev, Infinite Impulse Response.

http://www.tjhsst.edu/~sherbst/electronics/octave/dsp.html
"""
# Number of samplepoints
N= n = 40001*int(h.tstop)/1000 
#dt sample spacing
dt = 0.025 
f_s=1000/dt
#time vector
t = np.linspace(0.0, dt*N, N)
dtds=1/0.005
#tds = np.linspace(0.0, dtds*(N-1), (N-1)/200)

#Frequency vector
xf = np.linspace(0.0, dt*N, N/2)
f = np.linspace(0.0, dt*N, N)



plt.figure(1)

#interactive(True)
#plt.plot((xf,yf1))
plt.title("Raster Plot")
plt.hold(True)
colors=array([[0.42,0.67,0.84],[0.50,0.80,1.00],[0.90,0.32,0.00],[0.34,0.67,0.67],[0.42,0.82,0.83],[0.90,0.59,0.00],[0.33,0.67,0.47],[0.42,0.83,0.59],[0.90,0.76,0.00],[1.00,0.85,0.00],[0.71,0.82,0.41],[0.57,0.67,0.33],[1.00,0.38,0.60]]) # Colors for each cell population
j=len(colors)-1
#plt.plot(tvec,intervec,'bs')#plt.plot(tvec,pyra,'g^')
plt.plot(tvec,intervec,'.',c=colors[j], markeredgecolor = 'none')
j-=1
#,label='inhibitory interneurons',
#label='principle cells (pyramid)',
#Need a legend to describe cell type.
plt.plot(tvec,pyra,'.',c=colors[j], markeredgecolor = 'none')

j-=1
plt.plot(tvec,zerovec,'.',c=colors[j], markeredgecolor = 'none')


plt.hold(False)
maxtime=int(h.tstop)
#manrnpxtime=h.max(h.tvec)
plt.xlim(0,maxtime) # To convert to seconds
plt.ylim(0,numcell+10) # Just larger than the number of cells in the model
plt.ylabel("Neuron number")
plt.legend()
#plt.axis.set_xticks(r_[0:maxtime+1]) # Don't show half seconds
plt.show()
# Plot rasters

