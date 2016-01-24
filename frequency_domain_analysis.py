


from numpy import *
import scipy
import numpy as np
import pylab as p1
from pylab import *
from matplotlib.pyplot import step, show, xlim
import re
#import neuron
#p = re.
from mpi4py import MPI

import time; tic=time.clock()
from scipy import loadtxt, size, shape, zeros, mod, floor, mean
from pylab import figure, plot, xlabel, ylabel, legend, xlim, ylim, show, hold, squeeze, sqrt
from bsmart import granger # Load the Granger calculation tool


datafft = abs(scipy.fft(lfp1))
datafft2 = datafft[1:100]
#p1.interactive(True)
ioff()  
title('now how much would you pay?')
xticklabels(fontsize=20, color='green')
draw(datafft2)
savefig('FFT', dpi=300)
close()


"""

ioff()      # turn updates off
title('now how much would you pay?')
xticklabels(fontsize=20, color='green')
draw()      # force a draw
savefig('FFT', dpi=300)
close()


print 'plottedFFT'

datafft = abs(scipy.fft(h.record_soma))
datafft2 = datafft[1:100]
p1.interactive(True)
p1.figure(figsize=(8, 6), dpi=80)
figure(1)
p1.plot(datafft2)
p1.show()

#ion()      # turn updating back on
#plot(rand(20), mfc='g', mec='r', ms=40, mew=4, ls='--', lw=3)
"""
