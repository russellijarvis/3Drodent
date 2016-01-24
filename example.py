
import scipy.io as sio
Matrixnp=np.array(Matrix)#Numpy matrix is very flexible. ie a=Matrix2[295,0:400]  
Weightsnp=np.array(Weights)#Numpy matrix is very flexible. ie a=Matrix2[295,0:400]  
print Matrixnp
print Weights
#sio.savemat
sio.savemat('adj.mat', {'Matrixnp':Matrixnp})
sio.savemat('weights.mat', {'Weights':Weights})
"""
import pylab as pl
import numpy as np
#import scipy 
#import loadtxt, size, shape, zeros, mod, floor, mean


#import time; tic=time.clock()
#from pylab #import figure, plot, xlabel, ylabel, legend, xlim, ylim, show, hold, squeeze, sqrt
#from bsmart #import granger # Load the Granger calculation tool


#datafft = abs(scipy.fft(lfp1))
#datafft2 = datafft[1:100]
y = np.arange(10)
pl.figure()
pl.plot(y)
#pl.show()
pl.savefig('y=' + str(10) + '.png')
pl.close()
"""

