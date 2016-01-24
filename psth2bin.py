from neuron import h
import pylab as p1
from pylab import *
from matplotlib.pyplot import step, show, xlim
#import re
size=h.tstop*1000
size=int(size)
binary_train=[0]*size
indexs=[]
print len(binary_train)


indexs=h.vecin.to_python() #.c.mul(1000)

indexs2=[0]*len(h.vecin)
print len(indexs2), "indexs2"
for i in xrange(0,(len(indexs2)-1)):
 print i 
 indexs2[i]=int(indexs[i]*1000)


for i in xrange(0,(len(indexs2)-1)): 
 binary_train[indexs2[i]]=1
 print binary_train[indexs2[i]]
 
# scatter plot from in here.
print indexs2
print len(binary_train)

plot(binary_train)
show()
