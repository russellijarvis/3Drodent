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
#import nrn
#import hoc
#h=hoc.HocObject()
from neuron import h

#x=zip(time_courses[int(indexi)], time	_courses[int(indexj)] )#x is the data for at least two channels, e.g. a 2x8000 array consisting of two LFP time series
nctrls=1	#ntrls is the number of trials (whatever that means -- just leave it at 1)
npts=40001	#npts is the number of points in the data (in this example, 8000)
p=20	#p is the order of the polynomial fit (e.g. 10 for a smooth fit, 20 for a less smooth fit)
fs=40000	#fs is the sampling rate (e.g. 200 Hz)
freq=fs/2	#freq is the maximum frequency to calculate (e.g. fs/2=100, which will return 0:100 Hz)

print "Plotting..."

#def granger(vec1,vec2,order=10,rate=200,maxfreq=0):

fin=0
execfile('pyhoc.py')

for indexi in xrange(0,numcell):
 for indexj in xrange(0,numcell):
   outi=downsample(time_courses[int(indexi)],oldrate=40000,newrate=200)
   outj=downsample(time_courses[int(indexj)],oldrate=40000,newrate=200)
   F,pp,cohe,Fx2y,Fy2x,Fxy=granger(outi,outj,20)#,200,100)# nfs2[0], nfs2[1] ,

#outout=downsample(time_courses[int(outdegree)],oldrate=40000,newrate=200)
#outin=downsample(time_courses[int(indegree)],oldrate=40000,newrate=200)
#synapse=downsample(recinsyn[int(indegree)],oldrate=40000,newrate=200)

#time_courses[int(outdegree)],recinsyn
###
#plt.hold(False) #Matlab style hold
sfin=str(fin)+'png' 
 fig.savefig(sfin) 
 fig=plt.figure(sfin)
fin=+1
###


bc=' between cell i= '+str(int(indexi)) +'and cell j= '+str(int(indexj))#+ 'ff= '+str(int(h.ff))+'prune net= '+str(int(h.prunenet))


# Plot Granger spectra
plt.hold(True) #Matlab style hold
labels=list()
labels=['Fx2y','Fy2x','Fxy','cohe','pp']
alldata=list()
colors=[[0,0,1],[1,0,0],[0,0.5,0.5],[0.5,1,0],[1,0.5,1]]
lb1='causality between cells i= '+str(int(indexi)) +' and cell j= '+str(int(indexj))
lb2='causality between cells i= '+str(int(indexj)) +' and cell j= '+str(int(indexi))
plt.plot(F,Fxy,label=lb1,linewidth=1.5,c=colors[0])
plt.plot(F,Fy2x,label=lb2,linewidth=1.5,c=colors[1])
#plt.xlim(0,2)
#plt.ylim(-0.4,0.4)
plt.xlabel("Hz")
plt.xlim(0,100)
plt.ylim(-2,2)
plt.title('GSC components'+bc)
plt.legend()#bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
#      ncol=2, mode="expand", borderaxespad=0.)
plt.hold(False) #Matlab style hold
fin=fin+1
###

sfin=str(fin)+'png' 
 fig.savefig(sfin) 
 fig=plt.figure(sfin)
###
plt.hold(True) #Matlab style hold
lb1='spectral power of cell '+str(int(indexi)) +' and cell j= '+str(int(indexj))
lb2='coherance between power spectrum of '+str(int(indexi)) +' and cell j= '+str(int(indexj))
plt.plot(F,pp,label=lb1,linewidth=1.5,c=colors[2])
plt.plot(F,cohe,label=lb2,linewidth=1.5,c=colors[3])


#plt.xlim(0,40)
#plt.ylim(0,0.5)
plt.xlim(0,100)
plt.ylim(0,1.5)
plt.ylabel("abs(Y)")
plt.xlabel("Hz")
plt.legend()
plt.title('GSC power and coherance'+bc)
#       ncol=2, mode="expand", borderaxespad=0.)
plt.hold(False) #Matlab style hold

fin=fin+1
###
sfin=str(fin)+'png' 
 fig.savefig(sfin) 
 fig=plt.figure(sfin)
F,pp,cohe,Fx2y,Fy2x,Fxy=granger(time_courses[int(indexi)],time_courses[int(indexj)],20)#,200,100)# nfs2[0], nfs2[1] ,

plt.plot(F,pp,linewidth=4,c=colors[4])
plt.title('Spectral power'+bc)
plt.xlim(0,100)
#plt.ylim(0,30)
plt.ylabel('log(abs(Y))')
plt.yscale('log')
plt.xlabel("Hz")


fin=fin+1
sfin=str(fin)+'png' 
 fig.savefig(sfin) 
 fig=plt.figure(sfin)
plt.hold(True) #Matlab style hold
#spectrum1=downsample(lfp2,oldrate=40000,newrate=400)
data = np.array(outi)
ps = np.abs(np.fft.fft(data))**2
time_step = h.dt#1 / 30
freqs = np.fft.fftfreq(data.size, 0.0025)
idx = np.argsort(freqs)


data = np.array(outj)
ps1 = np.abs(np.fft.fft(data))**2
time_step = h.dt#1 / 30
freqs1 = np.fft.fftfreq(data.size, 0.0025)
idx = np.argsort(freqs)


plt.title('downsampled memb potential FFT')
#f, axarr = plt.subplots(2, sharex=True)
lb3='cell'+str(j)
plt.plot(freqs1[idx], ps1[idx],linewidth=1,c=colors[2],label=lb3)
plt.title('downsampled memb potential FFT')





plt.plot(freqs[idx], ps[idx],linewidth=1,c=colors[1],label='cell'+str(i))
plt.ylabel('log(abs(Y))')
plt.yscale('log')
plt.xlim(0,100)
plt.xlabel('Hz')

fin=fin+1
sfin=str(fin)+'png' 
 fig.savefig(sfin) 
 fig=plt.figure(sfin)
plt.legend()#bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
    #   ncol=2, mode="expand", borderaxespad=0.)       
plt.plot(F,cohe,linewidth=1,c=colors[3],label='coherance')
plt.title('coherance between cell FFT')
plt.ylabel('log(abs(Y))')
plt.yscale('log')
plt.xlim(0,100)
plt.xlabel('Hz')
#plt.show()
plt.hold(False)

