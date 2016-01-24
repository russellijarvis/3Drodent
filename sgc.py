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
from neuron import h
#import nrn
#import hoc
#h=hoc.HocObject()
lfp1=np.array(h.ecpv.to_python())
lfp2=np.array(h.ecpv2.to_python())

vsum1=np.array(h.vsum1.to_python())
vsum2=np.array(h.vsum2.to_python())

bc='ff= '+str(int(h.ff))+'prune net= '+str(int(h.prunenet))


# Number of samplepoints
N= n = 40001*int(1000)/1000 
#dt sample spacing
dt = 0.025 
f_s=1000/dt
#time vector

dtds=1/0.005
#tds = np.linspace(0.0, dtds*(N-1), (N-1)/200)

#Frequency vector
xf = np.linspace(0.0, dt*N, N/2)
f = np.linspace(0.0, dt*N, N)



sfin=str(fin)+'png' 
fig.savefig(sfin) 
fig=plt.figure(sfin)
fin+=1
#interactive(True)
#plt.plot((xf,yf1))
#plt.title("Raster Plot")
plt.hold(True)
colors=array([[0.42,0.67,0.84],[0.50,0.80,1.00],[0.90,0.32,0.00],[0.34,0.67,0.67],[0.42,0.82,0.83],[0.90,0.59,0.00],[0.33,0.67,0.47],[0.42,0.83,0.59],[0.90,0.76,0.00],[1.00,0.85,0.00],[0.71,0.82,0.41],[0.57,0.67,0.33],[1.00,0.38,0.60]]) # Colors for each cell population
j=len(colors)-1#12
plt.plot(tvec,intervec,'bo',label='inhibitory interneuron')#,c=colors[j], markeredgecolor = 'none')
j-=1
plt.plot(tvec,pyra,'g^')#o',c=colors[j], markeredgecolor = 'none')

plt.plot(tvec,zerovec,'g^', label='pyramidal cell')#,c=colors[j], markeredgecolor = 'none')
j-=1
plt.plot(record_inputs,vecin,'ro',linewidth=15, label='synapse input stim')#, markeredgecolor = 'none')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=9,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(1000)
#manrnpxtime=h.max(h.tvec)
plt.xlim(0,maxtime) # To convert to seconds
plt.ylim(-2,numcell) # Just larger than the number of cells in the model
plt.ylabel("Cell number")
plt.xlabel("spike time (ms)")
#Save fig hangs the program for some reason.
#plt.savefig("ff.png") # generates 'libpng error: zlib error' under nrniv
#plt.savefig("ff.eps")
#plt.axis.set_xticks(r_[0:maxtime+1]) # Don't show half seconds

plt.hold(False)
sfin=str(fin)+'png' 
fig.savefig(sfin) 
fig=plt.figure(sfin)

fig.savefig('RasterPlot'+'.png')

spc=0 # Only one column, so pick it

	
execfile('pyhoc.py')




sfin=str(fin)+'png' 
fig.savefig(sfin) 
fig=plt.figure(sfin)
fin+=1
##
plt.hold(True)
tc=np.array(time_courses[int(1)])
N=len(tc)
t = np.linspace(0.0, 0.025*N, N)
t=np.array(t)

tc=np.array(time_courses[int(10)])
str3='cell number= '+str(10)
plt.plot(t[0:7555],tc[0:7555],linewidth=1.5,label=str3)
#plt.plot(tc[0:25],out1[0:25])
plt.title('pyramidal neuron membrane potential')
plt.xlabel("ms")
plt.ylabel("mV")
plt.legend()

"""
for i in xrange(0,int(len(time_courses))):
 string='voltages'
 bc=' cells i= '+str(int(i)) +'ff= '+str(int(h.ff))+'prune net= '+str(int(h.prunenet))
 string=string+str(i)+bc 
 tc=np.array(time_courses[int(i)])
 plt.plot(t[500:1500],tc[500:1500],linewidth=3)
 #plt.plot(tc[0:25],out1[0:25])
 plt.title('cells index membrane potential')
 plt.xlabel("ms")
 plt.ylabel("mV")

 #out1=downsample(time_courses[int(i)],oldrate=40000,newrate=200)
 #downsampled.append(out1)





numcell=int(h.ncell)

plt.hold(False)
vsum1 = np.array(vsum1)
vsum2 = np.array(vsum2)
#vtotal=vtotal/2 #so its an average again.
in_trace = np.array(in_trace)
out_trace = np.array(in_trace)
in_trace=in_trace/numcell

"""

###
sfin=str(fin)+'png' 
fig.savefig(sfin) 
fig=plt.figure(sfin)
fin+=1
###
plt.title("Entropy Plot")
plt.hold(True)
j=len(colors)-1

across=np.arange(0,numcell,1)
#plt.plot(tvec,intervec,'bs')#plt.plot(tvec,pyra,'g^')
plt.plot(across,input_marke,'r-',linewidth=2,label='H(x) of synapse input stim')#,c=colors[j], markeredgecolor = 'none')
j-=1
plt.plot(across,pspke,'go', label='H(x) of cell num spike train')
plt.xlabel('cell number')
plt.ylabel('bits/sec')
plt.hold(False)
plt.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)



j=len(colors)-1

across=np.arange(0,numcell,1)
#plt.plot(tvec,intervec,'bs')#plt.plot(tvec,pyra,'g^')
divin=np.zeros(numcell)
#need to use numpy to remove inf
np.divide(input_marke,ratein,divin)
winfs=isinf(divin) # replace zeros with infinity
divin[winfs]=40
j-=1
divout=np.zeros(numcell)
np.divide(pspke,rates,divout)

plt.hold(True)

##
sfin=str(fin)+'png' 
fig.savefig(sfin) 
fig=plt.figure(sfin)
fin+=1
##
plt.plot(across,divin,'r-',linewidth=2,label='H(x) of synapse input stim')#,'b|',c=colors[j], markeredgecolor = 'none')


plt.plot(across,divout,'go',label='H(x) of cell num spike train')
plt.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)
plt.title("Entropy divided by rate")
plt.xlabel('cell number')
plt.ylabel('bits/sec')

###
plt.hold(False)
sfin=str(fin)+'png' 
fig.savefig(sfin) 
fig=plt.figure(sfin)
fin+=1
plt.hold(True)
###

plt.title("Lempel Ziv")
j=len(colors)-1

plt.plot(across,input_markl,'r-',linewidth=2,label='H(x) of synapse input stim')#,'|',c=colors[j],linewidth=5)
j-=1
plt.plot(across,pspkl,'go',label='H(x) of cell num spike train')
plt.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(h.tstop)

plt.xlabel("Neuron number")
plt.ylabel("Bits/sec")
plt.hold(False)
plt.title("Lempel Ziv Divided by rate")
plt.hold(True)
j=len(colors)-1

divin=np.zeros(numcell)
#need to use numpy to remove inf
np.divide(input_markl,ratein,divin) #dividing a large num by small number creates
#approaching infinitely large number.
winfs=isinf(divin) # replace zeros with infinity
divin[winfs]=40
j-=1
divout=np.zeros(numcell)
np.divide(pspkl,rates,divout)


plt.plot(across,divin,'r-',linewidth=2,label='H(x) of synapse input stim')#,'|',c=colors[j],linewidth=5)
j-=1
plt.plot(across,divout,'go',label='H(x) of cell num spike train')
plt.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(h.tstop)

plt.xlabel("Neuron number")
plt.ylabel("Bits/sec")

###
plt.hold(False)
sfin=str(fin)+'png' 
fig.savefig(sfin) 
fig=plt.figure(sfin)
fin+=1
#plt.hold(True)
# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(across,crlts,'go', markeredgecolor = 'none')
axarr[0].set_title('correlations between input and cell number')
axarr[1].plot(across,fhv,'go', markeredgecolor = 'none')
axarr[1].set_title('nTE between input and cell number')
plt.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)
plt.xlabel('cell number')


plt.hold(False)
maxtime=int(h.tstop)
#plt.xlabel("Neuron number")
#plt.ylabel("nTE")
#plt.show()
#axarr[1].xlabel('neuron number')
""" 
"""
plt.title("nTE and Correlations")
###
j=len(colors)-1
plt.plot(across,crlts,'o',c=colors[j], markeredgecolor = 'none')
j-=1
plt.plot(across,fhv,'o',c=colors[j], markeredgecolor = 'none')
"""



print "Plotting..."
F,pp,cohe,Fx2y,Fy2x,Fxy=granger( tc , recinsyn , order=20)# nfs2[0], nfs2[1] ,
###
sfin=str(fin)+'png' 
 fig.savefig(sfin) 
 fig=plt.figure(sfin) 
###

# Plot Granger spectra
plt.hold(True) #Matlab style hold
labels=list()
labels=['Fx2y','Fy2x','Fxy','cohe','pp']
alldata=list()
colors=[[0,0,1],[1,0,0],[0,0.5,0.5],[0.5,1,0],[1,0.5,1]]
F1=F
plt.plot(F,Fx2y,label='causality of channel X to channel Y',linewidth=2,c=colors[0])
plt.plot(F,Fy2x,label='causality of channel Y to channel X',linewidth=2,c=colors[1])

plt.plot(F,Fxy,label='causality of channel X to channel Y',linewidth=2,c=colors[2])
plt.plot(F,Fy2x,label='causality of channel Y to channel X',linewidth=2,c=colors[3])
plt.title('GSC components of average membrane potential compared with cell of high indegree')


#what if I had no such command about the axis would it set limits itself?
#plt.set_autoscale_on(True)
#plt.autoscale_view(True,True,True)
plt.ylim(-3,6)

#plt.ylim(-1,2)
#plt.xlim(0,2)
plt.hold(False) #Matlab style hold

plt.figure(6)
plt.hold(True) #Matlab style hold

plt.plot(F,pp,label='causality of channel X to channel Y',linewidth=2,c=colors[2])
plt.plot(F,cohe,label='coherance between two signals Y to channel X',linewidth=2,c=colors[3])
#plt.ylim(0,9)
#plt.xlim(0,4)
plt.hold(False) #Matlab style hold
plt.legend()
plt.title('GSC coherance and causality LFPs')
figure(9)
F,pp,cohe,Fx2y,Fy2x,Fxy=granger( in_trace , vsum1 , order=20)# nfs2[0], nfs2[1] ,
plt.hold(True) #Matlab style hold




vsum1[:]=vsum2[:]+in_trace[:]/int(numcell) #add back the contribution of indegree cell
vsum1[:]=vsum2[:]-out_trace[:]/int(numcell) #remove the contribution of outdegree cell



#vtotal=vtotal+in_trace/numcell#add in degree cells contribution back.

#vtotal=vtotal-out_trace/numcell # subtract out degrees contribution


#F,pp,cohe,Fx2y,Fy2x,Fxy=granger( vsum1 ,out_trace , order=20)# nfs2[0], nfs2[1] ,

###
plt.figure(7)
###

# Plot Granger spectra
plt.hold(True) #Matlab style hold
labels=list()
labels=['Fx2y','Fy2x','Fxy','cohe','pp']
alldata=list()
colors=[[0,0,1],[1,0,0],[0,0.5,0.5],[0.5,1,0],[1,0.5,1]]
F1=F
plt.plot(F,Fx2y,label='causality of channel X to channel Y',linewidth=2,c=colors[0])
plt.plot(F,Fy2x,label='causality of channel Y to channel X',linewidth=2,c=colors[1])

plt.plot(F,Fxy,label='causality of channel X to channel Y',linewidth=2,c=colors[2])
plt.plot(F,Fy2x,label='causality of channel Y to channel X',linewidth=2,c=colors[3])
plt.hold(False) #Matlab style hold
plt.legend()
plt.ylim(-10,10)
#plt.xlim(0,2)
plt.title('GSC components of averaged membrane with cell membrane time course of cell high out degree')



F is the frequency vector for the remaining quantities
pp is the spectral power
cohe is the coherence
Fx2y is the causality of channel X to channel Y
Fy2x is the causality of channel Y to channel X
Fxy is the "instantaneous" causality (cohe-Fx2y-Fy2x I think)


#F,pp,cohe,Fx2y,Fy2x,Fxy=granger( in_trace , vsum1 , order=20)# nfs2[0], nfs2[1] ,
plt.figure(8)
plt.hold(True) #Matlab style hold

plt.plot(F,pp,label='causality of channel X to channel Y',linewidth=2,c=colors[2])
plt.plot(F,cohe,label='coherance between two signals Y to channel X',linewidth=2,c=colors[3])
#plt.ylim(-0,9)
#plt.xlim(0,4)
plt.legend()
plt.hold(False) #Matlab style hold
plt.title('GSC components of coherance and causality for average membrane potential')

plt.show()


#plt.title('LFP signals hippocampus and neocortex')
#LFP is a sum of synchronous movement of APs down adjacent dendrites with the same orientation
#Filtered with a median filter to remove spikes.


plt.figure(11)
plt.hold(True)
#for i in xrange(0,int(len(time_courses))): 
#for i=0,len(time_courses):
tc=np.array(trains[int(0)])
print len(tc)
print len(t)
#plt.plot(t,tc)
plt.plot(vfr)

#plt.title('cells index thresholded')

#plt.title('LFP signals hippocampus and neocortex')
#LFP is a sum of synchronous movement of APs down adjacent dendrites with the same orientation
#Filtered with a median filter to remove spikes.

plt.hold(False)

plt.figure(2)
plt.hold(True)
#recinsynp=plt.plot(t,recinsyn)
plt.plot(recinsyn2)



#plt.plot(vsum1)
#plt.plot(vsum2)
plt.title('post synapse input onto cell[0]')

#plt.title('averaged membrane potentials hippocampus and neocortex')

plt.hold(True)






#datafft =abs(scipy.fft(vectorslist[indegree]))# abs(scipy.fft(lfp1))
#yf1 = abs(scipy.fft(recinsynp))# ,#h.record_soma))
#datafft = dfft[1:250]
#yf2 = abs(scipy.fft(tcp))# ,#h.record_soma))



#FFT = abs(scipy.fft(recinsyn))

#fy = np.fft.fft(recinsyn)


dt = 0.025
n = 40001
sf=np.loadtxt('samp_freq')
#sf=np.loadtxt('samp_freq')
components=[]
components=sf

figure(3)
amps=abs(scipy.fft(recinsyn))
#freqs2=np.loadtxt("freqs.np")
#freqs = scipy.fftpack.fftfreq(n, d=dt) # Frequencies associated with each samples
#plt.plot(sf,amps)

plt.plot(amps,f)
#plt.xlim(0,25)
#plt.ylim(0,35)

#plt.plot(sf)



#freqs = scipy.fftfreq(4001,dt)

figure(4)
#plt.subplot(211)
plt.plot(t, recinsyn)

#plt.subplot(3)
amps2=abs(scipy.fft(tc))
figure(5)
plt.plot(amps2,f)
#
#plt.plot(sf, abs(scipy.fft(tc)))
#plt.xlim(0,25)
#plt.ylim(0,35)

plt.show()

plt.figure(6)

#plt.plot(yf1)
#plt.hold(True)



plt.title("FFT power spectrum of the cell[0] memb pot")
#filtered=scipy.filters.median_filter(lfp1, size=4001, footprint=None, output=None, mode='reflect', cval=0.0, origin=0)

#plt.plot(filtered)
#plt.title('synapse g FFT, and membrane potential FFT')
#recinsyn.plot(recin,0.1)

plt.grid()
plt.hold(False)

#plt.plot(datafft)
#dfft2 = abs(scipy.fft(lfp2))# ,#h.record_soma))
#datafft2 = dfft2[1:250]
"""
"""
for i in xrange(0,int(len(time_courses))): 
 for j in xrange(int(len(time_courses)),0): 
  tci=np.array(trains[int(i)])
  tci=np.array(trains[int(j)])

  #print len(tc)
  #print len(t)
  F,pp,cohe,Fx2y,Fy2x,Fxy=granger( tcp , tci , order=20)# nfs2[0], nfs2[1] ,
  plt.plot(F,Fx2y,label='causality of channel X to channel Y',linewidth=2,c=colors[0])
  plt.plot(F,Fy2x,label='causality of channel Y to channel X',linewidth=2,c=colors[1])

  plt.plot(F,Fxy,label='causality of channel X to channel Y',linewidth=2,c=colors[2])
  plt.plot(F,Fy2x,label='causality of channel Y to channel X',linewidth=2,c=colors[3])
"""

"""
print "does not get here"
print os.getcwd()
print h.getcwd()
print matplotlib.use('Agg')
#print Agg
print plt
plt.savefig(str(change_iterator)+str(weight_value)+'SGCLFP.png')

#print "does not get here"
plt.hold(False)


#alldata.append(results)
plt.xlim(0,150)
plt.ylim(-0.10,1) # Need to define explicitly since otherwise plot tries to approximate infinity
weight_value=0
weight_value=h.i#h.w



plt.xlabel('Frequency (Hz)')
plt.ylabel('power')
run_iterator=str(h.run_iterator)
plt.title(str(run_iterator)+str(weight_value)+' SGC of LFP ch1,ch2')
plt.legend()

#plt.savefig(str(change_iterator)+str(weight_value)+'SGCofLFP.png')
plt.show()
plt.close
#End first SGC plot

plt.hold(False)
print 'plotted GC'

F,pp,cohe,Fx2y,Fy2x,Fxy=granger(vectorslist[indegree], vectorslist[outdegree],order=20) 
#vectorslist[int(h.i)]


print "Plotting..."


###
plt.figure(2)
###
# Plot Granger spectra
plt.hold(True) #Matlab style hold
labels=list()
labels=['Fx2y','Fy2x','Fxy','cohe','pp']
alldata=list()
colors=[[0,0,0.5],[0.5,0,0],[0,0.5,0.5],[0.5,0.5,0],[0.5,0.5,0.5]]
plt.plot(F,Fx2y,label='causality of channel X to channel Y',linewidth=2,c=colors[0])
plt.plot(F,Fy2x,label='causality of channel Y to channel X',linewidth=2,c=colors[1])
#plot(F,Fxy,label='instantaneous causalit',linewidth=2,c=colors[2])
#plot(F,cohe,label='coherance',linewidth=2,c=colors[3])
#plot(F,pp,label='spectral power',linewidth=2,c=colors[4])
alldata.append(results)
plt.xlim(0,50)
plt.ylim(0,4) # Need to define explicitly since otherwise plot tries to approximate infinity
weight_value=0
weight_value=h.w
change_iterator=h.change_iterator
plt.xlabel('Frequency (Hz)')
plt.ylabel('power')
plt.title(str(change_iterator)+str(weight_value)+' SGC of neuron 24 and 9')
plt.legend()
#show()
#plt.savefig(str(change_iterator)+str(weight_value)+'SGCofMembtimecourse.png')
plt.show()
plt.close


hold(False)
print 'plotted GC'



datafft =abs(scipy.fft(vectorslist[indegree]))# abs(scipy.fft(lfp1))
datafft2 =abs(scipy.fft(vectorslist[outdegree]))# abs(scipy.fft(lfp2))
datafft2 = datafft2[1:200]
datafft = datafft[1:200]

###
plt.figure(3)
###
hold(True)

plt.plot(datafft,label='indegree cell')
plt.plot(datafft2,label='outdegree cell')
#The MATLAB convention is to counter intuitively supply labels after the plot command is given.
plt.xlabel('Frequency (Hz)')
plt.ylabel('power')
plt.title('Fourier transform of high indgree and high outdegree cell')
plt.show()
plt.hold(False)

#plt.savefig(str(change_iterator)+str(weight_value)+'FFT.png')
plt.close()


##
# The question that FFT spectrum answers is what kind of filter is 
# The neural network.
##



print 'plottedFFT'
""" 
"""
###
###
plt.plot(datafft2,label='FFT of LFP channel 2')
plt.plot(datafft,label='FFT of LFP channel 1')
plt.title('Fourier LFP chan1 and LFP chan2')

#show()
plt.hold(False)
#plt.savefig(str(change_iterator)+str(weight_value)+'FFTofLFP2.png')
plt.close()
plt.show()
##
plt.figure(5)

# Prefer to make these plots before binning of the vectors.

dend_tvec = h.dend_tvec.to_python()
ones=[1]*len(dend_tvec)
plt.scatter(dend_tvec,ones)
plt.show()
#plt.savefig(str(change_iterator)+str(weight_value)+'spikes_dendrite.png')
plt.close()

plt.figure(6)


ones=[1]*len(axon_tvec)

plt.scatter(axon_tvec,ones)
plt.show()
#plt.savefig(str(change_iterator)+str(weight_value)+'spikes_axon.png')
plt.close()

toc=time.clock()
print toc
print 'Done; elapsed time was %0.1f seconds.' % (toc-tic)


corrbd2=stats.spearmanr(lfp1 , lfp2)
corrbd1=stats.spearmanr(vectorslist[indegree] , vectorslist[outdegree])


plt.figure(7)
im=imshow(Matrix2)
plt.colorbar(im)
plt.xlabel('columns = targets')
plt.ylabel('rows = sources')
plt.title('Adjacency matrix')
#plt.savefig('connection_matrix3.png')  

plt.savefig(str(change_iterator)+str(weight_value)+'adjacency_matrix.png')
plt.close()
  
"""

#t = np.arange(0.0, 1.0+0.01, 0.01)# declare a 3d np array. for plotting time samples.


#os.chdir(h.results_dir)
#os.chdir('results') changingdir should be unncessary. Time_courses are already imported.
#time_courses


#These 4 lines convert from NEURON vector to python list, to np array.
"""

becpv=np.array(h.becpv.to_python())
becpv2=np.array(h.becpv2.to_python())

bvsum1=np.array(h.bvsum1.to_python())
bvsum2=np.array(h.bvsum2.to_python())

plt.figure(1)
plt.hold(True)
dfft = abs(scipy.fft(becpv))# ,#h.record_soma))
datafft = dfft[1:250]
plt.plot(datafft)

dfft = abs(scipy.fft(becpv2))# ,#h.record_soma))
datafft = dfft[1:250]
plt.plot(datafft)


#plt.plot(becpv2)
plt.title('LFP signals hippocampus and neocortex before median filter')
#LFP is a sum of synchronous movement of APs down adjacent dendrites with the same orientation
#Filtered with a median filter to remove spikes.
plt.hold(False)

plt.figure(2)
plt.hold(True)

dfft = abs(scipy.fft(bvsum1))# ,#h.record_soma))
datafft = dfft[1:250]
plt.plot(datafft)
dfft = abs(scipy.fft(bvsum2))# ,#h.record_soma))
datafft = dfft[1:250]
plt.plot(datafft)

plt.title('averaged membrane potential signals hippocampus and neocortex before median filter')
#LFP is a sum of synchronous movement of APs down adjacent dendrites with the same orientation
#Filtered with a median filter to remove spikes.
plt.hold(False)

"""
#print lfp1
#print lfp2
#os.chdir(h.workingdir)

"""
#if vsum1!=0:
plt.figure(6)
plt.hold(True)
#datafft =abs(scipy.fft(vectorslist[indegree]))# abs(scipy.fft(lfp1))
dfft = abs(scipy.fft(vsum1))# ,#h.record_soma))
datafft = dfft[1:250]
plt.plot(datafft)
dfft2 = abs(scipy.fft(vsum2))# ,#h.record_soma))
datafft2 = dfft2[1:250]
#interactive(True)
#figure(figsize=(8, 6), dpi=80)
plt.plot(datafft2)
plt.title('FFT power spectrum average membrane polential')
plt.hold(False)

print lfp1 
results=[] # Initialize list to store results            
print 'Calculating Granger...'

lfp1=[]
lfp2=[]

#outfile = TemporaryFile()

change_iterator=h.run_iter
weight_value=h.run_iter

f1 = open(str(change_iterator)+str(weight_value)+'savelfp1', 'w')
np.save(f1, np.array(lfp1))

f2 = open(str(change_iterator)+str(weight_value)+'savelfp2', 'w')
np.save(f2, np.array(lfp2))

#I may as well send these files to matlab here too.


#f1 = open(str(change_iterator)+str(weight_value)+'savelfp1', 'w')
np.save(f1, np.array(lfp1))

from scipy import stats
    
#print len(lfp2), len(lfp1)
print lfp1
print lfp2    
""" 


