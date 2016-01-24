

import matplotlib
matplotlib.use('Agg')# //Must be before importing matplotlib.pyplot or pylab!
import Image
import matplotlib.pyplot as plt

import pylab as p1
from pylab import *

import numpy as np
import scipy

import time; tic=time.clock()
from scipy import loadtxt, size, shape, zeros, mod, floor, mean
from pylab import figure, plot, xlabel, ylabel, legend, xlim, ylim, show, hold, squeeze, sqrt
from bsmart import granger # Load the Granger calculation tool

t = np.arange(0.0, 1.0+0.01, 0.01)


os.chdir(h.graph_dir)


"""
data = loadtxt("HAM_P5R1_pvsoma.dat")
datafft = abs(scipy.fft(data))
datafft2 = datafft[1:100]

p1.figure(figsize=(8, 6), dpi=80)
figure(1)
p1.plot(datafft2)
show()

import os
os.system('ls *.dat > swc_names.txt')
f = open('swc_names.txt')
nfs = [line.strip() for line in open('swc_names.txt', 'r')] # a python list comprehension
nfs2 = [line.strip() for line in open('swc_names.txt', 'r')] # a python list comprehension
#avg_sig = [line.strip() for line in open('swc_names.txt', 'r')]
avg_sig=0
string=""
array = []
subplotarg=910

for x in range(0,len(nfs)-1):


 nfs2[x] = loadtxt(nfs[x])
 np.array(nfs2[x])
 avg_sig=avg_sig+nfs2[x] # add all the signals up.
 """
#np.array=nfs2

 #datafft = abs(scipy.fft(nfs2[x]))
 #datafft2 = datafft[1:100]
 #subplotarg+=x
 #subplot(subplotarg)
 
 #plot(figsize=(8, 6), dpi=80)
 #plot(datafft2)
 #figure(x)
 #plot(nfs2[x])
#avg_sig=avg_sig/len(nfs2) #divide by the number of signals.
#figure(len(nfs)) 
#plot(avg_sig)
#show()
 
 #for y in range(0,len(nfs2[x])-1):
  #print nfs2[x], x, y
#print nfs2[0] 
#string = "%s%d" % ('names',x)
#print string
#print nfs[x]
#exec "string=nfs[x]"
#eval(string,"=nfs[x]")
#String = 'nfs[x] = ' + nfs[x]
#print String
#exec String
#nfs[x]
#show()

#import bsmart *

#matrix1 = zip(nfs2[0])  

#os.chdir('/home/zaza3/ncdemo/ncdemo/graphs')
os.system("pwd")
#print h.system("pwd") 
#print "pwd"


            
#import pdb
            
results=[] # Initialize list to store results            
print 'Calculating Granger...'
#F,pp,cohe,Fx2y,Fy2x,Fxy=granger( nfs2[0] , nfs2[1] , order=20) # Pick out LFPs from layers (default 2/3 and 5)
lfp1=[]
lfp2=[]

lfp1=h.l1lfp
lfp2=h.l2lfp
#outfile = TemporaryFile()
f1 = open('savelfp1', 'w')
np.save(f1, np.array(lfp1))
f2 = open('savelfp2', 'w')
np.save(f2, np.array(lfp2))

f1 = open('savelfp1', 'w')
np.save(f1, np.array(lfp1))

from scipy import stats
#F,pp,cohe,Fx2y,Fy2x,Fxy=granger( lfp1 , lfp2, order=20)# nfs2[0], nfs2[1] , order=20) # Pick out LFPs from layers (default 2/3 and 
    
    
F,pp,cohe,Fx2y,Fy2x,Fxy=granger( lfp1 ,lfp2, order=20)# nfs2[0], nfs2[1] ,

#pdb.set_trace()
#pdb.set_trace()

#results.append(Fx2y)

print "Plotting..."
#p1=figure(figsize=(10,6)) # Open figure and store its handle

###
figure(1)
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
plt.xlim(0,150)
plt.ylim(-0.10,1) # Need to define explicitly since otherwise plot tries to approximate infinity
weight_value=0
weight_value=h.i#h.w

plt.savefig('SGCofLFP.png')
plt.close

#End first SGC plot


xlabel('Frequency (Hz)')
ylabel('power')
title(str(weight_value)+' SGC of LFP ch1,ch2')
legend()

#show()
plt.hold(False)
print 'plotted GC'

F,pp,cohe,Fx2y,Fy2x,Fxy=granger(vectorslist[indegree], vectorslist[outdegree],order=20) 
#vectorslist[int(h.i)]


print "Plotting..."
#p1=figure(figsize=(10,6)) # Open figure and store its handle

###
figure(2)
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
plt.xlabel('Frequency (Hz)')
plt.ylabel('power')
plt.title(str(weight_value)+' SGC of neuron 24 and 9')
plt.legend()
#show()
plt.savefig('SGCofMembtimecourse.png')
plt.close


hold(False)
print 'plotted GC'



datafft =abs(scipy.fft(vectorslist[indegree]))# abs(scipy.fft(lfp1))
datafft2 =abs(scipy.fft(vectorslist[outdegree]))# abs(scipy.fft(lfp2))
datafft2 = datafft2[1:200]
datafft = datafft[1:200]

###
figure(3)
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

plt.savefig('FFT.png')
plt.close


#p1.savefig('FFT.png')
##
# The question that FFT spectrum answers is what kind of filter is 
# The neural network.
##



print 'plottedFFT'
 
# datafft =abs(scipy.fft(vectorslist[indegree]))# abs(scipy.fft(lfp1))
dfft = abs(scipy.fft(lfp1))# ,#h.record_soma))
datafft = dfft[1:100]
dfft2 = abs(scipy.fft(lfp2))# ,#h.record_soma))
datafft2 = dfft2[1:100]
#interactive(True)
#figure(figsize=(8, 6), dpi=80)
plt.figure(4)
plt.hold(True)
###
###
plt.plot(datafft2,label='FFT of LFP channel 2')
plt.plot(datafft,label='FFT of LFP channel 1')
plt.title('Fourier LFP chan1 and LFP chan2')

#show()
plt.hold(False)
plt.savefig('FFTofLFP2.png')
plt.close()
##
plt.figure(5)

dend_tvec = h.dend_tvec.to_python()
ones=[1]*len(dend_tvec)
plt.scatter(dend_tvec,ones)
#plt.plot(dend_tvec,ones)
plt.savefig('spikes_dendrite.png')
plt.close()

plt.figure(6)


ones=[1]*len(axon_tvec)
#plt.plot(dend_tvec,ones)
plt.scatter(axon_tvec,ones)
plt.savefig('spikes_axon.png')
plt.close()


#nrnpython("execfile('import_data.py')")
toc=time.clock()
print toc
print 'Done; elapsed time was %0.1f seconds.' % (toc-tic)


#start=size(alldata[0],0)/12 # Start at 5 Hz
#corrbd=stats.spearmanr(alldata[0][startp1:endp1],alldata[1][startp1:endp1]) # Correlation between baseline and damage
#corrbp=stats.spearmanr(alldata[0][startp1:endp1],alldata[2][startp1:endp1]) # Correlation between baseline and prosthesis
corrbd2=stats.spearmanr(lfp1 , lfp2)
#corrbd1=stats.spearmanr(lfp2 , lfp1)
corrbd1=stats.spearmanr(vectorslist[indegree] , vectorslist[outdegree])

os.chdir(h.workingdir)


"""
print corrbd1; print corrbd2

print "Plotting..."
figh=figure(figsize=(8,11)) # Open figure and store its handle
figh.subplots_adjust(left=0.16) # Less space on left
figh.subplots_adjust(right=0.98) # Less space on right
figh.subplots_adjust(top=0.95) # Less space on top
figh.subplots_adjust(bottom=0.06) # Less space on bottom
figh.subplots_adjust(wspace=0.25) # More space between
figh.subplots_adjust(hspace=0.2) # More space between
figtext(0.03,0.85,'A',size='xx-large')
figtext(0.03,0.62,'B',size='xx-large')
figtext(0.03,0.39,'C',size='xx-large')
figtext(0.03,0.16,'D',size='xx-large')


# Plot rasters
colors=array([[0.42,0.67,0.84],[0.50,0.80,1.00],[0.90,0.32,0.00],[0.34,0.67,0.67],[0.42,0.82,0.83],[0.90,0.59,0.00],[0.33,0.67,0.47],[0.42,0.83,0.59],[0.90,0.76,0.00],[1.00,0.85,0.00],[0.71,0.82,0.41],[0.57,0.67,0.33],[1.00,0.38,0.60]]) # Colors for each cell population
spc=0 # Only one column, so pick it
#for i in range(nfiles):
thisaxis=figh.add_subplot(4,1,1)
hold(True)
for j in range(len(vectorslist)):
   # if size(popdata[i][j][:,3]==spc):
       # thisdata=popdata[i][j][popdata[i][j][:,3]==spc,:] # Pull out spikes from column spc
        # Reorder spikes so the E2 is at the top of the plot and I6 is at the bottom
        scatter(dend_tvec[:,0]/1000.,axon_tvec[:,1],s=5,c=colors[j],edgecolors='none')
        xlim(0,maxtime) # To convert to seconds
        ylim(0,ncells+10) # Just larger than the number of cells in the model
        ylabel("Neuron number")
thisaxis.set_xticks(r_[0:maxtime+1]) # Don't show half seconds

# Plot LFPs
lfpcolors=[[0,0,1],[1,0,0],[0,0.5,0]] # Use colors from comparecausality
for i in range(nfiles):
    offset=20 # Don't show the first few points
    thisaxis=figh.add_subplot(4,1,4)
    hold(True)
    timeaxis=r_[0:maxtime*fs]/float(fs) # Set time axis
    npts=size(timeaxis,0)
    tmp=plot(timeaxis,lfpdata[i][0+offset:npts+offset,whichlfp,0],linewidth=1,c=lfpcolors[i]) # The middle index -- default 4 -- indicates which LFP to use
    xlim(0,maxtime)
    ylim(0.5,1.0) # Cortex: big voltage change
    if i==nfiles-1: xlabel('Time (s)')
    ylabel('Normalized voltage')
    thisaxis.set_xticks(r_[0:maxtime+1]) # Don't show half seconds


toc=time.clock()
print 'Done; elapsed time was %0.1f seconds.' % (toc-tic)
show()

for i in range(len(h.cells)):#.count)):
 for j in range(len(h.cells)):#
    print i 
    print j
    corrbd2=stats.spearmanr(vectorslist[i]  ,  vectorslist[j])
    #corrbd1=stats.spearmanr(lfp2 , lfp1)
"""    

