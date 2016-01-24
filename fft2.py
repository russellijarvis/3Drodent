#!/usr/bin/python
# -*- coding: utf-8 -*-
#import py_compile
#py_compile.compile('mymodule.py')
#import matplotlib

#matplotlib.use('Agg')# //Must be before importing matplotlib.pyplot or p1!
# import Image
# import matplotlib.pyplot as p1

# from __future__ import division


import matplotlib.pyplot as p1
import pylab as p1
from pylab import *

# import nrn
# import hoc
# h=hoc.HocObject()

from neuron import h
import numpy as np

# import scipy does not work.

import scipy as scipy
import time
tic = time.clock()
import pylab as p1

# from scipy import loadtxt, size, shape, zeros, mod, floor, mean

from pylab import figure, plot, xlabel, ylabel, legend, xlim, ylim, \
    show, hold, squeeze, sqrt
from bsmart import pwcausalr  # Load the Granger calculation tool

from numpy import sin, linspace, pi
from pylab import plot, show, title, xlabel, ylabel, subplot
from scipy import fft, arange

from scipy import fftpack
dt = 0.025
N=len(ecpv2)
#N=tstop/dt + 1 # N = 40.001*tstop
# dt sample spacing
#dt = 0.025
#f_s = 1000 / dt
# time vector
t = np.linspace(0.0, dt * N, N)


fin=0

# import np.ff
ecpv=h.ecpv.to_python()
ecpv2=h.ecpv2.to_python()

"""
Octave code,convert to Python
s1=x(1:50000);
M=500000
y = (fft(s1,length(s1)));
%clf(2)
M=1024
figure(2);%subplot(2,1,1);
%hold on
%ωk =2π/N*k
#f_k=(f_s/2*pi)
hold on;
plot((1:length(y)/16)*(Fs/length(y)), abs(y(1:length(y)/16))) ;%

s1=ecpv


dt=int(h.dt)
xf = np.linspace(0.0, dt * N, N / 2)
f = np.linspace(0.0, dt * N, N)

p1.hold(True)
p1.plot(np.array(splfp2d.to_python()))
p1.plot(np.array(splfp1d.to_python()))
p1.hold(False)

p1.title('power_spectrum_of_LFP_via_HOC')
p1.xlabel('Hz')
p1.ylabel('|H(S)|')
sfin = 'power_spectrum_of_LFP_via_HOC' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(numcell) \
    + str(int(h.minr)) + str(int(fin)) + '.png'
p1.xlim(0,200)
p1.yscale('log')
p1.show()
#p1.savefig(sfin, format='png')
"""

#ecpv=thresholded_arr(ecpv,50,10000,0)
#ecpv=thresholded_arr(ecpv,50,10000,0)
#ecpv2=thresholded_arr(ecpv2,-50,-10000,0)
#ecpv2=thresholded_arr(ecpv2,-50,-10000,0)

fin += 1
fig = p1.figure()
fig.clf()
p1.hold(True)
##
p1.plot(t, ecpv2)  # ,linewidth=1.5)#,labels='inhibitory')
p1.plot(t, ecpv)  # ,linewidth=1.5)#,labels='excitatory')
##
p1.hold(False)
p1.title('LFPs')
p1.xlabel('ms')
p1.ylabel('mV')
sfin = 'LFP' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
    + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
    + str(int(fin)) + '.png'
p1.savefig(sfin)
"""
"""
def get_lim_max(data, ind): 	# return the position of the maximal element in data located in the postion indexed by ind
		return  ind[data[ind].argmax()]

def calc_psd(self, sfin):
	"""
	Unoriginal Code from Neymotin 2013
	Ih tunes oscillations in an In Silico CA3 model (Neymotin et al. 2013)
	Accession: 151282
	from network.py
	http://senselab.med.yale.edu/ModelDB/ShowModel.asp?model=151282&file=\ca3ihdemo\network.py
	"""

		#self.calc_lfp()
	t0   = 200 # reject first ms of the signal
	fmax = 200 # upper limit for a periodogram frequency
	div  = int(1000/h.dt/(2*fmax)) # downsample the signal
	tr = [3,  12] # Theta frequency range
	gr = [30, 80] # Gamma frequency range
	t0i = int(t0/h.dt)
	if t0i > len(self):
		print "LFP is too short! (<200 ms)"
		return 0,0,0,0,0,0

	p1.figure()
	p1.clf()

	p1.subplot(2,1,1) # plot LFP
	p1.plot(np.array(range(len(self)))*h.dt, self)

	p1.subplot(2,1,2) # plot periodogram
	data = np.array(self[t0i::div]) # downsample data
	Pxx, freqs = p1.psd(data-data.mean(), Fs=1000/h.dt/div) # calculate FFT
	tind = np.where((freqs>=tr[0]) & (freqs<=tr[1]))[0] # index where for theta frequences  
	gind = np.where((freqs>=gr[0]) & (freqs<=gr[1]))[0] # index where for gamma frequences
	tp = Pxx[tind].mean() * np.diff(tr) # integral over theta power
	gp = Pxx[gind].mean() * np.diff(gr) # integral over gamma power
	ti = get_lim_max(Pxx, tind) # index of the frequency with a maximal power in theta range  
	gi = get_lim_max(Pxx, gind) # index of the frequency with a maximal power in gamma range
	tf = freqs[ti]
	gf = freqs[gi]
	p1.hold(True)
	p1.scatter(tf, 10*np.log10(Pxx[ti]), 100, 'b','o')
	#print Pxx[ti]
	p1.scatter(gf, 10*np.log10(Pxx[gi]), 100, 'r','o')
	p1.hold(False)
	p1.xlim(0,fmax) 
	p1.savefig(sfin)

sfin = 'PSD_of_LFP_1_including_gamma_and_theta' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
        + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
        + str(int(fin)) + '.png'

calc_psd(ecpv,sfin)
#calc_psd(np.array(h.ecpv.to_python()),sfin)
sfin = 'PSD_of_LFP_2_including_gamma_and_theta' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
        + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
        + str(int(fin)) + '.png'
calc_psd(ecpv2,sfin)
#calc_psd(h.ecpv2.to_python(),sfin)
sfin = 'PSD_of_av_mem_ex_including_gamma_and_theta' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
        + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
        + str(int(fin)) + '.png'
calc_psd(h.vsum1.to_python(),sfin)
sfin = 'PSD_of_av_mem_in_including_gamma_and_theta' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
        + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
        + str(int(fin)) + '.png'
calc_psd(h.vsum2.to_python(),sfin)


t0   = 200 # reject first ms of the signal
fmax = 200 # upper limit for a periodogram frequency
div  = int(1000/h.dt/(2*fmax)) # downsample the signal
tr = [3,  12] # Theta frequency range
gr = [30, 80] # Gamma frequency range
t0i = int(t0/h.dt)


p1.figure()
p1.clf()
		
# plot periodogram
data = np.array(h.ecpv.to_python()[t0i::div]) # downsample data
Pxx, freqs = p1.psd(data-data.mean(), Fs=1000/h.dt/div) # calculate FFT
tind = np.where((freqs>=tr[0]) & (freqs<=tr[1]))[0] # index where for theta frequences  
gind = np.where((freqs>=gr[0]) & (freqs<=gr[1]))[0] # index where for gamma frequences
tp = Pxx[tind].mean() * np.diff(tr) # integral over theta power
gp = Pxx[gind].mean() * np.diff(gr) # integral over gamma power
ti = get_lim_max(Pxx, tind) # index of the frequency with a maximal power in theta range  
gi = get_lim_max(Pxx, gind) # index of the frequency with a maximal power in gamma range
tf = freqs[ti]
gf = freqs[gi]
print gf, tf, freqs
p1.plot(freqs, 10*np.log10(Pxx))#, 100, 'b','o')
print Pxx[ti]
p1.plot(freqs, Pxx)#, 100, 'r','o')
p1.xlim(0,fmax)
sfin = 'FFT of LFP linear scale'  + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
        + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
        + str(int(fin)) + '.png'
p1.savefig(sfin)




    
vsum2 = np.array(h.vsum2.to_python())


vsum1 = np.array(h.vsum1.to_python())


fin += 1
fig = p1.figure()
fig.clf()
p1.hold(True)
p1.plot(t, vsum2)  # ,linewidth=1.5)#,labels='inhibitory')
p1.plot(t, vsum1)  # ,linewidth=1.5)#,labels='excitatory')
p1.xlim(0, 1000)
p1.hold(False)

# p1.plot(tc[0:25],out1[0:25])

p1.title('averaged membrane potentials')
p1.xlabel('ms')
p1.ylabel('mV')
sfin = 'averaged_membrane_potentials' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(numcell) \
    + str(int(h.minr)) + str(int(fin)) + '.png'
#p1.savefig(sfin, format='png')

"""

# import np.fft as fft

# Number of samplepoints
"""
N = 40001

# dt sample spacing

dt = 0.025
#f_s = 1000 / dt

# time vector

t = np.linspace(0.0, dt * N, N)

# Frequency vector

#from bsmart import armorf, spectrum_AR

#ntrls = 1
#npts = 40001
#p = 20
#fs = 40000  # is the sampling rate (e.g. 200 Hz)

# 200Hz, is used because of downsampling which should be thought of as frequency binnning.

#freq = fs / 2  # is the maximum frequency to calculate (e.g. fs/2=100, which will return 0:100 Hz)
#fmax = 200

# psd(x,40000,axis='0,200')

"""
"""
# rate=1/200

rate = 200

# psd(x, NFFT=256, Fs=2, Fc=0, detrend=mlab.detrend_none,
#     window=mlab.window_hanning, noverlap=0, pad_to=None,
#     sides='default', scale_by_freq=None, **kwargs)
# ------------------- Ploting Power spectra ---------------------
#
fin=0
ecpv=h.ecpv.to_python()
ecpv2=h.ecpv2.to_python()
lfp1 = ecpv
lfp2 = ecpv2


def plotSpectrum(y,Fs): #putting a variable as an argument makes it local.
	global fin #not doing, but listing it here, tells it to use the variable that was declared outside of the
 #procedure definition.

	#Acknoledgement, from the glowing Python http://glowingpython.blogspot.com.au/
	#Plots a Single-Sided Amplitude Spectrum of y(t)

	n = len(y) # length of the signal
	k = arange(n)
	T = n/Fs
	frq = k/T # two sides frequency range
	frq = frq[range(n/2)] # one side frequency range

	Y = fft(y)/n # fft computing and normalization
	Y = Y[range(n/2)]
	fin += 1
	fig = p1.figure()
	fig.clf() 
	p1.plot(frq,abs(Y),'r') # plotting the spectrum
 #p1.autoscale(True)
 #p1.yscale('log')
 #p1.xscale('log')
	
	p1.xlabel('Freq (Hz)')
	p1.ylabel('|Y(freq)|')

	sfin = 'naive fourier transform no scale' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
	    + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
	    + str(int(fin)) + '.png'
	p1.savefig(sfin, format='png')
plotSpectrum(ecpv,40000) 

 
 

"""
if sum(ecpv) != 0:
    spectrumecpv = downsample(ecpv, oldrate=40000, newrate=200)
    spectrumecpv2 = downsample(ecpv2, oldrate=40000, newrate=200)
    fin += 1
    fig = p1.figure()
    fig.clf()
    p1.hold(True)
    p1.plot(t, ecpv2)  # ,linewidth=1.5)#,labels='inhibitory')
    p1.plot(t, ecpv)  # ,linewidth=1.5)#,labels='excitatory')
    p1.hold(False)

 # p1.plot(tc[0:25],out1[0:25])

    p1.title('averaged membrane potential')
    p1.xlabel('ms')
    p1.ylabel('mV')
    sfin = 'psd' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
        + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
        + str(int(fin)) + '.png'
    p1.savefig(sfin, format='png')

vsum2 = np.array(h.vsum2.to_python())
spectrumvs2 = np.array(downsample(vsum2, oldrate=40000, newrate=200))

vsum1 = np.array(h.vsum1.to_python())
spectrumvs1 = np.array(downsample(vsum1, oldrate=40000, newrate=200))

execfile('pyhoc.py')
(  # ,200,100)# nfs2[0], nfs2[1] ,
    F,
    pp,
    cohe,
    Fx2y,
    Fy2x,
    Fxy,
    ) = granger(spectrumvs1, spectrumvs2, 1)

# F,pp1,cohe,Fx2y,Fy2x,Fxy=granger(vsum1,vsum2,1)#,200,100)# nfs2[0], nfs

colors = [[0, 0, 1], [1, 0, 0], [0, 0.5, 0.5], [0.5, 1, 0], [1, 0.5, 1]]

fin += 1
fig = plt.figure(fin)
plt.plot(F, pp, linewidth=1, c=colors[0])  # ,labels='pp')

# plt.plot(F,pp,label=lb2,linewidth=1.5,c=colors[1])
# plt.xlim(0,2)
# plt.ylim(-0.4,0.4)

plt.xlabel('Hz')
plt.xlim(0, 100)
p1.yscale('log')
p1.xscale('log')
# plt.ylim(0,20)

plt.title('Power Spectrum')

# plt.legend()#bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
#      ncol=2, mode="expand", borderaxespad=0.)

###

fig.savefig('power_spec_via_sgc')

fin += 1
fig = p1.figure()
fig.clf()
p1.hold(True)
p1.plot(t, vsum2)  # ,linewidth=1.5)#,labels='inhibitory')
p1.plot(t, vsum1)  # ,linewidth=1.5)#,labels='excitatory')
p1.xlim(0, 1000)
p1.hold(False)

# p1.plot(tc[0:25],out1[0:25])

p1.title('averaged membrane potentials')
p1.xlabel('ms')
p1.ylabel('mV')
sfin = 'averaged_membrane_potentials' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(numcell) \
    + str(int(h.minr)) + str(int(fin)) + '.png'
p1.savefig(sfin, format='png')

# p1.plot(tc[0:25],out1[0:25])

fin += 1
fig = p1.figure()
fig.clf()

# p1.hold(True)

p1.plot(t, np.array(ecpv2))  # ,linewidth=1.5)#,labels='inhibitory')

p1.ylim(np.min(ecpv2), np.max(ecpv2))
p1.xlim(0, 1000)

# p1.hold(False)
# p1.plot(tc[0:25],out1[0:25])

p1.title('local_field_potentials1')
p1.xlabel('ms')
p1.ylabel('mV')
sfin = 'local_field_potentials1' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
    + str(int(fin)) + '.png'
p1.savefig(sfin, format='png')

fin += 1
fig = p1.figure()
fig.clf()

# p1.hold(True)

p1.plot(t, np.array(ecpv))  # ,linewidth=1.5)#,labels='inhibitory')

p1.ylim(np.min(ecpv), np.max(ecpv))
p1.xlim(0, 1000)

# p1.hold(False)
# p1.plot(tc[0:25],out1[0:25])

p1.title('local_field_potentials2')
p1.xlabel('ms')
p1.ylabel('mV')
sfin = 'local_field_potentials2' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
    + str(int(fin)) + '.png'
p1.savefig(sfin, format='png')

fin += 1
fig = p1.figure()
fig.clf()
p1.hold(True)
p1.plot(np.array(lfp1s))  # ,linewidth=1.5)#,labels='inhibitory')

# p1.plot(np.array(vfr))#,linewidth=1.5)#,labels='excitatory')

p1.xlim(0, 1000)

# p1.ylim(0,900)

#p1.yscale('log')
#p1.xscale('log')
#p1.hold(False)

# p1.plot(tc[0:25],out1[0:25])

p1.title('power_spectrum_of_LFP_via_HOC')
p1.xlabel('Hz')
p1.ylabel('|H(S)|')
sfin = 'power_spectrum_of_LFP_via_HOC' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(numcell) \
    + str(int(h.minr)) + str(int(fin)) + '.png'
p1.savefig(sfin, format='png')

fin += 1
fig = p1.figure()
fig.clf()
p1.hold(True)

# p1.plot(np.array(lfp1s))#,linewidth=1.5)#,labels='inhibitory')

p1.plot(np.array(vfr))  # ,linewidth=1.5)#,labels='excitatory')
p1.xlim(0, 200)

Re: How to plot power spectrum

Postby ted » Fri Aug 19, 2011 10:59 am
spctrm generates results identical to the procedure of the same name described in Numerical Recipes in C (Press et al., Cambridge University Press), with Bartlett window and overlapping data segments. NEURON's spctrm differs in the following ways:
1. It takes only one argument--a vector of data sampled at regular intervals, which may contain any number of samples.
2. It handles zero padding and selects the number of data segments automatically. Consequently the number of frequency bins, and the number of data segments that are averaged to calculate the power density in these bins, are not subject to user control--see Relationship between sample length and results below.

Given that the sample interval is dt, the Nyquist frequency is 1/(2*dt).

numf=	8191 =splfp1d.size-1

Given that the number of frequency bins is numf ( == length of result vector), then the frequency bin centers are at i/(2*dt*numf) where i = 0..numf-1. That is, the result vector contains signal power in numf frequency bins centered at 0, 1/(2*dt*numf), . . (numf-1)/(2*dt*numf).

Relationship between sample length and results

I ran some tests and found that the data vector must contain at least 16 samples in order to obtain a result vector with at least two elements (power at 0 Hz (DC) and one nonzero frequency). Furthermore, the number of frequencies in the result vector is a monotonically increasing staircase function of the number of data samples in the argument vector. The relationship between argument length and result length is summarized by this table, which lists the argument l engths (number of samples or "# S") at which the result length (number of frequencies or "# F") increases:

Code: Select all
    # S    # F
       16      2  | 2^1
  +8   24      4  | 2^2
  +16  40      8
       72     16
      136     32
      264     64
      520    128
    1032    256
    40000

One last note: there are many different ways to analyze the spectral content of signals, each with its own particular strengths and weaknesses. Almost certainly some of these have been implemented with Python, and should be callable directly from hoc.

# p1.ylim(0,900)
p1.xscale('log')
p1.yscale('log')
p1.hold(False)


dt = 0.0025  # 0.1ms = 0.0001 sec
nextpow2 = 32768  # here 2000ms of length

# nextpow2=4096 # only 250ms of length
# nextpow2=len(ecpv)

if sum(ecpv) != 0:
    fin += 1
    p1.figure(fin)
    rate = 200  # since down sampling
    psd(spectrumecpv, nextpow2, rate)  # signal, number of sample points.

  # sampling frequency.

    sfin = 'psd' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
        + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
        + str(int(fin)) + '.png'
    p1.title('psd of ecpv')

  # fig.savefig(sfin)

    p1.savefig(sfin, format='png')

    fin += 1
    p1.figure(fin)
    psd(spectrumecpv2, nextpow2, rate)
    sfin = 'psd' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
        + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
        + str(int(fin)) + '.png'
    p1.title('psd of ecpv2')

  # fig.savefig(sfin)

    p1.savefig(sfin, format='png')

fin += 1
p1.figure(fin)
psd(np.array(spectrumvs1), nextpow2, rate)
sfin = 'psd' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
    + str(int(h.ff)) + str(numcell) + str(int(h.minr)) + str(int(fin)) \
    + '.png'
p1.title('psd of vsum1')
p1.savefig(sfin, format='png')

fin += 1
p1.figure(fin)
psd(np.array(spectrumvs2), nextpow2, rate)
sfin = 'psd' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
    + str(int(h.ff)) + str(numcell) + str(int(h.minr)) + str(int(fin)) \
    + '.png'
p1.title('psd of vsum2')
p1.savefig(sfin, format='png')

fin += 1
fig = p1.figure()
fig.clf()
if int(h.get_dist) == 1:
    across = np.arange(0, int(h.adjnqs.v[8].size), 1)
    dist = h.adjnqs.v[8].to_python()
    transfer = h.adjnqs.v[7].to_python()

  # p1.plot(across,transfer)

    p1.scatter(dist, distance)
    p1.title('arc length versus attenuation ratio')
    p1.xlabel('um')
    p1.ylabel('Ohm')
    sfin = 'arc length ' + str(int(h.prunenet)) + str(tstop) \
        + str(h.plastic) + str(int(h.ff)) + str(numcell) \
        + str(int(fin)) + '.png'
    p1.savefig(sfin, format='png')

print 'fourier sum of power:', sum(lfp1s), ' ', sum(vfr)

 # p1.xlim(100)
# p1.hold(False)
"""

