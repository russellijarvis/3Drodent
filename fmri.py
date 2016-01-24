import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import csv2rec

import nitime
import nitime.algorithms as alg

import nitime.analysis as nta
import nitime.timeseries as ts
import nitime.utils as tsu
from nitime.viz import drawmatrix_channels

from nitime.viz import draw_graph

"""numcell=ncell
"roi_names=[0 for x in xrange(0,int(numcell))]
"for i in xrange(0,int(numcell-1)):
" roi_names[i]=str(i)+' '+str(h.cells.o(i).nametype)
" print roi_names[i]"""


TR = h.dt
f_ub = 150 #highest frequency
f_lb = 0.02 #lowest frequency.

roi_names = np.array(roi_names)
nseq = len(roi_names)-1

n_samples = int(len(time_courses[int(indexi)]))
data = np.zeros((nseq, n_samples))

for n_idx in xrange(0,int(len(time_courses)-1)):
    data[n_idx] = time_courses[n_idx]
    #if(np.isnan(data[n_idx])): 
    #   data[n_idx]=0
#We normalize the data in each of the ROIs to be in units of % change and initialize the TimeSeries object:
#figure(1)
pdata = tsu.percent_change(data)

#NormalizationAnalyzer(data.fir).z_score.data[0]
time_series = ts.TimeSeries(data, sampling_interval=TR)




#p1.plt(pdata)
#if(np.isnan(pdata)): print pdata

## shold I try to plot the time series to see what is wrong with it?
##Does this data contain NANs?


#We initialize the GrangerAnalyzer object, while specifying the order of the autoregressive model to be 1 (predict the current behavior of the time-series based on one time-point back).


#w, f_x2y, f_y2x, f_xy, Sw = alg.granger_causality_xy(time_series)#, ecov,                 n_freqs=n_freqs)

G = nta.GrangerAnalyzer(time_series, order=1) #number of lags should be two 5ms sample width lags.



#For comparison, we also initialize a CoherenceAnalyzer and a CorrelationAnalyzer, with the same TimeSeries object

C1 = nta.CoherenceAnalyzer(time_series)
C2 = nta.CorrelationAnalyzer(time_series)
#We are only interested in the physiologically relevant frequency band (approximately 0.02 to 0.15 Hz).

#The spectral resolution is different in these two different analyzers. In the CoherenceAnalyzer, the spectral resolution depends on the size of the window used for calculating the spectral density and cross-spectrum, whereas in the GrangerAnalyzer it is derived, as determined by the user, from the MAR model used.

#For this reason, the indices used to access the relevant part of the spectrum will be different in the different analyzers.

#freq_idx_G = np.where((G.frequencies > f_lb) * (G.frequencies < f_ub))[0]
freq_idx_C = np.where((C1.frequencies > f_lb) * (C1.frequencies < f_ub))[0]

coh = np.mean(C1.coherence[:, :, freq_idx_C], -1)  # Averaging on the last dimension
#g1 = np.m

#p1.show()
cohgraph=nx.to_networkx_graph(coh,create_using=nx.DiGraph())

fig03 = drawmatrix_channels(coh, roi_names, size=[10., 10.], color_anchor=0)
fig03.savefig(str(h.plastic)+str(h.ff)+str(tstop)+str(ncell)+str(h.prunenet)+str(fin)+'coherance.png')

fig05 = drawmatrix_channels(C2.corrcoef, roi_names, size=[10., 10.], color_anchor=0)
fig05.savefig(str(h.prunenet)+str(fin)+'correlation.png')

import pyentropy as pe
from pyentropy import DiscreteSystem
#These functions need a package not yet on Trifid.
#The functions need to be developed to include time binning.
tbinned=np.array(tbinned).astype(int)
trains=np.array(trains).astype(int)
def get_entropies():
  econt=[0 for x in xrange(0,int(num_cells))]
  econtub=[0 for x in xrange(0,int(num_cells))]
  mcont=[[0 for x in xrange(0,int(num_cells))] for x in xrange(0,int(num_cells))]
  econtb=[0 for x in xrange(0,int(num_cells))]

  for k in xrange(0,len(trains)):
    #quantised=pe.quantise(trains[k],3, uniform='sampling', minmax=None, centers=True)
    #10/0.025=10ms/f_{s}
    quantised=pe.quantise_discrete(trains[k],2) #np.array(quantised).astype(int)
    #pe.quantise_discrete(trains[k],10)
    sentbinned=DiscreteSystem(quantised,(1,2),quantised,(1,2) )
    sentbinned.calculate_entropies(method='plugin', calc=['HX'])
    econtb.append(sentbinned.H['HX'])          
    trains[k]=[int(i) for i in trains[k]]
    sentunbinned=DiscreteSystem(np.array(trains[k]),(1,2), np.array(trains[k]), (1,2))
    sentunbinned.calculate_entropies(method='plugin', calc=['HX'])
    econtub.append(sentunbinned.H['HX'])
    sent=DiscreteSystem((tbinned[k]),(1,5), (tbinned[k]), (1,5))
    sent.calculate_entropies(method='plugin', calc=['HX'])
    econt.append(sent.H['HX'])  #append only for every new k.   
    print sent.H['HX']," ",pspke[k]," ", k," ", sentbinned.H['HX'], sentunbinned.H['HX'],
    print sum(tbinned[k]), sum(quantised), "difference in spk count"
 
   
  
  across2=np.arange(0,int(ncell)*2+1,1)
  print pspke
  print econt
  print pspkl 
  print input_marke[0], econt[0]  


  return econt
  
get_entropies()
#sent.maxent() 
def get_MI():
  for k in xrange(0,len(trains)):
    for j in xrange(0,len(trains)):
      if(k!=j):
        trains[j]=[int(i) for i in trains[j]]
        trains[k]=[int(i) for i in trains[k]]
        if((sum(trains[j])!=0)|(sum(trains[k])!=0)):
          #if((sum(trains[j]))!=0)||(sum(trains[k])!=0)):
          #if(sum(trains[k])!=0):
          sent=DiscreteSystem((tbinned[j]),(1,5), (tbinned[k]), (1,5))
          sent=DiscreteSystem(np.array(trains[j]),(1,2), np.array(trains[k]), (1,2))
          #sent.calculate_entropies(method='plugin', calc=['HX'])
          #print sent.H, j, k 
          sent.calculate_entropies(method='plugin', calc=['HX', 'HXY'])
          #print sent.I(),j, k
      mcont[k][j]=sent.I['HXY']
      #sent.H['HX']
  return mcont  




print 'End the program the end!'

#quit()

#End the program the end!
#

#fig04 = drawgraph_channels(coh)
#fig04.savefig(str(h.prunenet)+str(fin)+'graphcoherance.png')



"""
def dB(x, out=None):
    if out is None:
        return 10 * np.log10(x)
    else:
        np.log10(x, out)
        np.multiply(out, 10, out)

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
import scipy.stats.distributions as dist

import nitime.algorithms as tsa
import nitime.utils as utils
from nitime.viz import winspect
from nitime.viz import plot_spectral_estimate
import nitime.algorithms as tsa
#f, adaptive_psd_mt, nu = tsa.multi_taper_psd(time_courses[j],  adaptive=True, jackknife=False)
f, adaptive_psd_mt, nu = tsa.multi_taper_psd(np.array(time_courses[j]),  adaptive=True, jackknife=False)


dB(adaptive_psd_mt, adaptive_psd_mt)


p975 = dist.chi2.ppf(.975, nu)
p025 = dist.chi2.ppf(.025, nu)

l1 = ln2db * np.log(nu / p975)
l2 = ln2db * np.log(nu / p025)

hyp_limits = (adaptive_psd_mt + l1, adaptive_psd_mt + l2)

fig06 = plot_spectral_estimate(freqs, psd, (adaptive_psd_mt,), hyp_limits,
                       elabels=('MT with adaptive weighting and 95% interval',))

ean(G.causality_xy[:, :, :,-1])#freq_idx_G], -1)
g1 = np.mean(G.causality_xy[:, :, freq_idx_G], -1)
g2 = np.mean(G.causality_xy[:, :, freq_idx_G] - G.causality_yx[:, :, freq_idx_G], -1)




fin+1
fig=p1.figure(fin)
fig.clf()
im=p1.imshow(g2,interpolation='nearest') 
p1.colorbar(im) 


fin+1
fig=p1.figure(fin)
fig.clf()
im=p1.imshow(g1,interpolation='nearest') 
p1.colorbar(im) 

fig01 = drawmatrix_channels(g1, roi_names, size=[10., 10.], color_anchor=0)
fig02 = drawmatrix_channels(g2, roi_names, size=[10., 10.], color_anchor=0)


#p1.show()


#fig03 = drawmatrix_channels(g1, roi_names, size=[10., 10.], color_anchor=0)
#fig02.savefig('grangerstuff.png')
"""
