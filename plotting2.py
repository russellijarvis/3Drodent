import nitime
import nitime.analysis as nta
import nitime.timeseries as ts
import nitime.utils as tsu
from nitime.viz import drawmatrix_channels
from nitime.viz import drawgraph_channels

import statsmodels as statsmodels
from matplotlib.colors import LogNorm 
import matplotlib 
#matplotlib.use('Agg') 
import pylab as p1
from scipy import signal 
from bsmart import granger # Load the Granger calculation tool
#and user make sure that param6.hoc has run.
# sfin=str(h.prunenet)+str(fin)+'.png' 
#fig.savefig(sfin) 


ioff()

fin=0 

fig=p1.figure(fin)

#interactive(True)
#p1.plot((xf,yf1))
#p1.title("Raster Plot")
p1.hold(True)
colors=array([[0.42,0.67,0.84],[0.50,0.80,1.00],[0.90,0.32,0.00],[0.34,0.67,0.67],[0.42,0.82,0.83],[0.90,0.59,0.00],[0.33,0.67,0.47],[0.42,0.83,0.59],[0.90,0.76,0.00],[1.00,0.85,0.00],[0.71,0.82,0.41],[0.57,0.67,0.33],[1.00,0.38,0.60]]) # Colors for each cell population
j=len(colors)-1#12
p1.plot(tvec,intervec,'bo',label='inhibitory interneuron')#,c=colors[j], markeredgecolor = 'none')
j-=1
p1.plot(tvec,pyra,'g^')#o',c=colors[j], markeredgecolor = 'none')

p1.plot(tvec,zerovec,'g^', label='pyramidal cell')#,c=colors[j], markeredgecolor = 'none')
j-=1
p1.plot(record_inputs,vecin,'ro',linewidth=15, label='synapse input stim')#, markeredgecolor = 'none')
p1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=9,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(1000)
#manrnpxtime=h.max(h.tvec)
p1.xlim(0,maxtime) # To convert to seconds
p1.ylim(-2,int(numcell)) # Just larger than the number of cells in the model
p1.ylabel("Cell number")
p1.xlabel("spike time (ms)")
#Save fig hangs the program for some reason.
#p1.savefig("ff.png") # generates 'libpng error: zlib error' under nrniv
#p1.savefig("ff.eps")
#p1.axis.set_xticks(r_[0:maxtime+1]) # Don't show half seconds

p1.hold(False)
sfin=str(h.prunenet)+str(fin)+'.png' 
draw()
fig.savefig(sfin)   
fin+=1
spc=0 # Only one column, so pick it

	
execfile('pyhoc.py')




sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin)
fin+=1
##
p1.hold(True)
tc=np.array(time_courses[int(1)])
N=len(tc)
t = np.linspace(0.0, 0.025*N, N)
t=np.array(t)

tc=np.array(time_courses[int(10)])
str3='cell number= '+str(10)
p1.plot(t[0:7555],tc[0:7555],linewidth=1.5,label=str3)
#p1.plot(tc[0:25],out1[0:25])
p1.title('pyramidal neuron membrane potential')
p1.xlabel("ms")
p1.ylabel("mV")
p1.legend()

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin)
fin+=1


for i in xrange(0,int(len(time_courses))):
 string='voltages'
 bc=' cells i= '+str(int(i)) +'ff= '+str(int(h.ff))+'prune net= '+str(int(h.prunenet))
 string=string+str(i)+bc 
 tc=np.array(time_courses[int(i)])
 p1.plot(t[500:1500],tc[500:1500],linewidth=3)
 #p1.plot(tc[0:25],out1[0:25])
 p1.title('cells index membrane potential')
 p1.xlabel("ms")
 p1.ylabel("mV")

 #out1=downsample(time_courses[int(i)],oldrate=40000,newrate=200)
 #downsampled.append(out1)





numcell=int(h.ncell)
"""
p1.hold(False)
vsum1 = np.array(vsum1)
vsum2 = np.array(vsum2)
#vtotal=vtotal/2 #so its an average again.
in_trace = np.array(in_trace)
out_trace = np.array(in_trace)
in_trace=in_trace/int(numcell)

"""

###
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
#belongs to above
fig=p1.figure(fin)
fin+=1
###
p1.title("Entropy Plot")
p1.hold(True)
j=len(colors)-1

across=np.arange(0,int(numcell),1)
#p1.plot(tvec,intervec,'bs')#p1.plot(tvec,pyra,'g^')
p1.plot(across,input_marke,'r-',linewidth=2,label='H(x) of synapse input stim')#,c=colors[j], markeredgecolor = 'none')
j-=1
p1.plot(across,pspke,'go', label='H(x) of cell num spike train')
p1.xlabel('cell number')
p1.ylabel('bits/sec')
p1.hold(False)
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)



j=len(colors)-1

across=np.arange(0,int(numcell),1)
#p1.plot(tvec,intervec,'bs')#p1.plot(tvec,pyra,'g^')
divin=np.zeros(int(numcell))
#need to use numpy to remove inf
np.divide(input_marke,ratein,divin)
winfs=isinf(divin) # replace zeros with infinity
divin[winfs]=40
j-=1
divout=np.zeros(int(numcell))
np.divide(pspke,rates,divout)

p1.hold(True)

##
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin)
fin+=1
##
p1.plot(across,divin,'r-',linewidth=2,label='H(x) of synapse input stim')#,'b|',c=colors[j], markeredgecolor = 'none')

p1.plot(input_marke,ratein,'go',label='H(x) of cell num spike train')

#p1.plot(across,divout,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)
p1.title("Entropy divided by rate")
p1.xlabel('cell number')
p1.ylabel('bits/sec')

###
p1.hold(False)
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin)
fin+=1
p1.hold(True)
###

p1.title("Lempel Ziv")
j=len(colors)-1

p1.plot(across,input_markl,'r-',linewidth=2,label='H(x) of synapse input stim')#,'|',c=colors[j],linewidth=5)
j-=1
p1.plot(across,pspkl,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(h.tstop)

p1.xlabel("Neuron number")
p1.ylabel("Bits/sec")
p1.hold(False)
p1.title("Lempel Ziv Divided by rate")
p1.hold(True)
j=len(colors)-1

divin=np.zeros(int(numcell))
#need to use numpy to remove inf
np.divide(input_markl,ratein,divin) #dividing a large num by small number creates
#approaching infinitely large number.
winfs=isinf(divin) # replace zeros with infinity
divin[winfs]=40
j-=1
divout=np.zeros(int(numcell))
np.divide(pspkl,rates,divout)


p1.plot(across,divin,'r-',linewidth=2,label='H(x) of synapse input stim')#,'|',c=colors[j],linewidth=5)
j-=1
p1.plot(divout,rates,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(h.tstop)

p1.xlabel("Neuron number")
p1.ylabel("Bits/sec")

###
p1.hold(False)
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fin+=1
fig=p1.figure(fin)
#p1.hold(True)
# Two subplots, the axes array is 1-d
f, axarr = p1.subplots(2, sharex=True)
axarr[0].plot(across,crlts,'go', markeredgecolor = 'none')
axarr[0].set_title('correlations between input and cell number')
axarr[1].plot(across,fhv,'go', markeredgecolor = 'none')
axarr[1].set_title('nTE between input and cell number')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)
p1.xlabel('cell number')


p1.hold(False)
maxtime=int(h.tstop)
#p1.xlabel("Neuron number")
#p1.ylabel("nTE")
#p1.show()
#axarr[1].xlabel('neuron number')
""" 
"""
p1.title("nTE and Correlations")
###
j=len(colors)-1
p1.plot(across,crlts,'o',c=colors[j], markeredgecolor = 'none')
j-=1
p1.plot(across,fhv,'o',c=colors[j], markeredgecolor = 'none')

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)
fin=fin+1
fig=p1.figure(sfin) 

#import os
#import numpy as np
#import matplotlib.pyplot as p1
#from matplotlib.mlab import csv2rec

#We then define a few parameters of the data: the TR and the bounds on the frequency band of interest.

TR = float(h.dt)#already is float, just acknowledgement.
f_ub = 100
f_lb = 0.02
#roi_names = []#np.array(allrows)
#We read in the resting state fMRI data into a recarray from a csv file:
#data_path = os.path.join(nitime.__path__[0], 'data')
#data_rec = csv2rec(os.path.join(data_path, 'fmri_timeseries.csv'))
#  storename= str(s[3])
#for i in xrange(1,int(len(allrows))):
# s=allrows[i]
# if(int(len(s))==9):
#  print i
#  roi_names[i] =  s[5]+str(i)


n_samples = int(tstop/h.dt)+2 #data_rec.shape[0]
#nseq=int(len(allrows))-1
nseq=int(numcell)
#data = np.zeros((int(numcell), int(numcell)))


roi_names=[0 for x in xrange(0,int(numcell)-1)]
for i in xrange(0,int(numcell)-1):
 roi_names[i]=str(h.cells.o(i).nametype)
 print roi_names[i]
#Make an empty container for the data
msgcv= [[0 for x in xrange(0,int(int(int(numcell))))] for x in xrange(0,int(int(int(numcell))))]
execfile('pyhoc.py')

for indexi in xrange(0,int(numcell)):
 for indexj in xrange(0,int(numcell)):
   outi=downsample(time_courses[int(indexi)],oldrate=40000,newrate=200)
   outj=downsample(time_courses[int(indexj)],oldrate=40000,newrate=200)
   #outi=percent_change(outi)
   #outi=percent_change(outj)
   F,pp,cohe,Fx2y,Fy2x,Fxy=granger(outi,outj,20)#,200,100)# nfs2[0], nfs2[1] ,
    #the two refers to a projection of x onto y, not an order of 2. Fill a matrix with
        #These values. 
    
        
   msgcv[indexi][indexj]=np.mean(Fx2y)
   
   
   np.shape(Fx2y)
 print 'pacificier', indexi 
msgcv=np.array(msgcv)
 	#Fx2y is the causality of channel X to channel Y
  	#Fy2x is the causality of channel Y to channel X
  	#Fxy is the "instantaneous" causality (cohe-Fx2y-Fy2x I think)
fig01 = drawmatrix_channels(msgcv, roi_names, size= [10., 10.], color_anchor=0)
sfin='granger_matrix1'+str(h.prunenet)+str(fin)+'.png' 
fig01.savefig(sfin)
fin=fin+1




#low_values_indices = np.mean(msgcv)+np.std(msgcv) < 0#///  # Where values are low
#msgcv[low_values_indices] = 0#//  # All low values set to 0
 
fig02 = drawgraph_channels(msgcv, roi_names)
sfin='grange+graph'+str(h.prunenet)+str(fin)+'.png' 
fig02.savefig(sfin)

fig=plt.imshow(msgcv, interpolation='nearest')
sfin='granger_matrix2'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)

  
msgcv=np.array(msgcv) 
low_values_indices = np.mean(msgcv)+np.std(msgcv) < 0#///  # Where values are low
msgcv[low_values_indices] = 0#//  # All low values set to 0

dirfinal =nx.MultiDiGraph() 
dirfinal.add_nodes_from(entfa) 
for i in xrange(0,len(entfa)):
 for j in xrange(0,len(entfa)):
  if entfa[i,j]>(np.mean(entfa)+np.std(entfa)):
    dirfinal.add_edge(i,j,weight=entfa[i,j]) 

fin+=1
fig = drawgraph_channels(dirfinal, roi_names)
sfin='grangergraph'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)



dirfinals =nx.MultiDiGraph() 
dirfinals.add_nodes_from(msgcv) 
for i in xrange(0,len(msgcv)):
 for j in xrange(0,len(msgcv)):
  if msgcv[i,j]>(np.mean(msgcv)+np.std(msgcv)):
    dirfinals.add_edge(i,j,weight=msgcv[i,j]) 

fin+=1
fig = drawgraph_channels(firfinals, roi_names)
sfin='grange+graph'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)

"""


fig01 = drawmatrix_channels(coh, roi_names, size=[10., 10.], color_anchor=0)
"""

fig01.savefig('sgc_matrix..png')

mbin2=np.array(mbin) 
mbin2=nx.to_networkx_graph(mbin2,create_using=nx.DiGraph())   #directed graph. 




"""		 
im=p1.imshow(WeightsAb,interpolation='nearest') 

p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Covariance of weight before training') 
# p1.autoscale(True) 
p1.grid(True)  """

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1

"""
im=p1.imshow(nTEin,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('The difference between weights before and after training') 
# p1.autoscale(True) 
p1.grid(True) 
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1
"""
"""
im=p1.imshow(nTE,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('nTE after running') 
p1.grid(True) 
#  p1.autoscale(True) 
#}
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1




im=p1.imshow(mbin,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Adjacency matrix both transmitters') 
p1.autoscale(True) 
p1.grid(True) 
# p1.figure(5) 
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1


im=p1.imshow(Matrix, norm=LogNorm(),interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix Excitatory and Inhibitory Connections') 
p1.autoscale(True) 
p1.grid(True) 
# p1 .grid(b='on',which='minor')   ##, which='major', axis='both', **kwargs)	
# p1.ax.grid(color='r',linewidth=1) 

# im=p1.imshow(Matrix)
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1
# im=p1.imshow(MatrixG) 
im=p1.imshow(MatrixG, norm=LogNorm(),interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix GABA') 
p1.grid(True) 
p1.autoscale(True) 

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1
# im=p1.imshow(MatrixA) 
im=p1.imshow(MatrixA, norm=LogNorm(),interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix AMPA') 
p1.grid(True) 
p1.autoscale(True) 


sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1
#import networkx as nx 
#print ma
#Ga=nx.from_numpy_matrix(ma) 
nx.draw_networkx(Ga,label='Structural connectivity, AMPA, degree n') 
p1.title('Structural connectivity, AMPA, degree n') 
p1.draw()

#nx.draw_networkx_labels(Ga,'Structural connectivity, AMPA, degree n') 


sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1
#   import networkx as nx 
Gg=nx.from_numpy_matrix(mg) 
nx.draw_networkx(Gg,label='Structural connectivity, GABA, degree n') 

p1.title('Structural connectivity, GABA, degree n')
p1.draw()
#p1.show()
# nx.draw_networkx(Gg) 
# p1.title('Structural connectivity, GABA, degree n') 
#sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
#fin=fin+1
#   import networkx as nx 
# Gg=nx.from_numpy_matrix(mg) 

#   draw_networkx_labels(Ga,'Structural connectivity, AMPA, degree n') 
#p1.draw()

#nx.draw_networkx_labels(mbin2,label='Structural connectivity,degree 1') 
# # # #/
#Show for every graph.
# #/




# p1.show() 
"""

ps='pruned by= '+str(h.prunenet)+' (ums) '+ 'feed forward (true/false)= '+str(h.ff)
# pos=nx.spring_layout(mbin) # positions for all nodes 
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1

#mbin2=nx.from_numpy_matrix(mbin) 
nx.draw_networkx(mbin2,label='Structural connectivity, degree 1') 
p1.title(ps+'Structural connectivity, degree 1') 
p1.draw() 
p1.savefig('Structural_connectivity_degree 1..png')



sfin=str(h.prunenet)+str(fin)+'.png' 
#sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1

#pos=nx.spring_layout(dir2) # positions for all nodes 
#nx.draw_networkx_nodes(dir2,pos,node_size=700) 
#nx.draw_networkx_edges(dir2,pos,width=2.5) 
#nx.draw_networkx_labels(dir2,pos,font_size=20,font_family='sans-serif')

p1.clf()
fig.clf()
fig=p1.figure(sfin)
nx.draw_networkx(dir2)#,label=ps+'Effective connectivity via pref direction, degree 1') 
p1.title(ps+'Effective connectivity, via pref direction, degree 1') 
p1.draw() 

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig('core'+sfin) 

fig.savefig('Effective_connectivity_degree_1.png')


execfile('gct.py')



# p1.show() 

"""

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1
im=p1.imshow(entfa,interpolation='nearest') 
im=p1.grid(b='True',color='white',which='both', axis='both') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
ff=h.ff
prunenet=h.prunenet
title_str='Between cell nTE across run'+'ff='+str(ff)+'prunenet='+str(prunenet) 
p1.grid(True) 
p1.title(title_str) 
#

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1
im=p1.imshow(entfa, norm=LogNorm(),interpolation='nearest') 
p1.autoscale(True) 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 

prunenet=prunenet
title_str='scaled Between cell nTE across run'+'ff='+str(ff)+'prunenet='+str(prunenet) 
p1.grid(True) 
p1.title(title_str) 

#print idvec.max*idvec.max, cells.count*cells.count, idvec.max*idvec.max
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1
im=p1.imshow(corr,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
ff=ff
prunenet=prunenet
title_str='Between cell correlations for run'+'ff='+str(ff)+'prunenet='+str(prunenet) 
p1.title(title_str) 

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1
im=p1.imshow(dir, norm=LogNorm(),interpolation='nearest') 
# p1.grid(True) 

p1.grid(b='True',color='black',which='both',axis='both',linestyle='--', linewidth=0.4) 
p1.autoscale(True) 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
ff=ff
prunenet=prunenet
title_str='Between cell prefered direction for run'+'ff='+str(ff)+'prunenet='+str(prunenet) 
p1.title(title_str) 

p1.show() 
"""

 
""" efc = nx.DiGraph() 
efc.add_nodes_from(range(0,int(numcell)+2)) 
efc.add_edges_from(entfa)  """
entfa2=nx.to_networkx_graph(entfa,create_using=nx.DiGraph())   #directed graph. 

nx.draw_networkx(entfa2,label='significance based Effective Connectivity') 
#p1.title('Effective Connectivity') 
p1.draw()

# this is a linear filter, I probably want gaussian smoothing, and butterworth bandpass.

"""
theta1 = signal.lfilter(0, 30, float(ecpv1))
gamma1 = signal.lfilter(0.80, 0.100, ecpv)


theta2 = signal.lfilter(0, 30, ecpv2)
gamma2 = signal.lfilter(80, 100, ecpv2)


fin+=1
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin)
p1.hold(True)
p1.plot(theta1)
p1.plot(gamma1)
p1.plot(theta2)
p1.plot(gamma2)
p1.hold(False)
p1.show()  

if(h.ff!=1):
outout=downsample(time_courses[int(outdegree)],oldrate=40000,newrate=200)
outin=downsample(time_courses[int(indegree)],oldrate=40000,newrate=200)
synapse=downsample(recinsyn[int(indegree)],oldrate=40000,newrate=200)
F,pp,cohe,Fx2y,Fy2x,Fxy=granger(outout,outin,20)#,200,100)# nfs2[0], nfs2[1] ,
fin+=1
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin) 
# Plot Granger spectra
p1.hold(True) #Matlab style hold
labels=list()
labels=['Fx2y','Fy2x','Fxy','cohe','pp']
alldata=list()
colors=[[0,0,1],[1,0,0],[0,0.5,0.5],[0.5,1,0],[1,0.5,1]]
lb1='causality between cells i= '+str(int(outdegree)) +' and cell j= '+str(int(indexj))
lb2='causality between cells i= '+str(int(indegree)) +' and cell j= '+str(int(indexi))
p1.plot(F,Fxy,label=lb1,linewidth=1.5,c=colors[0])
p1.plot(F,Fy2x,label=lb2,linewidth=1.5,c=colors[1])
#p1.xlim(0,2)
#p1.ylim(-0.4,0.4)
p1.xlabel("Hz")
p1.xlim(0,100)
p1.ylim(-2,2)
p1.title('GSC components'+bc)
p1.legend()#bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
#      ncol=2, mode="expand", borderaxespad=0.)
p1.hold(False)

F,pp,cohe,Fx2y,Fy2x,Fxy=granger(outout,synapse,20)#,200,100)# nfs2[0], nfs2[1] ,
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin) 
# Plot Granger spectra
p1.hold(True) #Matlab style hold
labels=list()
labels=['Fx2y','Fy2x','Fxy','cohe','pp']
alldata=list()
colors=[[0,0,1],[1,0,0],[0,0.5,0.5],[0.5,1,0],[1,0.5,1]]
lb1='causality between cells i= '+str(int(outdegree)) +' and cell j= '+str(int(indexj))
lb2='causality between cells i= '+str(int(indegree)) +' and cell j= '+str(int(indexi))
p1.plot(F,Fxy,label=lb1,linewidth=1.5,c=colors[0])
p1.plot(F,Fy2x,label=lb2,linewidth=1.5,c=colors[1])
#p1.xlim(0,2)
#p1.ylim(-0.4,0.4)
p1.xlabel("Hz")
p1.xlim(0,100)
p1.ylim(-2,2)
p1.title('GSC components'+bc)
p1.legend()#bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
#      ncol=2, mode="expand", borderaxespad=0.)
p1.hold(False)

#F,pp,cohe,Fx2y,Fy2x,Fxy=granger(outout,synapse,20)#,200,100)# nfs2[0], nfs2[1] ,
p1.show() 
"""	



"""
def getdirm():	  
 dirsgc=nx.DiGraph()#directed graph. 
 dirsgc=np.array(sgcdir)")
 indexi= 0
 indexj=0
 for i in xrange(0,int(h.ntefa.v[0].size-1)): 

   sources=ntefa.v[0].x[iter]
   targets=ntefa.v[1].x[iter]

  ent=ntefa.v[2].x[iter] //get the nte convert to pythonic form
  //print py.ent
  
  if(ntefa.v[2].x[iter]>thr):
      print py.ent
      cnt+=1
      py.indexi= py.sources
      py.indexj= py.targets
      py.fingure_cnt=py.int(cnt)
  

   entfa[int(sources)][int(targets)] = ent

  if(ent>thr/10):      
      dir2.add_edge(int(sources),int(targets),weight=ent)
  


get_dirm()
"""
