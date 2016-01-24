"""
os.chdir('nitime')
import nitime
import nitime.analysis as nta
import nitime.timeseries as ts
import nitime.utils as tsu
from nitime.viz import drawmatrix_channels
#from nitime.viz import drawmatrix_channels2
from nitime.viz import drawgraph_channels
from nitime.viz import draw_graph
os.chdir('../')
"""

import pyhoc    
import statsmodels as statsmodels
from matplotlib.colors import LogNorm 
import matplotlib 
#matplotlib.use('Agg') 
import pylab as p2
from scipy import signal 
from bsmart import granger # Load the Granger calculation tool
#and user make sure that param6.hoc has run.
# sfin=str(h.prunenet)+str(fin)+'.png' 
#fig.savefig(sfin) 

import matplotlib
matplotlib.use('Agg') 



# --------brainx utils------------------------------------------------------
# These utils were copied over from brainx - needed for viz


def threshold_arr(cmat, threshold=0.0, threshold2=None):
    """Threshold values from the input array.

    Parameters
    ----------
    cmat : array

    threshold : float, optional.
      First threshold.

    threshold2 : float, optional.
      Second threshold.

    Returns
    -------
    indices, values: a tuple with ndim+1

    Examples
    --------
    >>> np.set_printoptions(precision=4)  # For doctesting
    >>> a = np.linspace(0,0.2,5)
    >>> a
    array([ 0.  ,  0.05,  0.1 ,  0.15,  0.2 ])
    >>> threshold_arr(a,0.1)
    (array([3, 4]), array([ 0.15,  0.2 ]))

    With two thresholds:
    >>> threshold_arr(a,0.1,0.2)
    (array([0, 1]), array([ 0.  ,  0.05]))
    """
    # Select thresholds
    if threshold2 is None:
        th_low = -np.inf
        th_hi = threshold
    else:
        th_low = threshold
        th_hi = threshold2

    # Mask out the values we are actually going to use
    idx = np.where((cmat < th_low) | (cmat > th_hi))
    vals = cmat[idx]

    return idx + (vals,)


def thresholded_arr(arr, threshold=0.0, threshold2=None, fill_val=np.nan):
    """Threshold values from the input matrix and return a new matrix.

    Parameters
    ----------
    arr : array

    threshold : float
      First threshold.

    threshold2 : float, optional.
      Second threshold.

    Returns
    -------
    An array shaped like the input, with the values outside the threshold
    replaced with fill_val.

    Examples
    --------
    """
    a2 = np.empty_like(arr)
    a2.fill(fill_val)
    mth = threshold_arr(arr, threshold, threshold2)
    idx, vals = mth[:-1], mth[-1]
    a2[idx] = vals

    return a2


def rescale_arr(arr, amin, amax):
    """Rescale an array to a new range.

    Return a new array whose range of values is (amin,amax).

    Parameters
    ----------
    arr : array-like

    amin : float
      new minimum value

    amax : float
      new maximum value

    Examples
    --------
    >>> a = np.arange(5)

    >>> rescale_arr(a,3,6)
    array([ 3.  ,  3.75,  4.5 ,  5.25,  6.  ])
    """

    # old bounds
    m = arr.min()
    M = arr.max()
    # scale/offset
    s = float(amax - amin) / (M - m)
    d = amin - s * m

    # Apply clip before returning to cut off possible overflows outside the
    # intended range due to roundoff error, so that we can absolutely guarantee
    # that on output, there are no values > amax or < amin.
    return np.clip(s * arr + d, amin, amax)

#fin is the figure index, for refering to unique figures. Initialise to zero.

fin=0


fig=p1.figure(fin)
fig.clf()
title='Raster Plot'+str(int(h.prunenet))+str(int(h.ff))+str(int(fin))
p1.title(title)
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
p1.hold(False)
sfin='raster_plot'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin)   
fin+=1



fig=p1.figure(fin)
fig.clf()
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
sfin='pyr_nrn_memb_potential'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin) 


fin+=1
fig=p1.figure(fin)
fig.clf()
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
sfin='traces'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin) 
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


fin+=1
fig=p1.figure(fin)
fig.clf()

###
title='Entropy Plot'+str(int(h.prunenet))+str(int(h.ff))+str(int(fin))
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
winfs=isinf(divin) # replace infinities with 40
divin[winfs]=40
j-=1
divout=np.zeros(int(numcell))
np.divide(pspke,rates,divout)
p1.hold(True)
sfin='entropy_plot'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin) 

##
##

fin+=1
fig=p1.figure(fin)
fig.clf()
p1.plot(across,divin,'r-',linewidth=2,label='H(x) of synapse input stim')#,'b|',c=colors[j], markeredgecolor = 'none')

p1.plot(ratein,input_marke,'go',label='H(x) of cell num spike train')

#p1.plot(across,divout,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)
p1.title("Entropy divided by rate")
p1.xlabel('cell number')
p1.ylabel('bits/sec')
p1.hold(False)
sfin='entropy_divided_by_rate'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin) 

###
###

fin+=1
fig=p1.figure(fin)
p1.clf()
p1.hold(True)
p1.title("Lempel Ziv")
j=len(colors)-1
p1.plot(across,input_markl,'r-',linewidth=2,label='H(x) of synapse input stim')#,'|',c=colors[j],linewidth=5)
j-=1
p1.plot(across,pspkl,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(h.tstop)
p1.title("Lempel Ziv")
p1.xlabel("Neuron number")
p1.ylabel("Bits/sec")
p1.hold(False)

sfin='Lempel_Ziv'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin) 

###
###


fin+=1
fig=p1.figure(fin)
p1.clf()
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
p1.plot(across,divout,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(h.tstop)

p1.xlabel("Neuron number")
p1.ylabel("Bits/sec")


p1.hold(False)
sfin='Lempel_Ziv_Divided_By_Rate'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 

###
###

fin+=1
fig=p1.figure(fin)
fig.clf()
p1.hold(True)
# Two subplots, the axes array is 1-d
f, axarr = p1.subplots(2, sharex=True)
axarr[0].plot(across,crlts,'go', markeredgecolor = 'none',label='correlations')
axarr[0].set_title('correlations between input and cell number')
axarr[1].plot(across,fhv,'go', markeredgecolor = 'none',label='transfer entropies')
axarr[1].set_title('nTE between input and cell number')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)
p1.xlabel('cell number')


p1.hold(False)
sfin='nTE_versus_correlations_sub_plot'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)

maxtime=int(h.tstop)
#p1.xlabel("Neuron number")
#p1.ylabel("nTE")
#p1.show()
#axarr[1].xlabel('neuron number')
""" 
"""

fin+=1
fig=p1.figure(fin)
fig.clf()
p1.title("nTE and Correlations relative to Input Stimulus")
###
j=len(colors)-1
p1.hold(True)
p1.plot(across,crlts,'o',c=colors[j], markeredgecolor = 'none',label='correlations')
j-=1
p1.plot(across,fhv,'o',c=colors[j], markeredgecolor = 'none',label='transfer entropies')
p1.legend()
sfin='nTE_versus_correlations'+str(h.prunenet)+str(fin)+'.png' 
p1.hold(False)
fig.savefig(sfin)


import pyentropy as pe
from pyentropy import DiscreteSystem
def get_entropies():
  econt=[0 for x in xrange(0,int(num_cells))]
  mcont=[[0 for x in xrange(0,int(num_cells))] for x in xrange(0,int(num_cells))]


  for k in xrange(0,len(trains)):
    for j in xrange(0,len(trains)):
      if(k!=j):
        trains[j]=[int(i) for i in trains[j]]
        trains[k]=[int(i) for i in trains[k]]
        if((sum(trains[j])!=0)|(sum(trains[k])!=0)):
          #if((sum(trains[j]))!=0)||(sum(trains[k])!=0)):
          #if(sum(trains[k])!=0):
          sent=DiscreteSystem(np.array(trains[j]),(1,2), np.array(trains[k]), (1,2))
          sent.calculate_entropies(method='plugin', calc=['HX'])
          econt.append(sent.H)  #append only for every new k.   
          print sent.H," ",pspke[j]," ",pspke[k]," ",j, k 
   
          #sent.calculate_entropies(method='plugin', calc=['HX', 'HXY'])
          #print sent.I(),j, k
      #mcont[k][j]=sent.I()      

    #This contribution now contains both H(X) and H(X|Y)
    #Need to do this in its own loop if wish to seperate from H(X) from H(X|Y)
  #"""
  #print input_marke
  across2=np.arange(0,int(numcell)*2+1,1)
  print pspke
  print econt
  print pspkl 
  print input_marke[0], econt[0]  


  return econt
  
def get_MI():
  for k in xrange(0,len(trains)):
    for j in xrange(0,len(trains)):
      if(k!=j):
        trains[j]=[int(i) for i in trains[j]]
        trains[k]=[int(i) for i in trains[k]]
        if((sum(trains[j])!=0)|(sum(trains[k])!=0)):
          #if((sum(trains[j]))!=0)||(sum(trains[k])!=0)):
          #if(sum(trains[k])!=0):
          sent=DiscreteSystem(np.array(trains[j]),(1,2), np.array(trains[k]), (1,2))
          #sent.calculate_entropies(method='plugin', calc=['HX'])
          #print sent.H, j, k 
          sent.calculate_entropies(method='plugin', calc=['HX', 'HXY'])
          #print sent.I(),j, k
      mcont[k][j]=sent.I()
  return mcont      
 



numcell=ncell    
n_samples = int(tstop/h.dt)+2 #data_rec.shape[0]
nseq=int(numcell)

roi_names=[0 for x in xrange(0,int(numcell))]
for i in xrange(0,int(numcell-1)):
 roi_names[i]=str(i)+' '+str(h.cells.o(i).nametype)
 print roi_names[i]
#Make an empty container for the data
msgcv= [[0 for x in xrange(0,int(int(int(numcell))))] for x in xrange(0,int(int(int(numcell))))]
msgcv2= [[0 for x in xrange(0,int(int(int(numcell))))] for x in xrange(0,int(int(int(numcell))))]

#cant even really use a basic example, because if its full of nans it wont evaluate subsequentely.
from pyhoc import downsample

#  for i=0,cells.count-1{//count up
#    for(j=cells.count-1;j>0;j-=1){//count down 
for indexi in xrange(0,int(numcell)):
  for indexj in xrange(0,int(numcell)):
    #indexj=0
    if(indexi!=indexj):    


          vec1=downsample(time_courses[int(indexi)],oldrate=40000,newrate=200)
          vec2=downsample(time_courses[int(indexj)],oldrate=40000,newrate=200)
          order=10
          rate=200
          maxfreq=0
          
       
          x=array([vec1.to_python(),vec2.to_python()])
          npts=len(vec1)
          fs=200
          n=len(vec1)
          freq=100
          p=15
          ntrls=1
          
          
          #F,pp,cohe,Fx2y,Fy2x,Fxy=pwcausalr(x,ntrls,npts,p,fs,fs/2);                      
          F,pp,cohe,Fx2y,Fy2x,Fxy=granger(vec1.to_python(),vec2.to_python(),1) #lag = 2 bins of 5ms each.
          #the order is the number of lags. It could be that using too many lags means relationships are found which don't mean anything in terms of causality. Strong correlations that arise from the fact too much data is being used.
                
          if(~np.isnan(np.mean(Fx2y))):
            msgcv2[indexi][indexj]=(np.mean(Fx2y)-np.mean(Fy2x))
            msgcv[indexi][indexj]=np.mean(Fx2y)
       
      #Since granger values are x->y and are not constrained, it's more informative to look at the values for x -> y minus y -> x (to reduce the confound of simultaneous influence) We'll create that difference matrix and plot it as well.
          else:
            msgcv[indexi][indexj]=0

          
#Below, remove Not A Numbers (NANs), replace nans with 0s
  
for indexi in xrange(0,int(numcell-1)):
  for indexj in xrange(0,int(numcell-1)):
     if(~np.isnan(msgcv[indexi][indexj])):
       msgcv[indexi][indexj]= msgcv[indexi][indexj]
     else:
       msgcv[indexi][indexj]=0
     if(~np.isnan(msgcv[indexi][indexj])):
       msgcv2[indexi][indexj]= msgcv[indexi][indexj]
     else:
       msgcv2[indexi][indexj]=0




#Now that I have limited the number of connections these colours for the degree matrices are realistic.

fin=fin+1
fig=p1.figure(fin) 
fig.clf()
im=p1.imshow(mbin,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Adjacency matrix both transmitters') 
p1.autoscale(True) 
p1.grid(True) 
# p1.figure(5) 
fig=p1.figure(fin) 
fin=fin+1
fig.savefig('Both_transmitters_adjacency_matrix.png')

fin=fin+1
fig=p1.figure(fin) 
fig.clf()
im=p1.imshow(Matrix, interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix Excitatory and Inhibitory Connections') 
p1.autoscale(True) 
p1.grid(True) 
fig.savefig('Both_transmitters_degree_matrix.png')
# p1 .grid(b='on',which='minor')   ##, which='major', axis='both', **kwargs)	
# p1.ax.grid(color='r',linewidth=1) 

# im=p1.imshow(Matrix)
fin=fin+1
fig=p1.figure(fin) 
fig.clf()
# im=p1.imshow(MatrixG) 
im=p1.imshow(MatrixG, interpolation='nearest') 
p1.autoscale(True) 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix GABA') 
p1.grid(True) 

fig.savefig('GABA_degree_matrix.png')

fin=fin+1
fig=p1.figure(fin) 
fig.clf()
# im=p1.imshow(MatrixA) 
im=p1.imshow(MatrixA, interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix AMPA') 
p1.grid(True) 
p1.autoscale(True) 
fig.savefig('AMPA_degree_matrix.png')




fin+1
fig=p1.figure(fin)
fig.clf()
#x=np.matrix(msgcv)
im=p1.imshow(msgcv,interpolation='nearest') # norm=LogNorm()
p1.colorbar(im) 
#p1.autoscale(True) 
p1.grid(True) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('SGC Effective Connectivity')
sfin='SGC_matrix_imshow'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin)



fin=fin+1
fig=p1.figure(fin)
fig.clf() 
im=p1.imshow(entfa,interpolation='nearest') 
p1.colorbar(im) 
#p1.autoscale(True)  norm=LogNorm() 
p1.grid(True) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('nTE Effective Connectivity') 
sfin='nTE_matrix_imshow'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)    
   


fin+1
fig=p1.figure(fin)
fig.clf()
im=p1.imshow(msgcv2,interpolation='nearest') 
p1.colorbar(im) 
#p1.autoscale(True)  norm=LogNorm() 
p1.grid(True) 
#p1.imshow(dir, interpolation='nearest')
sfin='SGC_Mean_Subtraction'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
p1.title('SGC_Mean_Subtraction')

fig.savefig(sfin)









##
#Effort to threshold. SHould use msu.threshold instead.
#

"""


  
fin+=1  
fig01.clf()
fig01 = drawgraph_channels(m2)
sfin='filtered_sgc_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png'   
fig01.savefig(sfin)

dirg =nx.DiGraph() 
dirg2 =nx.DiGraph() 
#dirg.add_nodes_from(msgcv) 
for i in xrange(0,len(msgcv)):
 for j in xrange(0,len(msgcv)):
   if (int(h.prunenet>0)):
     if msgcv[i][j]>(np.mean(msgcv)+2.5*np.std(msgcv)):
     
     
     
     
       dirg.add_edge(i,j,weight=msgcv[i][j]) 
     if msgcv2[i][j]>(np.mean(msgcv)+2.5*np.std(msgcv)):
       dirg2.add_edge(i,j,weight=msgcv2[i][j]) 
   else:
     if msgcv[i][j]>(np.mean(msgcv)+1.5*np.std(msgcv)):
       dirg.add_edge(i,j,weight=msgcv[i][j]) 
     if msgcv2[i][j]>(np.mean(msgcv)+1.5*np.std(msgcv)):
       dirg2.add_edge(i,j,weight=msgcv2[i][j]) 
"""   
    

m2 = thresholded_arr(np.array(entfa),0,0.2,0) #replace with zeros 
#values between 0, and 0.4
m2 = thresholded_arr(m2,-0.2,0,0) #replace with zeros values bewteen -0.4 and 0

m1 = thresholded_arr(np.array(msgcv2),0,0.2,0) #replace with zeros 
#values between 0, and 0.4
m1 = thresholded_arr(m2,-0.2,0,0) #replace with zeros values bewteen -0.4 and 0


"""
  # Create the actual plot where the graph will be displayed
    figsize = np.array(figsize, float)
    figsize *= stretch_factor

    fig = p1.figure(figsize=figsize)
    ax_graph = fig.add_subplot(1, 1, 1)
    fig.sca(ax_graph)

    if layout is None:
"""

dirg=nx.to_networkx_graph(m2,create_using=nx.DiGraph())
dirg2=nx.to_networkx_graph(m1,create_using=nx.DiGraph())

  


fin+=1
fig = p1.figure(fin)
p1.clf()

nx.draw_networkx(dir2)#,nodelist=dictc)#, roi_names)
p1.title('Effective Connectivity via dir2')
p1.draw() 
sfin='sgc_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin)  
  
fin+=1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(dirg)#,nodelist=dictc)#, roi_names)
p1.title('Effective Connectivity via SGC')
p1.draw() 
sfin='sgc_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin)
#fig05 = draw_graph(dirg)
#fig05.savefig('experiment.png')

fin+=1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(dirg2)#, roi_names)
p1.title('Effective Connectivity via nTE')
p1.draw() 
sfin='nte_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin)



def highest_centrality(cent_dict):
# Create ordered tuple of centrality data
 cent_items=[(b,a) for (a,b) in cent_dict.iteritems()]
 # Sort in descending order
 cent_items.sort()
 #cent_items.reverse()
 return tuple((cent_items[:5]))


out_degrees = dirg.out_degree()
# dictionary node:degree
out_values = sorted(set(out_degrees.values()))
out_hist = [out_degrees.values().count(x) for x in out_values]

in_degrees = dirg.in_degree()
# dictionary node:degree
in_values = sorted(set(in_degrees.values()))
in_hist = [in_degrees.values().count(x) for x in in_values]
fin+=1
p1.figure(fin)
p1.clf()
p1.hold(True)
p1.plot(in_values,in_hist,'ro-')
# in-degree
p1.plot(out_values,out_hist,'bv-')
p1.hold(False)
# out-degree
p1.legend(['In-degree','Out-degree'])
p1.xlabel('Degree')
p1.ylabel('Number of Neurons that have this degree value')
p1.title('Degree Distribution of Network')
p1.savefig('dirg_degree_distribution.png')
#p1.close()

def trim_nodes(G,d):
  """ returns a copy of G without 
  the nodes with a degree less than d """
  Gt = G.copy()
  #dn=Gt.out_degree()
  #dn = 
  dn=nx.degree(Gt)
  for n in Gt.nodes():
    print n, dn[n], d
    if dn[n] <= d:

      Gt.remove_node(n)
  return Gt

# drawing the network without
# nodes with degree less than 10
Gt = trim_nodes(dirg,4)
fin+=1
p1=figure(fin)
p1.clf()
nx.draw(Gt,node_size=5,node_color='r',edge_color='b',alpha=.2)
p1.savefig('dirg_trimmed.png')

# nodes with degree less than 10
Gt = trim_nodes(dirg2,4)
fin+=1
p1=figure(fin)
p1.clf()
nx.draw(Gt,node_size=5,node_color='r',edge_color='b',alpha=.2)
p1.savefig('dirg2_trimmed.png')


def centrality_scatter(dict1,dict2,path="",
  ylab="",xlab="",title="",line=False):
  # Create figure and drawing axis
  fin+=1
  fig = p1.figure(fin)#figsize=(7,7))
  #ax1 = fig.add_subplot(111)
  # Create items and extract centralities
  items1 = sorted(dict1.items())
  items2 = sorted(dict2.items())
  xdata=[b for a,b in items1]
  ydata=[b for a,b in items2]
  # Add each actor to the plot by ID
  for p in xrange(len(items1)):
    p1.text(x=xdata[p], y=ydata[p],s=str(items1[p][0]), color="b")
 #   Basic network analysis - plotting results
  if line:
    # use NumPy to calculate the best fit
    p1.scatter(xdata,ydata)
    p1.hold(True)
    slope, yint = p1.polyfit(xdata,ydata,1)
    xline = p1.xticks()[0]
    yline = map(lambda x: slope*x+yint,xline)
    ax1.plot(xline,yline,ls='--',color='b')
    # Set new x- and y-axis limits
    p1.xlim((0.0,max(xdata)+(.15*max(xdata))))
    p1.ylim((0.0,max(ydata)+(.15*max(ydata))))
    # Add labels and save
    ax1.set_title(title)
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
  p1.savefig(path)
  return 0    
dirg_ud = dirg.to_undirected()    
bet_cen = nx.betweenness_centrality(dirg_ud)
print bet_cen, 'betweeness centrality'
# Closeness centrality
clo_cen = nx.closeness_centrality(dirg_ud)
print clo_cen, 'closeness centrality'

centrality_scatter(clo_cen,bet_cen,path='scatter_centrality.png',ylab="closeness centrality",xlab="betweeness centrality",title="spectral granger causality",line=True)
"""
def cdc(graph):
g=graph
dc=nx.degree_centrality(g)
nx.set_node_attributes(g,'degree_cent',dc)
dcsorted=sorted(dc.items(),key=itemgetter(1),reverse=True)
#for key,value in dcsorted[0:10]:
#print "highest degree", key, value
return graph, dcsorted
"""


#sent=DiscreteSystem(np.array(trains[j]),(1,2), np.array(trains[k]), (1,2))
#sent.calculate_entropies(method='plugin', calc=['HX'])
get_entropies()

def summary():

 dicts=nx.degree_centrality(dirg2)    
 dicte=nx.degree_centrality(dirg)

 

 nodelists=highest_centrality(dicts)
 nodeliste=highest_centrality(dicte)

 

 print "what are the key values to calculate entropy for?"

 #dirg_components =nx.connected_component_subgraphs(dirg)
 #dirg_mc = dirg_components[0]
 #print dirg_mc, 'connected components'
 # Betweenness centrality
 bet_cen = nx.betweenness_centrality(dirg_ud)
 print bet_cen, 'betweeness centrality'
 # Closeness centrality
 clo_cen = nx.closeness_centrality(dirg_ud)
 # Eigenvector centrality
 eig_cen = nx.eigenvector_centrality(dirg_ud)


 print 'summary:'
 print 'centrality, sgc'
 print nodelists
 print 'centrality, ent'
 print nodeliste
 return 0
 
 
h.summary()
summary()

"""

dicts=nx.degree_centrality(dirg2)    
dicte=nx.degree_centrality(dirg)


nodelists=highest_centrality(dicts)
nodeliste=highest_centrality(dicte)


layout = nx.circular_layout
# Compute positions for all nodes - nx has several algorithms
#if callable(layout):
#nodes = np.array(sorted(nodeliste))
nodes = np.array(sorted(dirg.nodes()))

nnod = len(nodes)
pos = layout(dirg)
edgelist = dirg.edges()
edge_pos = np.asarray([(pos[e[0]], pos[e[1]]) for e in edgelist])

for nod in nodes:
   nx.draw_networkx_nodes(dirg,
   pos,
   nodelist=[nod],
   node_color=None, #ode_colors[nod],
   node_shape=None, #node_shapes[nod],
   node_size=None) #node_sizes[nod])

nx.draw_networkx_edges(dirg,
pos,
edgelist=[(u, v)],
width=min_width + alpha * max_width,
edge_color=edge_color,
style=edge_style)


fig01 = drawmatrix_channels(m1, np.array(roi_names), size= [10., 10.], color_anchor=0)
sfin='granger_channels_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig01.savefig(sfin)
fig01.clf()
#fin=fin+1  
fig01 = drawmatrix_channels(m2, np.array(roi_names), size= [10., 10.], color_anchor=0)
sfin='nte_channels_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig01.savefig(sfin)

#fig2=drawgraph_channels(np.shape(msgcv2), channel_names=None,cmap=p1.cm.RdBu_r, node_shapes=None, node_colors=None, title=None, layout=None, threshold=0.4)
#p1.show()
#mx=msgcv2
#low_values_indices = np.mean(msgcv)> 0#///  # Where values are low
#mx[low_values_indices] = 0#
"""



"""
fin+1
fig=p1.figure(fin)
fig.clf()
im=p1.imshow(entfa,interpolation='nearest') 
p1.colorbar(im) 
p1.grid(True) 
p1.title('nTE Effective Connectivity')
fig.savefig(sfin)
"""

#challenge Matrix has to be square. Its possible that I could figure out the code to only get the nodes with the highest clustering coefficient. But I would still need a sqaure matrix. This would be too hard and take too long.
#fig04 = drawgraph_channels(np.array(entfa))


#idx = np.hstack(entfa[index])
#idx1 = np.vstack([[idx[i]] * 4 for i in range(4)]).ravel()
#idx2 = np.hstack(4 * [idx])
#fig04 = drawgraph_channels(entfa[index], roi_names[index])


#interactive(True)
#p1.plot((xf,yf1))

#nx.draw(dirg2[dictc[0,10]])#, roi_names)
   


"""


fig01 = drawmatrix_channels(coh, roi_names, size=[10., 10.], color_anchor=0)
"""

#fig01.savefig('sgc_matrix..png')

#mbin2=np.array(mbin) 
#mbin2=nx.to_networkx_graph(mbin2,create_using=nx.DiGraph())   #directed graph. 




"""		 
im=p1.imshow(WeightsAb,interpolation='nearest') 

p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Covariance of weight before training') 
# p1.autoscale(True) 
p1.grid(True)  """

#sfin=str(h.prunenet)+str(fin)+'.png' 
#fig.savefig(sfin) 
#fig=p1.figure(fin) 
#fin=fin+1

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
fig=p1.figure(fin) 
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
fig=p1.figure(fin) 
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
fig=p1.figure(fin) 
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
fig=p1.figure(fin) 
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
fig=p1.figure(fin) 
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
fig=p1.figure(fin) 
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
fig=p1.figure(fin) 
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
fig=p1.figure(fin) 
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


ps='pruned by= '+str(h.prunenet)+str(h.ff)+' (ums) '+ 'feed forward (true/false)= '+str(h.ff)
# pos=nx.spring_layout(mbin) # positions for all nodes 
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin) 
fin=fin+1

#mbin2=nx.from_numpy_matrix(mbin) 
nx.draw_networkx(mbin2,label='Structural connectivity, degree 1') 
p1.title(ps+'Structural connectivity, degree 1') 
p1.draw() 
p1.savefig('Structural_connectivity_degree 1..png')



sfin=str(h.prunenet)+str(fin)+'.png' 
#sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin) 
fin=fin+1


#execfile('gct.py')

"""

# p1.show() 

"""

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin) 
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
fig=p1.figure(fin) 
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
fig=p1.figure(fin) 
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
fig=p1.figure(fin) 
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
#entfa2=nx.to_networkx_graph(entfa,create_using=nx.DiGraph())   #directed graph. 

#nx.draw_networkx(entfa2,label='significance based Effective Connectivity') 
#p1.title('Effective Connectivity') 
#p1.draw()

# this is a linear filter, I probably want gaussian smoothing, and butterworth bandpass.

"""
theta1 = signal.lfilter(0, 30, float(ecpv1))
gamma1 = signal.lfilter(0.80, 0.100, ecpv)


theta2 = signal.lfilter(0, 30, ecpv2)
gamma2 = signal.lfilter(80, 100, ecpv2)


fin+=1
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin)
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


#low_values_indices = np.mean(msgcv)+np.std(msgcv) < 0#///  # Where values are low
#msgcv[low_values_indices] = 0#//  # All low values set to 0
 
#fig02 = drawgraph_channels(np.matrix(msgcv), np.array(roi_names), size= [10., 10.], color_anchor=0)
#fig02.savefig('chanel_type_graph.png')

#fig03 = draw_graph(nx.to_networkx_graph(np.matrix(msgcv),create_using=nx.DiGraph()))
#sfin='grange_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
#fig02.savefig(sfin)
#sfin='grange_graphfinal'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 

#fig03.savefig(sfin)
fig=p1.figure(fin)
fig.clf()
p1.imshow(msgcv, interpolation='nearest')
sfin='granger_matrix_imshow'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin)

  
#msgcv=np.array(msgcv) 
#low_values_indices = np.mean(msgcv)+np.std(msgcv) < 0#///  # Where values are low
#msgcv[low_values_indices] = 0#//  # All low values set to 0


fin+=1
fig = p1.figure(fin)
nx.draw(dir2)#, roi_names)
p1.title('add edge if threshold sgc')
p1.draw() 

#fig.savefig(sfin) 
fig=p1.figure(fin)
fig.savefig('core Effective_connectivity_degree_1'+sfin)
#s = g.subgraph(nx.shortest_path(g.to_undirected(),55))
#dir2.subgraph(
#core=nx.strongly_connected_component_subgraphs(dir2)

#for h in nx.connected_component_subgraphs(dir2.to_undirected()):
#   nx.draw(h,pos,node_color='red')
 
#if target_node in h:
# else:
#  nx.draw(h,pos,node_color='white')

#p1.show()
#pos=nx.circular_layout(core) # positions for all nodes 
#nx.draw_networkx_nodes(core,pos,node_size=700) 
#nx.draw_networkx_edges(core,pos,width=2.5) 
#nx.draw_networkx_labels(core,pos,font_size=20,font_family='sans-serif')

#nx.draw_networkx(dir2,nodes=nx.strongly_connected_component_subgraphs(dir2),label=sfin+'Effective connectivity via pref direction, degree 1') 
#sfin=str(h.prunenet)+str(fin)+'.png' 
#p1.title(sfin+'Effective connectivity, via pref direction, degree 1') 


dictc=nx.degree_centrality(dirg2)
#nx.dr aw(dirg2[dictc[0,10]])#, roi_names)
    
fin+=1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(dirg)#,nodelist=dictc)#, roi_names)
sfin='nowsgc'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
p1.title('add edge if threshold sgc')
fig.savefig(sfin)
#fig05 = draw_graph(dirg)
#fig05.savefig('experiment.png')



# This method recentely deleted.
#fin+=1
#fig = p1.figure(fin)
#nx.draw(dirg2)#, roi_names)
#p1.title('add edge if threshold sgc diff x-y')
#p1.draw() 
#fig.savefig(sfin)


fig = drawgraph_channels(np.array(entfa))
fig.clf()
p1.title('Using drawgraph_channels, to be skillfully employed later to Make a Subset Network.')
sfin='nowsgc'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin)


fin+=1
fig = p1.figure(fin)
fig.clf()
nx.draw(dir2)#, roi_names)
p1.title('add edge threshold ent ')
#from post_analysis4.hoc
p1.draw() 
sfin='nowsgc'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin)   






"""

"""

fin+1
fig=p1.figure(fin)
fig.clf()
#Dir matrix not constructed properly.
im=p1.imshow(dir2,interpolation='nearest') 
p1.colorbar(im) 
#p1.autoscale(True)   norm=LogNorm()
p1.grid(True) 
sfin='entfa_pref_dir_imshow'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(int(fin))+'.png' 
fig.savefig(sfin)   
"""

