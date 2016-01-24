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
#import statsmodels as statsmodels
from matplotlib.colors import LogNorm 
import matplotlib 
#matplotlib.use('Agg') 
import pylab as p2
import pylab as p1
#from scipy import signal 
from bsmart import granger # Load the Granger calculation tool
#and user make sure that param6.hoc has run.
# sfin=str(h.prunenet)+str(fin)+'.png' 
#fig.savefig(sfin) 

import matplotlib
matplotlib.use('Agg') 

def drawmatrix_channels(in_m, channel_names=None, fig=None, x_tick_rot=0,
                        size=None, cmap=plt.cm.RdBu_r, colorbar=True,
                        color_anchor=None, title=None):
    r"""Creates a lower-triangle of the matrix of an nxn set of values. This is
    the typical format to show a symmetrical bivariate quantity (such as
    correlation or coherence between two different ROIs).

    Parameters
    ----------

    in_m: nxn array with values of relationships between two sets of rois or
    channels

    channel_names (optional): list of strings with the labels to be applied to
    the channels in the input. Defaults to '0','1','2', etc.

    fig (optional): a matplotlib figure

    cmap (optional): a matplotlib colormap to be used for displaying the values
    of the connections on the graph

    title (optional): string to title the figure (can be like '$\alpha$')

    color_anchor (optional): determine the mapping from values to colormap
        if None, min and max of colormap correspond to min and max of in_m
        if 0, min and max of colormap correspond to max of abs(in_m)
        if (a,b), min and max of colormap correspond to (a,b)

    Returns
    -------

    fig: a figure object

    """
    N = in_m.shape[0]
    ind = np.arange(N)  # the evenly spaced plot indices

    def channel_formatter(x, pos=None):
        thisind = np.clip(int(x), 0, N - 1)
        return channel_names[thisind]

    if fig is None:
        fig = plt.figure()

    if size is not None:

        fig.set_figwidth(size[0])
        fig.set_figheight(size[1])

    w = fig.get_figwidth()
    h = fig.get_figheight()

    ax_im = fig.add_subplot(1, 1, 1)

    #If you want to draw the colorbar:
    if colorbar:
        divider = make_axes_locatable(ax_im)
        ax_cb = divider.new_vertical(size="10%", pad=0.1, pack_start=True)
        fig.add_axes(ax_cb)

    #Make a copy of the input, so that you don't make changes to the original
    #data provided
    m = in_m.copy()

    #Null the upper triangle, so that you don't get the redundant and the
    #diagonal values:
    #idx_null = triu_indices(m.shape[0])
    #m[idx_null] = np.nan

    #Extract the minimum and maximum values for scaling of the
    #colormap/colorbar:
    max_val = np.nanmax(m)
    min_val = np.nanmin(m)

    if color_anchor is None:
        color_min = min_val
        color_max = max_val
    elif color_anchor == 0:
        bound = max(abs(max_val), abs(min_val))
        color_min = -bound
        color_max = bound
    else:
        color_min = color_anchor[0]
        color_max = color_anchor[1]

    #The call to imshow produces the matrix plot:
    im = ax_im.imshow(m, origin='upper', interpolation='nearest',
       vmin=color_min, vmax=color_max, cmap=cmap)

    #Formatting:
    ax = ax_im
    ax.grid(True)
    #Label each of the cells with the row and the column:
    if channel_names is not None:
        for i in xrange(0, m.shape[0]):
            if i < (m.shape[0] - 1):
                ax.text(i - 0.3, i, channel_names[i], rotation=x_tick_rot)
            if i > 0:
                ax.text(-1, i + 0.3, channel_names[i],
                        horizontalalignment='right')

        ax.set_axis_off()
        ax.set_xticks(np.arange(N))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(channel_formatter))
        fig.autofmt_xdate(rotation=x_tick_rot)
        ax.set_yticks(np.arange(N))
        ax.set_yticklabels(channel_names)
        ax.set_ybound([-0.5, N - 0.5])
        ax.set_xbound([-0.5, N - 1.5])

    #Make the tick-marks invisible:
    for line in ax.xaxis.get_ticklines():
        line.set_markeredgewidth(0)

    for line in ax.yaxis.get_ticklines():
        line.set_markeredgewidth(0)

    ax.set_axis_off()

    if title is not None:
        ax.set_title(title)

    #The following produces the colorbar and sets the ticks
    if colorbar:
        #Set the ticks - if 0 is in the interval of values, set that, as well
        #as the maximal and minimal values:
        if min_val < 0:
            ticks = [color_min, min_val, 0, max_val, color_max]
        #Otherwise - only set the minimal and maximal value:
        else:
            ticks = [color_min, min_val, max_val, color_max]

        #This makes the colorbar:
        cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal',
                          cmap=cmap,
                          norm=im.norm,
                          boundaries=np.linspace(color_min, color_max, 256),
                          ticks=ticks,
                          format='%.2f')

    # Set the current figure active axis to be the top-one, which is the one
    # most likely to be operated on by users later on
    fig.sca(ax)

    return fig
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
p1.ylim(-2,int(ncell)) # Just larger than the number of cells in the model
p1.ylabel("Cell number")
p1.xlabel("spike time (ms)")
p1.hold(False)
sfin='raster_plot'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
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
sfin='pyr_nrn_memb_potential'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
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
sfin='traces'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin) 
ncell=int(h.ncell)
"""
p1.hold(False)
vsum1 = np.array(vsum1)
vsum2 = np.array(vsum2)
#vtotal=vtotal/2 #so its an average again.
in_trace = np.array(in_trace)
out_trace = np.array(in_trace)
in_trace=in_trace/int(ncell)

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



f
fig.savefig(sfin) 

j=len(colors)-1

across=np.arange(0,int(ncell),1)
#p1.plot(tvec,intervec,'bs')#p1.plot(tvec,pyra,'g^')
divin=np.zeros(int(numcell))
#need to use numpy to remove inf
np.divide(input_marke,ratein,divin)
winfs=isinf(divin) # replace infinities with 40
divin[winfs]=40
j-=1
divout=np.zeros(int(ncell))
np.divide(pspke,rates,divout)



##
##

fin+=1
fig=p1.figure(fin)
fig.clf()
p1.hold(True)
p1.plot(across,divin,'r-',linewidth=2,label='H(x) of synapse input stim')#,'b|',c=colors[j], markeredgecolor = 'none')

p1.plot(ratein,input_marke,'go',label='H(x) of cell num spike train')

#p1.plot(across,divout,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)
p1.title("Entropy divided by rate")
p1.xlabel('cell number')
p1.ylabel('bits/sec')
p1.hold(False)
sfin='entropy_divided_by_rate'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
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

sfin='Lempel_Ziv'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin) 

###
###


fin+=1
fig=p1.figure(fin)
p1.clf()
p1.title("Lempel Ziv Divided by rate")
p1.hold(True)
j=len(colors)-1

divin=np.zeros(int(ncell))
#need to use numpy to remove inf
np.divide(input_markl,ratein,divin) #dividing a large num by small number creates
#approaching infinitely large number.
winfs=isinf(divin) # replace zeros with infinity
divin[winfs]=40
j-=1
divout=np.zeros(int(ncell))
np.divide(pspkl,rates,divout)

#xdata=transfer[np.where(targets==indegree)]
#ydata=dist[np.where(targets==indegree)]


p1.plot(across,divin,'r-',linewidth=2,label='H(x) of synapse input stim')#,'|',c=colors[j],linewidth=5)
j-=1
p1.plot(across,divout,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(h.tstop)

p1.xlabel("Neuron number")
p1.ylabel("Bits/sec")


p1.hold(False)
sfin='Lempel_Ziv_Divided_By_Rate'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin) 



###
###

"""
fin+=1
fig=p1.figure(fin)
fig.clf()
p1.hold(True)
# Two subplots, the axes array is 1-d
f, axarr = p1.subplots(2, sharex=True)
axarr[0].plot(crlts,across,'go', markeredgecolor = 'none',label='correlations')
axarr[0].set_title('correlations between input and cell number')
axarr[1].plot(fhv,across,'go', markeredgecolor = 'none',label='transfer entropies')
axarr[1].set_title('nTE between input and cell number')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)
p1.xlabel('cell number')
"""

p1.hold(False)
sfin='nTE_versus_correlations_sub_plot'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
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

#import pyentropy as pe
#from pyentropy import DiscreteSystem
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
  across2=np.arange(0,int(ncell)*2+1,1)
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
 



ncell=ncell    
n_samples = int(tstop/h.dt)+2 #data_rec.shape[0]
nseq=int(ncell)

roi_names=[0 for x in xrange(0,int(ncell))]
for i in xrange(0,int(ncell-1)):
 roi_names[i]=str(i)+' '+str(h.cells.o(i).nametype)
 print roi_names[i]
#Make an empty container for the data
msgcv= [[0 for x in xrange(0,int(int(int(ncell))))] for x in xrange(0,int(int(int(ncell))))]
msgcv2= [[0 for x in xrange(0,int(int(int(ncell))))] for x in xrange(0,int(int(int(ncell))))]

#cant even really use a basic example, because if its full of nans it wont evaluate subsequentely.
from pyhoc import downsample

#  for i=0,cells.count-1{//count up
#    for(j=cells.count-1;j>0;j-=1){//count down 
for indexi in xrange(0,int(ncell)):
  for indexj in xrange(0,int(ncell)):
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
  
for indexi in xrange(0,int(ncell-1)):
  for indexj in xrange(0,int(ncell-1)):
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


""" clutter
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
sfin='Both_transmitters_adjacency_matrix'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin)
"""


fin=fin+1
fig=p1.figure(fin) 
fig.clf()
#im=p1.imshow(Matrix, interpolation='nearest') 

fig=drawmatrix_channels(Matrix, channel_names=roi_names, fig=fin, x_tick_rot=0,
                        size=None, cmap=plt.cm.RdBu_r, colorbar=True,
                        color_anchor=None, title='Degree matrix Excitatory and Inhibitory Connections')
#p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
#p1.title('Degree matrix Excitatory and Inhibitory Connections') 
#p1.autoscale(True) 
#p1.grid(True) 
sfin='Both_transmitters_degree_matrix'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin)
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
sfin='GABA_degree_matrix'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin)

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
sfin='AMPA_degree_matrix'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin)




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
sfin='SGC_matrix_imshow'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
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
sfin='SGC_Mean_Subtraction'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
p1.title('SGC_Mean_Subtraction')

fig.savefig(sfin)
#
#
fin+=1
fig=p1.figure(fin)
fig.clf()
im=p1.imshow(corr,interpolation='nearest')
p1.colorbar(im) 
#p1.autoscale(True)  norm=LogNorm() 
p1.grid(True) 
#p1.imshow(dir, interpolation='nearest')
sfin='correlations'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
p1.title('correlations_Subtraction')
fig.savefig(sfin)






##
#Effort to threshold. SHould use msu.threshold instead.
#

"""


  
fin+=1  
fig01.clf()
fig01 = drawgraph_channels(m2)
sfin='filtered_sgc_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png'   
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
m1 = thresholded_arr(m1,-0.2,0,0) #replace with zeros values bewteen -0.4 and 0


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
fig=p1.figure(fin)
fig.clf()
im=p1.imshow(m1,interpolation='nearest') 
p1.colorbar(im) 
p1.grid(True) 
p1.title('subnetwork via sgc')
sfin='reduced_sgc_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin)  


fin+=1
fig=p1.figure(fin)
fig.clf()
im=p1.imshow(m2,interpolation='nearest') 
p1.colorbar(im) 
p1.grid(True) 
p1.title('subnetwork via nte Effective Connectivity')
sfin='reduced_sgc_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin)
  
"""

fin+=1
fig = p1.figure(fin)
p1.clf()

nx.draw_networkx(dirg2)#,nodelist=dictc)#, roi_names)
p1.title('Effective Connectivity via nte')
p1.draw() 
sfin='sgc_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin)  
"""
 
fin+=1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(dirg)#,nodelist=dictc)#, roi_names)
p1.title('Effective Connectivity via SGC')
p1.draw() 
sfin='sgc_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin)
#fig05 = draw_graph(dirg)
#fig05.savefig('experiment.png')


fin+=1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(dirg2)#, roi_names)
p1.title('Effective Connectivity via nTE 2')
p1.draw() 
sfin='nte_graph'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png' 
fig.savefig(sfin)



def highest_graph_list(graph_dict):
# Create ordered tuple of centrality data
 list_items=[(b,a) for (a,b) in graph_dict.iteritems()]
 # Sort in descending order
 list_items.sort()
 #cent_items.reverse()
 return tuple((list_items[:5]))



entfa2=nx.to_networkx_graph(np.array(entfa),create_using=nx.DiGraph())   #directed graph. 
msgcv4=nx.to_networkx_graph(np.array(msgcv),create_using=nx.DiGraph())   #directed graph. 
print 'is structure different to function?'
print 'left function, right structure'
print highest_graph_list(entfa2.in_degree()), indegree, ' indegree'
print highest_graph_list(entfa2.out_degree()), outdegree, ' outdegree'

print highest_graph_list(msgcv4.in_degree()), indegree, ' indegree'
print highest_graph_list(msgcv4.out_degree()), outdegree, ' outdegree'


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
sfin='dirg_degree_distribution.png'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png'
p1.savefig(sfin)
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
Gt = trim_nodes(dirg,10)
p1.figure(fin)
p1.clf()

nx.draw(Gt,node_size=5,node_color='r',edge_color='b',alpha=.2)

sfin='evt_trimmed'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png'
p1.savefig(sfin)

# nodes with degree less than 10
Gt = trim_nodes(dirg2,10)
fin+=1
p1.figure(fin)
p1.clf()
#nx.draw(Gt,node_size=0,node_color='w',edge_color='b',alpha=.2)
nx.draw(Gt,node_size=5,node_color='b',edge_color='b',alpha=.2)
sfin='sgc2_trimmedx'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png'
p1.savefig(sfin)

Gt = trim_nodes(dirg,10)
fin+=1
p1.figure(fin)
p1.clf()
nx.draw_networkx(Gt)#,node_size=5,node_color='r',edge_color='b',alpha=.2)
sfin='ent_trimmed3x'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png'
p1.savefig(sfin)

# nodes with degree less than 10
Gt = trim_nodes(dirg2,10)
fin+=1
p1.figure(fin)
p1.clf()
nx.draw_networkx(Gt)#,node_size=5,node_color='r',edge_color='b',alpha=.2)
sfin='sgc2_trimmed3x'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png'
p1.savefig(sfin)

fin+=1
def centrality_scatter(dict1,dict2,path="",
  ylab="",xlab="",title="",line=False,fin=0):
  # Create figure and drawing axis

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
    p1.plot(xline,yline,ls='--',color='b')
    # Set new x- and y-axis limits
    p1.xlim((0.0,max(xdata)+(.15*max(xdata))))
    p1.ylim((0.0,max(ydata)+(.15*max(ydata))))
    # Add labels and save
    p1.title(title)
    p1.xlabel(xlab)
    p1.ylabel(ylab)
    path='plot_centrality'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(ncell)+str(int(fin))+'.png'
  p1.savefig(path)
  return 0    
dirg_ud = dirg.to_undirected()    
bet_cen = nx.betweenness_centrality(dirg_ud)
print bet_cen, 'betweeness centrality'
# Closeness centrality
clo_cen = nx.closeness_centrality(dirg_ud)
print clo_cen, 'closeness centrality'

centrality_scatter(clo_cen,bet_cen,path='scatter_centrality.png',ylab="closeness centrality",xlab="betweeness centrality",title="spectral granger causality",line=True,fin=fin)
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
#get_entropies()
def sum_diff():
  accumx=0
  accumy=0


  for i in xrange(len(entfa)-1):
   accumx+=input_marke[i] - np.mean(entfa[:][i])
   accumy+=input_marke[i]- np.mean(entfa[i][:])
  return accumx

def sum_diff2(neuron):
  accumx=0
  accumy=0
  accumx+=input_marke[neuron] - np.mean(entfa[:][neuron])
  accumy+=input_marke[neuron] - np.mean(entfa[neuron][:])
  print accumx
  print accumy
  return accumx
  
def summary():
 hvc=divout[np.where(divout>divin)]#high variability cells.
 itemindex=np.where(hvc!=inf)
 cell_index=numpy.where(hvc!=inf)
 print 'these cells fired with high variability after adjusting for rate', cell_index

 #print dirg 
 #print dirg2
 print np.sum(m1), "sum m1", np.sum(m2), "sum m2"
 dicts=nx.degree_centrality(dirg2)    
 dicte=nx.degree_centrality(dirg)



 nodelists=highest_graph_list(dicts)
 nodeliste=highest_graph_list(dicte)

 print 'summary:'
 print 'centrality, sgc'
 print nodelists
 print 'centrality, ent'
 print nodeliste

 print 'effective connectivity based in degree outdegree'
 #print dirg.in_degree()
 #print dirg.out_degree()
 oute= dirg.out_degree()
 #print dirg2.in_degree()
 oute2= dirg2.out_degree()
 effectout2=highest_graph_list(oute2)
 effectout=highest_graph_list(oute)
 oute3= dirg.in_degree()
 #print dirg2.in_degree()
 oute4= dirg2.in_degree()
 effectout3=highest_graph_list(oute3)
 effectout4=highest_graph_list(oute4)



 # Betweenness centrality
 bet_cen = nx.betweenness_centrality(dirg_ud)
 #print bet_cen, 'betweeness centrality'
 # Closeness centrality
 clo_cen = nx.closeness_centrality(dirg_ud)

 close=highest_graph_list(clo_cen)
 between=highest_graph_list(bet_cen) 

 print 'before processing'
 print 'is structure different to function?'
 print 'left function, right structure'
 print highest_graph_list(entfa2.in_degree()), indegree, ' indegree'
 print highest_graph_list(entfa2.out_degree()), outdegree, ' outdegree'

 print highest_graph_list(msgcv4.in_degree()), indegree, ' indegree'
 print highest_graph_list(msgcv4.out_degree()), outdegree, ' outdegree'

 print 'after processing'
 print 'is structure different to function?'
 print 'left function, right structure'
 print 'close, between', close, between 
 print 'out degree', effectout, effectout2, outdegree, ' outdegree'
 print 'in degree', effectout3, effectout4, indegree, ' indegree'
 # print clo_cen, 'closeness centrality'
 # Eigenvector centrality
 #eig_cen = nx.eigenvector_centrality(dirg_ud)




 return 0
 
 
h.summary()
summary()


print 'summary 3'
print sum_diff2(int(highest_graph_list(msgcv4.out_degree())[0][0])), 'diff between ent-nte for cell', int(highest_graph_list(msgcv4.out_degree())[0][0]), 'which was the cell for which structure is different to function.'
print sum_diff2(int(highest_graph_list(msgcv4.in_degree())[0][0]))
print "sum entropies", np.sum(entfa), "sum of granger causalities", np.sum(msgcv), "sum of sgc matrix", sum_diff()
function_cell=int(highest_graph_list(msgcv4.out_degree())[0][0])
print 'difference between impedences, between structure and function'
h('print "This is the difference between input impedences at 20Hz Impedence(Functional_Hub)-Impedence(Structural_Hub)",cells.o(py.int(py.function_cell)).inimp20-cells.o(py.int(py.outdegree)).inimp20, "MOhm ", "If this was negative it implies the structual hub was easier to depolarise due to larget input impedence"')

sum(entfa), "sum of nte matrix"
sum_diff()


