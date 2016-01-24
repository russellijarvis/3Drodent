#!/usr/bin/python
# -*- coding: utf-8 -*-

#import py_compile
#py_compile.compile('plotting8.py')

import pyhoc
import brian as brian
import information_theory as infth


# import statsmodels as statsmodels

from matplotlib.colors import LogNorm
import matplotlib

# matplotlib.use('Agg')
# import pylab as p2

import pylab as p1

# from scipy import signal

from bsmart import granger  # Load the Granger calculation tool

# and user make sure that param6.hoc has run.
# sfin=str(h.prunenet)+str(fin)+'.png'
# fig.savefig(sfin,dpi=fig.dpi)

import matplotlib
matplotlib.use('Agg')

input_marke = []
input_markl = []
pspke = []
pspkl = []
rates = []
ratein = []
crlts = []
fhv = []
input_marke = h.input_mark.to_python()
input_markl = h.input_markl.to_python()
pspke = h.pspke.to_python()
pspkl = h.pspkl.to_python()
rates = h.rates.to_python()
ratein = h.ratein.to_python()
across = np.arange(0, numcell, 1)

crlts = h.crlts.to_python()
fhv = h.fhv.to_python()

from mpl_toolkits.axes_grid import make_axes_locatable
import matplotlib.ticker as ticker



""" --------brainx utils------------------------------------------------------
 These utils were copied over from brainx - needed for viz
Copyright (c) 2006-2009, NIPY Developers
All rights reserved.

https://github.com/fperez/brainx/blob/master/LICENSE

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

    * Neither the name of the NIPY Developers nor the names of any
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
     "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
     LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
     A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
     OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
     SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
     LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
     DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
     THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
     (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


def mutual_information(d1, d2):
    """Mutual information between two graph partitions.<-In Brain structure.

Read in two dictionaries of sets (i.e. a graph partition) and assess how
similar they are using mutual information as in Danon, Diaz-Guilera, Duch &
Arenas, J Statistical Mechanics 2005.

Parameters
----------
d1 : dict
dictionary of 'real communities'
d2 : dict
dictionary of 'found communities'

Returns
-------
mi : float
Value of mutual information between the two partitions.
"""
    log = np.log
    nansum = np.nansum

    N = confusion_matrix(d1, d2)

    nsum_row = N.sum(0)[np.newaxis, :]
    nsum_col = N.sum(1)[:, np.newaxis]

    # Sanity checks: a zero in either of these can only happen if there was an
    # empty module in one of the input partitions. Rather than manually check
    # the entire partitions, we look for this problem at this stage, and bail
    # if there was an empty module.
## if (nsum_row==0).any():
## EmptyModuleError("Empty module in second partition.")
## if (nsum_col==0).any():
## EmptyModuleError("Empty module in first partition.")

    # nn is the total number of nodes
    nn = nsum_row.sum()
    num = nansum(N*log(N*nn/(nsum_row*nsum_col)))
    den = nansum(nsum_row*log(nsum_row/nn)) + nansum(nsum_col*log(nsum_col/nn))

    return -2*num/den



def find_unconnected_nodes(self):
    """ checks for nodes in graph with no edges """
    graph = nx.from_numpy_matrix(self.graph_adj_matrix)
    unconnected = [ n for n,d in graph.degree_iter() if d==0 ]
    return unconnected

def drawmatrix_channels(
    in_m,
    channel_names=None,
    fig=None,
    x_tick_rot=0,
    size=None,
    cmap=p1.cm.RdBu_r,
    colorbar=True,
    color_anchor=None,
    title=None,
    ):
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
        fig = p1.figure()

    if size is not None:

        fig.set_figwidth(size[0])
        fig.set_figheight(size[1])

    w = fig.get_figwidth()
    h = fig.get_figheight()

    ax_im = fig.add_subplot(1, 1, 1)

    # If you want to draw the colorbar:

    if colorbar:
        divider = make_axes_locatable(ax_im)
        ax_cb = divider.new_vertical(size='10%', pad=0.1,
                pack_start=True)
        fig.add_axes(ax_cb)

    # Make a copy of the input, so that you don't make changes to the original
    # data provided

    m = in_m.copy()

    # Null the upper triangle, so that you don't get the redundant and the
    # diagonal values:
    # idx_null = triu_indices(m.shape[0])
    # m[idx_null] = np.nan

    # Extract the minimum and maximum values for scaling of the
    # colormap/colorbar:

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

    # The call to imshow produces the matrix plot:

    im = ax_im.imshow(
        m,
        origin='upper',
        interpolation='nearest',
        vmin=color_min,
        vmax=color_max,
        cmap=cmap,
        )

    # Formatting:

    ax = ax_im
    ax.grid(True)

    # Label each of the cells with the row and the column:

    if channel_names is not None:
        for i in xrange(0, m.shape[0]):

             if i < (m.shape[0] - 1):
                ax.text(i - 0.3, 0.3, channel_names[i], rotation=x_tick_rot)
            # I suspect this is the cause of the x axis ticks.
            #This one, increments in x, set in y.
            if i > 0:
            #This one, set in x, incrememnts in y
                ax.text(-1, i + 0.3, channel_names[i],
                        horizontalalignment='right')

        ax.set_axis_off()
        ax.set_xticks(np.arange(N))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(channel_formatter))
        fig.autofmt_xdate(rotation=x_tick_rot)
        ax.set_yticks(np.arange(N))
        ax.set_yticklabels(channel_names)
        ax.set_ybound([-0.50, N - 0.50])
        ax.set_xbound([-0.50, N - 1.5])

    # Make the tick-marks invisible:

    for line in ax.xaxis.get_ticklines():
        line.set_markeredgewidth(0)

    for line in ax.yaxis.get_ticklines():
        line.set_markeredgewidth(0)

    ax.set_axis_off()

    if title is not None:
        ax.set_title(title)

    # The following produces the colorbar and sets the ticks

    if colorbar:

        # Set the ticks - if 0 is in the interval of values, set that, as well
        # as the maximal and minimal values:

        if min_val < 0:
            ticks = [color_min, min_val, 0, max_val, color_max]
        else:

        # Otherwise - only set the minimal and maximal value:

            ticks = [color_min, min_val, max_val, color_max]

        # This makes the colorbar:

        cb = fig.colorbar(
            im,
            cax=ax_cb,
            orientation='horizontal',
            cmap=cmap,
            norm=im.norm,
            boundaries=np.linspace(color_min, color_max, 256),
            ticks=ticks,
            format='%.2f',
            )

    # Set the current figure active axis to be the top-one, which is the one
    # most likely to be operated on by users later on

    fig.sca(ax)

    return fig


def drawmatrix_channels_nullh(
    in_m,
    channel_names=None,
    fig=None,
    x_tick_rot=0,
    size=None,
    cmap=p1.cm.RdBu_r,
    colorbar=True,
    color_anchor=None,
    title=None,
    ):
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
        fig = p1.figure()

    if size is not None:

        fig.set_figwidth(size[0])
        fig.set_figheight(size[1])

    w = fig.get_figwidth()
    h = fig.get_figheight()

    ax_im = fig.add_subplot(1, 1, 1)

    # If you want to draw the colorbar:

    if colorbar:
        divider = make_axes_locatable(ax_im)
        ax_cb = divider.new_vertical(size='10%', pad=0.1,
                pack_start=True)
        fig.add_axes(ax_cb)

    # Make a copy of the input, so that you don't make changes to the original
    # data provided

    m = in_m.copy()

    # Null the upper triangle, so that you don't get the redundant and the
    # diagonal values:
    idx_null = triu_indices(m.shape[0])
    m[idx_null] = np.nan

    # Extract the minimum and maximum values for scaling of the
    # colormap/colorbar:

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

    # The call to imshow produces the matrix plot:

    im = ax_im.imshow(
        m,
        origin='upper',
        interpolation='nearest',
        vmin=color_min,
        vmax=color_max,
        cmap=cmap,
        )

    # Formatting:

    ax = ax_im
    ax.grid(True)

    # Label each of the cells with the row and the column:

    if channel_names is not None:
        for i in xrange(0, m.shape[0]):

            if i < (m.shape[0] - 1):
                ax.text(i - 0.3, -0.3, channel_names[i], 
                        horizontalalignment='right')
               # I suspect this is the cause of the x axis ticks.

            if i > 0:
                ax.text(-1, i + 0.3, channel_names[i],
                        horizontalalignment='right')

        ax.set_axis_off()
        ax.set_xticks(np.arange(N))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(channel_formatter))
        fig.autofmt_xdate(rotation=x_tick_rot)
        ax.set_yticks(np.arange(N))
        ax.set_yticklabels(channel_names)
        ax.set_ybound([-0.50, N - 0.50])
        ax.set_xbound([-0.50, N - 1.5])

    # Make the tick-marks invisible:

    for line in ax.xaxis.get_ticklines():
        line.set_markeredgewidth(0)

    for line in ax.yaxis.get_ticklines():
        line.set_markeredgewidth(0)

    ax.set_axis_off()

    if title is not None:
        ax.set_title(title)

    # The following produces the colorbar and sets the ticks

    if colorbar:

        # Set the ticks - if 0 is in the interval of values, set that, as well
        # as the maximal and minimal values:

        if min_val < 0:
            ticks = [color_min, min_val, 0, max_val, color_max]
        else:

        # Otherwise - only set the minimal and maximal value:

            ticks = [color_min, min_val, max_val, color_max]

        # This makes the colorbar:

        cb = fig.colorbar(
            im,
            cax=ax_cb,
            orientation='horizontal',
            cmap=cmap,
            norm=im.norm,
            boundaries=np.linspace(color_min, color_max, 256),
            ticks=ticks,
            format='%.2f',
            )

    # Set the current figure active axis to be the top-one, which is the one
    # most likely to be operated on by users later on

    fig.sca(ax)

    return fig


# --------brainx utils------------------------------------------------------
# These utils were copied over from brainx - needed for viz
fin = 0
fig=p1.figure(fin)

fig03 = drawmatrix_channels(np.matrix(Matrix[:len(Matrix)-1]), roi_names[:len(Matrix)-1],
                            size=[10., 10.], color_anchor=0)

# ROI names and the matrix need to be the same length.
# I think I am justified using

fig03.savefig(str(h.plastic) + str(h.ff) + str(tstop) + str(ncell)
              + str(h.prunenet) + str(fin) + 'Degree_Matrix.png')



fig03 = drawmatrix_channels_nullh(np.matrix(corr[:len(corr)-1]), roi_names[:len(corr)-1],
                            size=[10., 10.], color_anchor=0)

# ROI names and the matrix need to be the same length.
# I think I am justified using

fig03.savefig(str(h.plastic) + str(h.ff) + str(tstop) + str(ncell)
              + str(h.prunenet) + str(fin) + 'Correlation_Matrix.png')


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

    return idx + (vals, )


def thresholded_arr(
    arr,
    threshold=0.0,
    threshold2=None,
    fill_val=np.nan,
    ):
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
    (idx, vals) = (mth[:-1], mth[-1])
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


# fin is the figure index, for refering to unique figures. Initialise to zero.

outdegree = int(h.outdegree)
indegree = int(h.indegree)
outdegreeg=int(h.outdegreeg)
outdegreeg=int(h.indegreeg)





pag = np.array(intervec)
outig = pag[np.where(pag == outdegreeg)]  # Index of the input cell in pyramidal cells.
tvg = np.array(tvec)
outtg = tvg[np.where(pag == outdegreeg)]

pa = np.array(pyra)
outi = pa[np.where(pa == outdegree)]  # Index of the input cell in pyramidal cells.
tv = np.array(tvec)
outt = tv[np.where(pa == outdegree)]
fig = p1.figure(fin)
fig.clf()
title = 'Raster Plot' + str(int(h.prunenet)) + str(int(h.ff)) \
    + str(int(fin))
p1.title(title)
p1.hold(True)
colors = array([  # Colors for each cell population
    [0.42, 0.67, 0.84],
    [0.50, 0.80, 1.0],
    [0.90, 0.32, 0.0],
    [0.34, 0.67, 0.67],
    [0.42, 0.82, 0.83],
    [0.90, 0.59, 0.0],
    [0.33, 0.67, 0.47],
    [0.42, 0.83, 0.59],
    [0.90, 0.76, 0.0],
    [1.0, 0.85, 0.0],
    [0.71, 0.82, 0.41],
    [0.57, 0.67, 0.33],
    [1.0, 0.38, 0.60],
    ])
j = len(colors) - 1  # 12
p1.plot(tvec, intervec, 'bo', label='inhibitory interneuron')  # ,c=colors[j], markeredgecolor = 'none')
j -= 1
p1.plot(tvec, pyra, 'g^')  # o',c=colors[j], markeredgecolor = 'none')

# p1.plot(tvec[outi],pyra[outi],'yo')#o',c=colors[j], markeredgecolor = 'none')
if h.ff==0:
    p1.plot(outt, outi, 'yo', label='input cell')  # o',c=colors[j], markeredgecolor = 'none')
    p1.plot(outtg, outig, 'yo', label='GABA outdegree')  # o',c=colors[j], markeredgecolor = 'none')

p1.plot(tvec, zerovec, 'g^', label='pyramidal cell')  # ,c=colors[j], markeredgecolor = 'none')
j -= 1
p1.plot(record_inputs, vecin, 'ro', linewidth=15,
        label='synapse input stim')  # , markeredgecolor = 'none')
p1.legend(bbox_to_anchor=(0.0, 1.02, 1.0, .102), loc=9, ncol=2,
          mode='expand', borderaxespad=0.0)

maxtime = int(1000)



"""
spikes = []
num_cells = 10
num_spikes_per_cell = 20
frequency = 20
"""

# Make the spike data. Use a simple Poisson-like spike generator 
# (just for illustrative purposes here. Better spike generators should 
# be used in simulations).
#for i in range(num_cells):
    #isi=numpy.random.poisson(frequency, num_spikes_per_cell)
#    spikes.append(numpy.cumsum(isi))
    
# spikes is now a list of lists where each cell has a list of spike
# times. Now, let's plot these spikes with the default parameters.

import neuronpy.util.spiketrain as st
import numpy
from neuronpy.graphics import spikeplot
st = spikeplot.SpikePlot(sth_ratio=0.3, savefig=True)
st.plot_spikes(tvec)


#st.get_isi_vec(tvec)

# linestyle, 
fin+=1
fig=p1.figure(fin)
p1.hist(st.get_histogram(tvec),ls = 'steps')
sfin='interspikeintervalhistogram'+ + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
fig.savefig(sfin)

fin+=1
fig=p1.figure(fin)
p1.plot(np.hist(st.get_isi_vec(tvec))
sfin='interspikeintervalhistogram2'+ + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
fig.savefig(sfin)
##Above works

def firing_rate(spikes):
    '''
    Rate of the spike train.
    '''
    if spikes==[]:
        return NaN
    return (len(spikes) - 1) / (spikes[-1] - spikes[0])

def CV(spikes):
    '''
    Coefficient of variation.
    '''
    if spikes==[]:
        return NaN
    ISI = diff(spikes) # interspike intervals
    return std(ISI) / mean(ISI)   

def ISI(spikes):
    if spikes==[]:
        return NaN
    return diff(spikes)
   
print ISI(tvec)
#send as ISI, as an argument to entrop

print firing_rate(tvec)
print CV(tvec)
###
# result=ent(prob(a1,len(a1)))
###


#,h.times[j].to_python())
      #brian.correlogram(T1,T2,width=20*ms,bin=1*ms,T=None)
# manrnpxtime=h.max(h.tvec)

p1.xlim(0, tstop)  # maxtime) # To convert to seconds
p1.ylim(-2, int(ncell) + 10)  # Just larger than the number of cells in the model
p1.ylabel('Cell number')
p1.xlabel('spike time (ms)')
p1.hold(False)
sfin = 'raster_plot' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'


##
fin += 1
fig = p1.figure(fin)
fig.clf()

p1.hold(True)
tc = np.array(time_courses[int(1)])
N = len(tc)

yint = np.linspace(1, 1, N)
t = np.linspace(0.0, 0.025 * N, N)
tc = np.array(time_courses[int(indegree)])
str3 = 'cell number= ' + str(10)

p1.plot(t, tc, linewidth=1.5, label=str3)
p1.plot(t, h.spk[int(indegree)].to_python(), linewidth=1.5)
#p1.scatter(yint,h.spk[int(indegree)].to_python(),'g^', linewidth=1.5)
p1.title('pyramidal neuron membrane potential_check')
p1.xlabel('ms')
p1.ylabel('mV')

p1.legend()
p1.hold(False)
sfin = 'pyr_nrn_memb_potential_in_degree_check' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin,dpi=fig.dpi)

fin += 1
fig = p1.figure(fin)
fig.clf()
for i in xrange(0, int(len(time_courses))):
    string = 'voltages'
    bc = ' cells i= ' + str(int(i)) + 'ff= ' + str(int(h.ff)) \
        + 'prune net= ' + str(int(h.prunenet))
    string = string + str(i) + bc
    tc = np.array(time_courses[int(i)])
    p1.plot(t[500:1500], tc[500:1500], linewidth=3)

   # p1.plot(tc[0:25],out1[0:25])

    p1.title('cells index membrane potential')
    p1.xlabel('ms')
    p1.ylabel('mV')

   # out1=downsample(time_courses[int(i)],oldrate=40000,newrate=200)
   # downsampled.append(out1)

sfin = 'traces' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
    + str(int(h.ff)) + str(ncell) + str(int(fin)) + str(int(h.minr)) \
    + '.png'
fig.savefig(sfin,dpi=fig.dpi)
ncell = int(h.ncell)

###

fin += 1
fig = p1.figure(fin)
fig.clf()

###

p1.title('Entropy Plot')
p1.hold(True)
j = len(colors) - 1

outa = across[np.where(across == outdegree)]  # Index of the input cell in pyramidal cells.

# p1.plot(tvec,intervec,'bs')#p1.plot(tvec,pyra,'g^')

p1.plot(across, input_marke, 'r-', linewidth=2,
        label='H(x) of synapse input stim')  # ,c=colors[j], markeredgecolor = 'none')
j -= 1
p1.plot(outa, pspke[int(outa)], 'y^', label='input cell')
p1.ylim(0, max(input_marke))
p1.xlabel('cell number')
p1.ylabel('bits/sec')
p1.hold(False)
p1.legend(bbox_to_anchor=(0.0, -0.106, 1.0, .102), loc=8, ncol=2,
          mode='expand', borderaxespad=0.0)

p1.ylim(0, max(pspke))

sfin = 'Entropy Plot' + str(int(h.prunenet)) + str(int(h.ff)) \
    + str(int(fin))
fig.savefig(sfin,dpi=fig.dpi)

j = len(colors) - 1

across = np.arange(0, int(ncell), 1)

# p1.plot(tvec,intervec,'bs')#p1.plot(tvec,pyra,'g^')

divin = np.zeros(int(numcell))

# need to use np to remove inf

np.divide(input_marke, ratein, divin)
winfs = isinf(divin)  # replace infinities with 40
divin[winfs] = 40
j -= 1
divout = np.zeros(int(ncell))
np.divide(pspke, rates, divout)

##
##

fin += 1
fig = p1.figure(fin)
fig.clf()
p1.hold(True)
p1.plot(across, divin, 'r-', linewidth=2,
        label='H(x) of synapse input stim')  # ,'b|',c=colors[j], markeredgecolor = 'none')

p1.plot(ratein, input_marke, 'go', label='H(x) of cell num spike train')
p1.plot(outa, input_marke[int(outa)], 'y^', label='input cell')

# p1.plot(across,divout,'go',label='H(x) of cell num spike train')

p1.legend(bbox_to_anchor=(0.0, -0.106, 1.0, .102), loc=8, ncol=2,
          mode='expand', borderaxespad=0.0)
p1.title('Entropy divided by rate')

# p1.ylim(0,max(input_marke))

p1.xlabel('cell number')
p1.ylabel('bits/sec')
p1.hold(False)
sfin = 'entropy_divided_by_rate' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
fig.savefig(sfin,dpi=fig.dpi)

###
###

fin += 1
fig = p1.figure(fin)
p1.clf()
p1.hold(True)
p1.title('Lempel Ziv')
j = len(colors) - 1
p1.plot(across, input_markl, 'r-', linewidth=2,
        label='H(x) of synapse input stim')  # ,'|',c=colors[j],linewidth=5)
j -= 1
p1.plot(outa, pspkl[int(outa)], 'y^', label='input cell')
p1.plot(across, pspkl, 'go', label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0.0, -0.106, 1.0, .102), loc=8, ncol=2,
          mode='expand', borderaxespad=0.0)

# p1.ylim(0,max(pspkl))

maxtime = int(h.tstop)
p1.title('Lempel Ziv')
p1.xlabel('Neuron number')
p1.ylabel('Bits/sec')
p1.hold(False)

sfin = 'Lempel_Ziv' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
fig.savefig(sfin,dpi=fig.dpi)

###
###

fin += 1
fig = p1.figure(fin)
p1.clf()
p1.title('Lempel Ziv Divided by rate')
p1.hold(True)
j = len(colors) - 1

divin = np.zeros(int(ncell))

# need to use np to remove inf

np.divide(input_markl, ratein, divin)  # dividing a large num by small number creates

# approaching infinitely large number.

winfs = isinf(divin)  # replace zeros with infinity
divin[winfs] = 40
j -= 1
divout = np.zeros(int(ncell))
np.divide(pspkl, rates, divout)

# xdata=transfer[np.where(targets==indegree)]
# ydata=dist[np.where(targets==indegree)]

p1.plot(across, divin, 'r-', linewidth=2,
        label='H(x) of synapse input stim')  # ,'|',c=colors[j],linewidth=5)
j -= 1
p1.plot(across, divout, 'go', label='H(x) of cell num spike train')
p1.plot(outa, divout[int(outa)], 'y^', label='input cell')

# p1.ylim(0,max(divout))

p1.legend()

maxtime = int(h.tstop)

p1.xlabel('Neuron number')
p1.ylabel('Bits/sec')

p1.hold(False)
sfin = 'Lempel_Ziv_Divided_By_Rate' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
fig.savefig(sfin,dpi=fig.dpi)

###
###

fin += 1
fig = p1.figure(fin)
fig.clf()
p1.title('nTE and Correlations relative to Input Stimulus')

###

j = len(colors) - 1
p1.hold(True)
p1.plot(outa, crlts[int(outa)], 'y^', linewidth=2, label='input cell')
p1.plot(
    across,
    crlts,
    'o',
    c=colors[j],
    markeredgecolor='none',
    label='correlations',
    )
j -= 1
p1.plot(
    across,
    fhv,
    'o',
    c=colors[j],
    markeredgecolor='none',
    label='transfer entropies',
    )
p1.ylim(0, 1.1)

# p1.legend()

p1.legend(bbox_to_anchor=(1.0, 0.50))
sfin = 'nTE_versus_correlations' + str(h.prunenet) + str(fin) + '.png'
p1.hold(False)
fig.savefig(sfin,dpi=fig.dpi)

ncell = ncell
n_samples = int(tstop / h.dt) + 2  # data_rec.shape[0]
nseq = int(ncell)

roi_names = [0 for x in xrange(0, int(ncell))]
for i in xrange(0, int(ncell - 1)):
    roi_names[i] = str(i) + ' ' + str(h.cells.o(i).reponame)
    print roi_names[i]

# Make an empty container for the data

msgcv = [[0 for x in xrange(0, int(int(int(ncell))))] for x in
         xrange(0, int(int(int(ncell))))]
msgcv2 = [[0 for x in xrange(0, int(int(int(ncell))))] for x in
          xrange(0, int(int(int(ncell))))]

# cant even really use a basic example, because if its full of nans it wont evaluate subsequentely.

from pyhoc import downsample

#  for i=0,cells.count-1{//count up
#    for(j=cells.count-1;j>0;j-=1){//count down

for indexi in xrange(0, int(ncell)):
    for indexj in xrange(0, int(ncell)):

    # indexj=0

        if indexi != indexj:

            vec1 = downsample(time_courses[int(indexi)], oldrate=40000,
                              newrate=200)
            vec2 = downsample(time_courses[int(indexj)], oldrate=40000,
                              newrate=200)
            order = 10
            rate = 200
            maxfreq = 0

            x = array([vec1.to_python(), vec2.to_python()])
            npts = len(vec1)
            fs = 200
            n = len(vec1)
            freq = 100
            p = 15
            ntrls = 1

          # F,pp,cohe,Fx2y,Fy2x,Fxy=pwcausalr(x,ntrls,npts,p,fs,fs/2);

            (  # lag = 2 bins of 5ms each.
                F,
                pp,
                cohe,
                Fx2y,
                Fy2x,
                Fxy,
                ) = granger(vec1.to_python(), vec2.to_python(), 1)

          # the order is the number of lags. It could be that using too many lags means relationships are found which don't mean anything in terms of causality. Strong correlations that arise from the fact too much data is being used.

            if ~np.isnan(np.mean(Fx2y)):
                msgcv2[indexi][indexj] = np.mean(Fx2y) - np.mean(Fy2x)
                msgcv[indexi][indexj] = np.mean(Fx2y)
            else:

      # Since granger values are x->y and are not constrained, it's more informative to look at the values for x -> y minus y -> x (to reduce the confound of simultaneous influence) We'll create that difference matrix and plot it as well.

                msgcv[indexi][indexj] = 0

# Below, remove Not A Numbers (NANs), replace nans with 0s

for indexi in xrange(0, int(ncell - 1)):
    for indexj in xrange(0, int(ncell - 1)):
        if ~np.isnan(msgcv[indexi][indexj]):
            msgcv[indexi][indexj] = msgcv[indexi][indexj]
        else:
            msgcv[indexi][indexj] = 0
        if ~np.isnan(msgcv[indexi][indexj]):
            msgcv2[indexi][indexj] = msgcv[indexi][indexj]
        else:
            msgcv2[indexi][indexj] = 0

# Now that I have limited the number of connections these colours for the degree matrices are realistic.

fin = fin + 1
fig = p1.figure(fin)
fig.clf()
im = p1.imshow(Matrix, interpolation='nearest')

# fig=drawmatrix_channels(Matrix, channel_names=roi_names, fig=fin, x_tick_rot=0,
#                        size=None, cmap=p1.cm.RdBu_r, colorbar=True,
#                        color_anchor=None, title='Degree matrix Excitatory and Inhibitory Connections')
# p1.colorbar(im)

p1.xlabel('columns = targets')
p1.ylabel('rows = sources')
p1.title('Degree matrix Excitatory and Inhibitory Connections')
p1.autoscale(True)
p1.grid(True)
sfin = 'Both_transmitters_degree_matrix' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin,dpi=fig.dpi)

# p1 .grid(b='on',which='minor')   ##, which='major', axis='both', **kwargs)....
# p1.ax.grid(color='r',linewidth=1)

# im=p1.imshow(Matrix)

fin = fin + 1
fig = p1.figure(fin)
fig.clf()

# im=p1.imshow(MatrixG)

# fig=drawmatrix_channels(MatrixG, channel_names=roi_names, fig=fin, x_tick_rot=0,
#                        size=None, cmap=p1.cm.RdBu_r, colorbar=True,
#                        color_anchor=None, title='Degree  Inhibitory Connections')

im = p1.imshow(MatrixG, interpolation='nearest')
p1.autoscale(True)
p1.colorbar(im)
p1.xlabel('columns = targets')
p1.ylabel('rows = sources')
p1.title('Degree matrix GABA')
p1.grid(True)

sfin = 'GABA_degree_matrix' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
fig.savefig(sfin,dpi=fig.dpi)

fin = fin + 1
fig = p1.figure(fin)
fig.clf()

# fig=drawmatrix_channels(MatrixA, channel_names=roi_names, fig=fin, x_tick_rot=0,
#                        size=None, cmap=p1.cm.RdBu_r, colorbar=True,
#                        color_anchor=None, title='Degree matrix Excitatory Connections')
# im=p1.imshow(MatrixA)

im = p1.imshow(MatrixA, interpolation='nearest')
p1.colorbar(im)
p1.xlabel('columns = targets')
p1.ylabel('rows = sources')
p1.title('Degree matrix AMPA')
p1.grid(True)
p1.autoscale(True)

sfin = 'AMPA_degree_matrix' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
fig.savefig(sfin,dpi=fig.dpi)

fin + 1
fig = p1.figure(fin)
fig.clf()

# x=np.matrix(msgcv)

# fig=drawmatrix_channels(msgcv, channel_names=roi_names, fig=fin, x_tick_rot=0,
#                        size=None, cmap=p1.cm.RdBu_r, colorbar=True,
#                        color_anchor=None, title='SGC Effective Connectivity')

im = p1.imshow(msgcv, interpolation='nearest')  # norm=LogNorm()
p1.colorbar(im)

# p1.autoscale(True)

p1.grid(True)
p1.xlabel('columns = targets')
p1.ylabel('rows = sources')
p1.title('SGC Effective Connectivity')

sfin = 'SGC_matrix_imshow' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
fig.savefig(sfin,dpi=fig.dpi)

fin = fin + 1
fig = p1.figure(fin)
fig.clf()

im = p1.imshow(entfa, interpolation='nearest')
p1.colorbar(im)

# p1.autoscale(True)  norm=LogNorm()

p1.grid(True)
p1.xlabel('columns = targets')
p1.ylabel('rows = sources')
p1.title('nTE Effective Connectivity')

sfin = 'nTE_matrix_imshow' + str(h.prunenet) + str(fin) + '.png'
fig.savefig(sfin,dpi=fig.dpi)

# im=p1.imshow(entfa,interpolation='nearest')
# ax = gca()

# ax.get_xaxis().set_major_formatter(FuncFormatter(lambda x,p :"%.2f"%(x/44100)))
# ax.get_yaxis().set_major_locator(FixedFormatter(roi_names))
# ax.get_yaxis().set_major_locator(LinearLocator(0.01))
# ax.get_yaxis().set_major_formatter(FixedFormatter(roi_names))
# draw()

fin + 1
fig = p1.figure(fin)
fig.clf()
im = p1.imshow(msgcv2, interpolation='nearest')
p1.colorbar(im)

# p1.autoscale(True)  norm=LogNorm()

p1.grid(True)

# p1.imshow(dir, interpolation='nearest')

sfin = 'SGC_Mean_Subtraction' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
p1.title('SGC_Mean_Subtraction')

fig.savefig(sfin,dpi=fig.dpi)

#
#

fin += 1
fig = p1.figure(fin)
fig.clf()
im = p1.imshow(corr, interpolation='nearest')
p1.colorbar(im)

# p1.autoscale(True)  norm=LogNorm()

p1.grid(True)

# p1.imshow(dir, interpolation='nearest')

sfin = 'correlations' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
p1.title('correlations_Subtraction')
fig.savefig(sfin,dpi=fig.dpi)

# fig03 = drawmatrix_channels(np.matrix(corr), roi_names, size=[10., 10.], color_anchor=0)
# fig03.savefig(str(h.plastic)+str(h.ff)+str(tstop)+str(ncell)+str(h.prunenet)+str(fin)+'correlation.png')

##
# Effort to threshold. SHould use msu.threshold instead.
#

m2 = thresholded_arr(np.array(entfa), 0, .2, 0)  # replace with zeros

# values between 0, and 0.4

m2 = thresholded_arr(m2, -.2, 0, 0)  # replace with zeros values bewteen -0.4 and 0

m1 = thresholded_arr(np.array(msgcv2), 0, .2, 0)  # replace with zeros

# values between 0, and 0.4

m1 = thresholded_arr(m1, -.2, 0, 0)  # replace with zeros values bewteen -0.4 and 0

dirg = nx.to_networkx_graph(m2, create_using=nx.DiGraph())
dirg2 = nx.to_networkx_graph(m1, create_using=nx.DiGraph())

fin += 1
fig = p1.figure(fin)
fig.clf()
im = p1.imshow(m1, interpolation='nearest')
p1.colorbar(im)
p1.grid(True)
p1.title('subnetwork via sgc')
sfin = 'reduced_sgc_graph' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
fig.savefig(sfin,dpi=fig.dpi)

fin += 1
fig = p1.figure(fin)
fig.clf()
im = p1.imshow(m2, interpolation='nearest')
p1.colorbar(im)
p1.grid(True)
p1.title('subnetwork via nte Effective Connectivity')
sfin = 'reduced_sgc_graph' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
fig.savefig(sfin,dpi=fig.dpi)

fin += 1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(dirg)  # ,nodelist=dictc)#, roi_names)
p1.title('Effective Connectivity via SGC')
p1.draw()
sfin = 'sgc_graph' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
    + str(int(h.ff)) + str(ncell) + str(int(fin)) + str(int(h.minr)) \
    + '.png'
fig.savefig(sfin,dpi=fig.dpi)

# fig05 = draw_graph(dirg)
# fig05.savefig('experiment.png')

fin += 1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(dirg2)  # , roi_names)
p1.title('Effective Connectivity via nTE 2')
p1.draw()
sfin = 'nte_graph' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
    + str(int(h.ff)) + str(ncell) + str(int(fin)) + str(int(h.minr)) \
    + '.png'
fig.savefig(sfin,dpi=fig.dpi)


def highest_graph_list(graph_dict):

# Create ordered tuple of centrality data

    list_items = [(b, a) for (a, b) in graph_dict.iteritems()]

 # Sort in descending order

    list_items.sort()

 # cent_items.reverse()

    return tuple(list_items[:5])


entfa2 = nx.to_networkx_graph(np.array(entfa),
                              create_using=nx.DiGraph())  # directed graph.
msgcv4 = nx.to_networkx_graph(np.array(msgcv),
                              create_using=nx.DiGraph())  # directed graph.
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

out_degrees2 = dirg2.out_degree()

# dictionary node:degree

out_values2 = sorted(set(out_degrees2.values()))
out_hist2 = [out_degrees2.values().count(x) for x in out_values2]

in_degrees = dirg.in_degree()

# dictionary node:degree

in_values = sorted(set(in_degrees.values()))
in_hist = [in_degrees.values().count(x) for x in in_values]
fin += 1
p1.figure(fin)
p1.clf()
p1.hold(True)
p1.plot(in_values, in_hist, 'ro-')

# in-degree

p1.plot(out_values, out_hist, 'bv-')
p1.xscale('log')
p1.yscale('log')
p1.hold(False)

# out-degree

p1.legend(['In-degree', 'Out-degree'])
p1.xlabel('Degree Scale free network')
p1.ylabel('Number of Neurons that have this degree value')
p1.title('Degree Distribution of Network')
sfin = 'dirg_degree_distribution.png' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + '.png'
p1.savefig(sfin)


# p1.close()

def trim_nodes(G, d):
    """ returns a copy of G without 
  the nodes with a degree less than d 
   Acknoledgement, from the glowing Python http://glowingpython.blogspot.com.au/
  Plots a Single-Sided Amplitude Spectrum of y(t)
 
  
  """

    Gt = G.copy()

  # dn=Gt.out_degree()
  # dn =

    dn = nx.degree(Gt)
    for n in Gt.nodes():
        print n, dn[n], d
        if dn[n] <= d:

            Gt.remove_node(n)
    return Gt


def trim_nodes_out(G, d):
    """ returns a copy of G without 
  the nodes with a degree less than d """

    Gt = G.copy()

  # dn=Gt.out_degree()
  # dn =

    dn = Gt.out_degree()

  # dn=nx.degree(Gt)

    for n in Gt.nodes():
        print n, dn[n], d
        if dn[n] <= d:

            Gt.remove_node(n)
    return Gt


# drawing the network without
# nodes with degree less than 10

if max(out_values) == 1:
    threshold = 1
else:
    threshold = max(out_values) - 1
if max(out_values2) == 1:
    threshold2 = 1
else:
    threshold2 = max(out_values2) - 1

Gt = trim_nodes(dirg, threshold)
p1.figure(fin)
p1.clf()

nx.draw(Gt, node_size=5, node_color='r', edge_color='b', alpha=.2)

sfin = 'evt_trimmed' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + '.png'
p1.savefig(sfin)

Gt = trim_nodes_out(dirg2, threshold2)
p1.figure(fin)
p1.clf()

nx.draw(Gt, node_size=5, node_color='r', edge_color='b', alpha=.2)

sfin = 'arguably the most influential cell' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + '.png'
p1.savefig(sfin)

# nodes with degree less than 10

Gt = trim_nodes(dirg2, 10)
fin += 1
p1.figure(fin)
p1.clf()

# nx.draw(Gt,node_size=0,node_color='w',edge_color='b',alpha=.2)

nx.draw(Gt, node_size=5, node_color='b', edge_color='b', alpha=.2)
sfin = 'sgc2_trimmedx' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + '.png'
p1.savefig(sfin)

Gt = trim_nodes(dirg, 10)
fin += 1
p1.figure(fin)
p1.clf()
nx.draw_networkx(Gt)  # ,node_size=5,node_color='r',edge_color='b',alpha=.2)
sfin = 'ent_trimmed3x' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + '.png'
p1.savefig(sfin)

# nodes with degree less than 10

Gt = trim_nodes(dirg2, 10)
fin += 1
p1.figure(fin)
p1.clf()
nx.draw_networkx(Gt)  # ,node_size=5,node_color='r',edge_color='b',alpha=.2)
sfin = 'sgc2_trimmed3x' + str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(ncell) + str(int(fin)) \
    + '.png'
p1.savefig(sfin)

fin += 1


def centrality_scatter(
    dict1,
    dict2,
    path='',
    ylab='',
    xlab='',
    title='',
    line=False,
    fin=0,
    ):

  # Create figure and drawing axis

    fig = p1.figure(fin)  # figsize=(7,7))

  # ax1 = fig.add_subplot(111)
  # Create items and extract centralities

    items1 = sorted(dict1.items())
    items2 = sorted(dict2.items())
    xdata = [b for (a, b) in items1]
    ydata = [b for (a, b) in items2]

  # Add each actor to the plot by ID

    for p in xrange(len(items1)):
        p1.text(x=xdata[p], y=ydata[p], s=str(items1[p][0]), color='b')

 #   Basic network analysis - plotting results

    if line:

    # use np to calculate the best fit

        p1.scatter(xdata, ydata)
        p1.hold(True)
        (slope, yint) = p1.polyfit(xdata, ydata, 1)
        xline = p1.xticks()[0]
        yline = map(lambda x: slope * x + yint, xline)
        p1.plot(xline, yline, ls='--', color='b')

    # Set new x- and y-axis limits

        p1.xlim((0.0, max(xdata) + .15 * max(xdata)))
        p1.ylim((0.0, max(ydata) + .15 * max(ydata)))

    # Add labels and save

        p1.title(title)
        p1.xlabel(xlab)
        p1.ylabel(ylab)
        path = 'plot_centrality' + str(int(h.prunenet)) + str(tstop) \
            + str(h.plastic) + str(int(h.ff)) + str(ncell) \
            + str(int(fin)) + '.png'
    p1.savefig(path)
    return 0


dirg_ud = dirg.to_undirected()
bet_cen = nx.betweenness_centrality(dirg_ud)
print bet_cen, 'betweeness centrality'

# Closeness centrality

clo_cen = nx.closeness_centrality(dirg_ud)
print clo_cen, 'closeness centrality'

centrality_scatter(
    clo_cen,
    bet_cen,
    path='scatter_centrality.png',
    ylab='closeness centrality',
    xlab='betweeness centrality',
    title='spectral granger causality',
    line=True,
    fin=fin,
    )


# sent=DiscreteSystem(np.array(trains[j]),(1,2), np.array(trains[k]), (1,2))
# sent.calculate_entropies(method='plugin', calc=['HX'])
# get_entropies()

def sum_diff():
    accumx = 0
    accumy = 0

    for i in xrange(len(entfa) - 1):
        accumx += input_marke[i] - np.mean(entfa[:][i])
        accumy += input_marke[i] - np.mean((entfa[i])[:])
    return accumx


def sum_diff2(neuron):
    accumx = 0
    accumy = 0
    accumx += input_marke[neuron] - np.mean(entfa[:][neuron])
    accumy += input_marke[neuron] - np.mean((entfa[neuron])[:])
    print accumx
    print accumy
    return accumx


def summary():
    hvc = divout[np.where(divout > divin)]  # high variability cells.
    itemindex = np.where(hvc != inf)
    cell_index = np.where(hvc != inf)
    print 'these cells fired with high variability after adjusting for rate', \
        cell_index

    print 'how correlated (related) was the network behavior (of the cells that fired) compared to the input?'
    print sum(crlts[:int(max(h.idvec.to_python())) - 1]) \
        / max(h.idvec.to_python())
    cell_indexc = np.where(crlts < 0.50)
    print 'these cells where uncorrelated with the output', cell_indexc

 # print dirg
 # print dirg2

    print np.sum(m1) / (numcell * numcell), 'sum m1', np.sum(m2) \
        / (numcell * numcell), 'sum m2'
    dicts = nx.degree_centrality(dirg2)
    dicte = nx.degree_centrality(dirg)

    nodelists = highest_graph_list(dicts)
    nodeliste = highest_graph_list(dicte)

    print 'summary:'
    print 'centrality, sgc'
    print nodelists
    print 'centrality, ent'
    print nodeliste

    print 'effective connectivity based in degree outdegree'

 # print dirg.in_degree()
 # print dirg.out_degree()

    oute = dirg.out_degree()

 # print dirg2.in_degree()

    oute2 = dirg2.out_degree()
    effectout2 = highest_graph_list(oute2)
    effectout = highest_graph_list(oute)
    oute3 = dirg.in_degree()

 # print dirg2.in_degree()

    oute4 = dirg2.in_degree()
    effectout3 = highest_graph_list(oute3)
    effectout4 = highest_graph_list(oute4)

    bet_cen = nx.betweenness_centrality(dirg_ud)

 # print bet_cen, 'betweeness centrality'
 # Closeness centrality

    clo_cen = nx.closeness_centrality(dirg_ud)

    close = highest_graph_list(clo_cen)
    between = highest_graph_list(bet_cen)

    print 'before processing'
    print 'is structure different to function?'
    print 'left function, right structure'
    print highest_graph_list(entfa2.in_degree()), indegree, ' indegree'
    print highest_graph_list(entfa2.out_degree()), outdegree, \
        ' outdegree'

    print highest_graph_list(msgcv4.in_degree()), indegree, ' indegree'
    print highest_graph_list(msgcv4.out_degree()), outdegree, \
        ' outdegree'

    print 'after processing'
    print 'is structure different to function?'
    print 'left function, right structure'
    print 'close, between', close, between
    print 'out degree', effectout, effectout2, outdegree, ' outdegree'
    print 'in degree', effectout3, effectout4, indegree, ' indegree'

    print sum_diff2(int(highest_graph_list(msgcv4.out_degree())[0][0])), \
        'diff between ent-nte for cell', \
        int(highest_graph_list(msgcv4.out_degree())[0][0]), \
        'which was the cell for which structure is different to function.'
    print sum_diff2(int(highest_graph_list(msgcv4.in_degree())[0][0]))
    print 'sum entropies', np.sum(entfa) / (numcell * numcell), \
        'sum of granger causalities', np.sum(msgcv) / (numcell
            * numcell), 'sum of sgc matrix', sum_diff()
    function_cell = int(highest_graph_list(msgcv4.out_degree())[0][0])

 # print int(function_cell)
 # print 'difference between impedences, between structure and function'
 # h('print "This is the difference between input impedences at 20Hz Impedence(Functional_Hub)-Impedence(Structural_Hub)",cells.o(py.int(py.function_cell)).inimp20-cells.o(py.int(py.outdegree)).inimp20, "MOhm ", "If this was negative it implies the structual hub was easier to depolarise due to larget input impedence"')

    (sum(entfa[:]) / (numcell * numcell), 'sum of nte matrix')
    sum_diff()
    h.make_times()
    print 'Above is a check for consistancy between  spk[j].sum()==times[j].sum(),'
    print 'ultimately I should make a cells spike train, an object attribute of the Cell class'

    return 0


