#!/bin/env python
#Unoriginal code concepts, with novel syntax.
'''

Comment by Thomas Kreuz:

This Python code (including all further comments) was written by Jeremy Fix (see http://jeremy.fix.free.fr/),
based on Matlab code written by Thomas Kreuz.

The SPIKE-distance is described in this paper:

Kreuz T, Chicharro D, Houghton C, Andrzejak RG, Mormann F:
Monitoring spike train synchrony.
J Neurophysiol 109, 1457-1472 (2013).

The Matlab codes as well as more information can be found at http://www.fi.isc.cnr.it/users/thomas.kreuz/sourcecode.html.

'''



import numpy as np
import pylab as p1

'''
Return the times (t1,t2) of the spikes in train[ibegin:]
such that t1 < t and t2 >= t
'''
def find_corner_spikes(t, train, ibegin, ti, te):
    if(ibegin == 0):
        tprev = ti
    else:
        tprev = train[ibegin-1]
    for idts, ts in enumerate(train[ibegin:]):
        if(ts >= t):
            return np.array([tprev, ts]), idts+ibegin
        tprev = ts
    return np.array([train[-1],te]), idts+ibegin


def bivariate_spike_distance(t1, t2, ti, te, N):
    '''Computes the bivariate SPIKE distance of Kreuz et al. (2012)
       t1 and t2 are 1D arrays with the spiking times of two neurons    
       It returns the array of the values of the distance
       between time ti and te with N samples.
       The arrays t1, t2 and values ti, te are unit less '''
    t = np.linspace(ti+(te-ti)/N, te, N)
    d = np.zeros(t.shape)

    t1 = np.insert(t1, 0, ti)
    t1 = np.append(t1, te)
    t2 = np.insert(t2, 0, ti)
    t2 = np.append(t2, te)

    # We compute the corner spikes for all the time instants we consider
    # corner_spikes is a 4 column matrix [t, tp1, tf1, tp2, tf2]
    corner_spikes = np.zeros((N,5))
 
    ibegin_t1 = 0
    ibegin_t2 = 0
    corner_spikes[:,0] = t
    for itc, tc in enumerate(t):
       corner_spikes[itc,1:3], ibegin_t1 = find_corner_spikes(tc, t1, ibegin_t1, ti, te)
       corner_spikes[itc,3:5], ibegin_t2 = find_corner_spikes(tc, t2, ibegin_t2, ti, te)

    #print corner_spikes
    xisi = np.zeros((N,2))
    xisi[:,0] = corner_spikes[:,2] - corner_spikes[:,1]
    xisi[:,1] = corner_spikes[:,4] - corner_spikes[:,3]
    norm_xisi = np.sum(xisi,axis=1)**2.0

    # We now compute the smallest distance between the spikes in t2 and the corner spikes of t1
    # with np.tile(t2,(N,1)) we build a matrix :
    # np.tile(t2,(N,1)) =    [   t2   ]        -   np.tile(reshape(corner_spikes,(N,1)), t2.size) = [                        ]
    #                        [   t2   ]                                                             [  corner  corner  corner]
    #                        [   t2   ]                                                             [                        ]
    dp1 = np.min(np.fabs(np.tile(t2,(N,1)) - np.tile(np.reshape(corner_spikes[:,1],(N,1)),t2.size)),axis=1)
    df1 = np.min(np.fabs(np.tile(t2,(N,1)) - np.tile(np.reshape(corner_spikes[:,2],(N,1)),t2.size)),axis=1)
    # And the smallest distance between the spikes in t1 and the corner spikes of t2
    dp2 = np.min(np.fabs(np.tile(t1,(N,1)) - np.tile(np.reshape(corner_spikes[:,3],(N,1)),t1.size)),axis=1)
    df2 = np.min(np.fabs(np.tile(t1,(N,1)) - np.tile(np.reshape(corner_spikes[:,4],(N,1)),t1.size)),axis=1)

    xp1 = t - corner_spikes[:,1]
    xf1 = corner_spikes[:,2] - t 
    xp2 = t - corner_spikes[:,3]
    xf2 = corner_spikes[:,4] - t

    S1 = (dp1 * xf1 + df1 * xp1)/xisi[:,0]
    S2 = (dp2 * xf2 + df2 * xp2)/xisi[:,1]

    d = (S1 * xisi[:,1] + S2 * xisi[:,0]) / (norm_xisi/2.0)
    
    # d is the  distance.
    # t is the time vector for plotting.  
    return t,d

def multivariate_spike_distance(spike_trains, ti, te, N):
    ''' t is an array of spike time arrays
    ti the initial time of the recordings
    te the end time of the recordings
    N the number of samples used to compute the distance
    spike_trains is a list of arrays of shape (N, T) with N spike trains
    The multivariate distance is the instantaneous average over all the pairwise distances
    '''
    d = np.zeros((N,))
    n_trains = len(spike_trains)
    t = 0
    for i, t1 in enumerate(spike_trains[:-1]):
        for t2 in spike_trains[i+1:]:
            tij, dij = bivariate_spike_distance(t1, t2, ti, te, N)
            if(i == 0):
                t = tij # The times are only dependent on ti, te, and N
            d = d + dij
    d = d / float(n_trains * (n_trains-1) /2)
    return t,d


#if __name__ == '__main__':
import matplotlib.pylab as p1
# We test the simulation of Kreuz(2012)

######################
# With 2 spikes trains


# We compute the multivariate SPIKE distance
list_spike_trains = [] # A misnomer since this is really a list of ISIs


#[list_spike_trains.append(h.cells.o(int(j)).isil.to_python()) for j in range(int(numcell-1))]
list_spike_trains= [ tvec[i] for i in xrange(0,np.max(gidvec))]

ti = 0
#tf = int(h.high_isi)# Below should be the same value
tf = np.max(np.array(list_spike_trains))

"""
t1 = np.arange(100, 1201, 100)
t2 = np.arange(100, 1201, 110)

"""




# Note just give the function bivariate_spike_distance ISI arguments, instead of times arguments.
#
#

'''
t, Sb = bivariate_spike_distance(h.cells.o(21).isil.to_python(), h.cells.o(23).isil.to_python(), ti, tf, 50)
print sum(Sb)

timespython=np.array(timespython).astype(int)
p1.figure(figsize=(20,6))



p1.subplot(211)
p1.hold(True)
#for i in range(np.shape(timespython)[1]-1):
#may need to use code from make_times, to detect when cells didn't fire.
# and add extra elements such that all isil vectors are the same as the isil instance that has the maximum number of elements (intervals [spike times]).
for i in xrange(int(len(h.cells.o(1).isil.to_python())-1)):
    p1.plot([int(h.cells.o(1).isil.x[i]),int(h.cells.o(1).isil.x[i])
], [0.5, 1.5], 'k')
for i in xrange(int(len(h.cells.o(3).isil.to_python())-1)):
    p1.plot([int(h.cells.o(3).isil.x[i]),int(h.cells.o(3).isil.x[i])], [1.5, 2.5], 'k')
p1.xlim([ti,tf])
p1.hold(False)
p1.ylim([0,3])
p1.title("Spike trains")

p1.subplot(212)
p1.plot(t, Sb,'k')
p1.xlim([ti,tf])
p1.ylim([0,1])
p1.xlabel("Time (ms)")
p1.title("Bivariate ISI distance")

sfin = 'kreuz_bivariateb_isi' + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
        + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
        + str(int(fin)) + '.png'
p1.savefig(sfin)
'''


#############################
# With multiple spikes trains
"""
num_trains = 50
num_spikes = 40 # Each neuron fires exactly 40 spikes
num_events = 5  # Number of events with increasing jitter
# spike_times is an array where each rows contains the spike times of a neuron
spike_times = np.zeros((num_trains, num_spikes))
# The first spikes are randomly spread in the first half of the simulation time
spike_times[:,range(num_spikes/2)] = tf/2.0 * np.random.random((num_trains, num_spikes/2))

# We now append the times for the events with increasing jitter
for i in range(1,num_events+1):
    tb = tf/2.0 * i / num_events 
    spike_times[:,num_spikes/2+i-1] = tb + (50 *(i-1) / num_events)* (2.0 * np.random.random((num_trains,)) - 1)


# And the second events with the decreasing jitter
num_last_events = num_spikes/2-num_events
for i in range(num_last_events):
    tb = tf/2.0 + tf/2.0 * i / (num_last_events-1)
    spike_times[:, -(i+1)] = tb + (50 - 50 *i / (num_last_events-1))* (2.0 * np.random.random((num_trains,)) - 1)


# Finally we sort the spike times by increasing time for each neuron
spike_times.sort(axis=1) 


"""


ti = 0
tf = int(h.high_isi)# Below should be the same value
tf = np.max(np.array(list_spike_trains))

#list_spike_trains=np.array(list_spike_trains).astype(int)
#print list_spike_trains
t, Sb = multivariate_spike_distance(list_spike_trains, ti, tf, 1000)
"""
print sum(sum(Sb)), 'sum of distances.'

Sbmatrix=np.matrix(Sb)
p1.imshow(Sbmatrix, interpolation='nearest')
#fig.savefig('distance matrix')
p1.show()
"""
# We plot the spike trains

"""
for i in range(np.shape(timespython)[1]-1):
    p1.plot([timespython[21][i], timespython[21][i]], [0.5, 1.5], 'k')
for i in range(np.shape(timespython)[1]-1):
    p1.plot([timespython[23][i], timespython[23][i]], [1.5, 2.5], 'k')
"""

"""

p1.figure(figsize=(20,6))
p1.subplot(211)
for i in range(int(h.times.size-1)): 
   for j in range(int(h.times[i].size-1)): 
       p1.hold(True)
        #nrnpython("i=int(0); j=int(1)")
       # nrnpython("p1.plot([a,a],[i, i+1],'k');p1.show()")
       a=h.times[int(i)].x[int(j)] 
       p1.plot([a,a],[i, i+1],'k')
      
      
p1.hold(False)      
p1.show()
"""
fig=p1.figure()
fig.clf()


p1.figure(figsize=(20,6))
#p1.subplot(211)
p1.hold(True)

for j in xrange(int(ncell-1)): #range arguments are lists.
    for i in xrange(int(h.cells.o(j).isil.size())):

        #p1.plot([timespython[j][i], timespython[j][i]],[i, i+1],'k')
        #
        p1.plot([h.cells.o(j).isil.x[i], h.cells.o(j).isil.x[i]],[j, j+1],'k')
        print [h.cells.o(j).isil.x[i], h.cells.o(j).isil.x[i]],[j, j+1]
        #int(h.cells.o(j).isil.x[i])
p1.hold(False)
p1.title('Inter Spike Intervals')
"""
p1.subplot(212)
p1.plot(t,Sb,'k')
p1.xlim([0, tf])
p1.ylim([0, 1])
p1.xlabel("Time (ms)")
"""

#title_str=''
title_str='Coefficient of Variation'+' '+str(h.get_CV())
p1.ylabel("Multivariate_ISI_distance")
p1.xlabel("Cell number")

p1.title(title_str)

sfin = 'kreuz_multivariatem_isi' + str(int(h.prunenet)) +str(int(h.iw)) + str(tstop) + str(h.plastic) \
        + str(int(h.ff)) + str(numcell) + str(int(h.minr)) \
        + str(int(fin)) + '.png'
        
p1.savefig(sfin)
#

"""
"""
