# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 17:01:22 2015

The parallel wiring related functions are written by Russell Jarvis rjjarvis@asu.edu
"""
import allensdk
from allensdk.api.queries.biophysical_perisomatic_api import \
    BiophysicalPerisomaticApi

from allensdk.api.queries.cell_types_api import CellTypesApi
import allensdk.core.swc as swc
import os
from allensdk.core.nwb_data_set import NwbDataSet
import glob
from allensdk.model.biophysical_perisomatic.utils import Utils
from allensdk.model.biophys_sim.config import Config
#import d3py
import pickle
bp = BiophysicalPerisomaticApi('http://api.brain-map.org')
#import unittest
from utils import Utils
import numpy as np
import pdb 
config = Config().load('config.json')
# The readin flag when set enables the wiring to be read in from pre-existing 
# pickled files with rank specific file names.
utils = Utils(config,NCELL=40,readin=1)
info_swc=utils.gcs(utils.NCELL)
utils.wirecells()#wire cells on different hosts.
utils.matrix_reduce()
utils.h('forall{ for(x,0){ uninsert xtra}}')   #mechanism only needed for wiring cells not for simulating them. 
from rigp import NetStructure
if utils.COMM.rank==0:
    hubs=NetStructure(utils,utils.global_ecm,utils.global_icm,utils.visited,utils.celldict)
    print 'experimental rig'
    #utils.plotgraph()
    hubs.save_matrix()
    #if utils.COMM.rank==0:
    #
    # A global analysis of hub nodes, using global complete adjacency matrices..
    #
    hubs.hubs()    
   
    amplitude=0.27 #pA or nA?
    delay=50 # was 1020.0 ms, as this was long enough to notice unusual rebound spiking
    duration=100.0 #was 750 ms, however this was much too long.

    hubs.insert_cclamp(hubs.outdegree,hubs.indegree,amplitude,delay,duration)
    utils.dumpjsongraph()


hubs=NetStructure(utils,utils.ecm,utils.icm,utils.visited,utils.celldict)
#
# A local analysis of hub nodes, using local incomplete adjacency matrices.
#
hubs.hubs()
amplitude=0.27 #pA or nA?
delay=50 # was 1020.0 ms, as this was long enough to notice unusual rebound spiking
duration=5.0 #was 750 ms, however this was much too long.

hubs.insert_cclamp(hubs.outdegree,hubs.indegree,amplitude,delay,duration)
vec = utils.record_values()
print 'setup recording'
#tstop=20
tstop=10
#tstop = 2150
utils.COMM.barrier()
utils.prun(tstop)

utils.global_vec = utils.COMM.gather(vec,root=0) # Results in a list of dictionaries on rank 0 called utils.global_vec
# Convert the list of dictionaries into one big dictionary called global_vec (type conversion).
if utils.COMM.rank==0:
    utils.global_vec = {key : value for dic in utils.global_vec for key,value in dic.iteritems()  } 


if utils.COMM.rank==0:        
    import matplotlib 
    import matplotlib.pyplot as plt
    matplotlib.use('Agg') 
    fig = plt.figure()
    fig.clf()
    plt.hold(True) #seems to be unecessary function call.
    #TODO outsource management of membrane traces to neo/elephant.
    #TODO use allreduce to reduce python dictionary to rank0
    #This is different to Allreduce.
    for gid,v in utils.global_vec['v'].iteritems():
        #print v.to_python()
        plt.plot(utils.global_vec['t'].to_python(),v.to_python())
    fig.savefig('membrane_traces_from_all_ranks'+str(utils.COMM.rank)+'.png')    
    plt.hold(False) #seems to be unecessary function call.
    plt.xlabel('time (ms)')
    plt.ylabel('Voltage (mV)')
    plt.title('traces')
    plt.grid(True)
    
utils.spike_reduce() #Call matrix_reduce again in order to evaluate global_spike
if utils.COMM.rank==0:        
    #tvec and gidvec are Local variable copies of utils instance variables.
    tvec=np.zeros_like(np.array(utils.tvec.to_python))
    gidvec=np.zeros_like(np.array(utils.gidvec.to_python))

    for i,j in utils.global_spike:# Unpack list of tuples
        #create local variables.
        i=np.array(i)
        j=np.array(j)
        tvec.extend(i)
        gidvec.extend(j)
        print utils.tvec
        print utils.gidvec
    #utils.global_spike(utils.gidvec,utils.tvec)
    
    utils.dumpjsongraph()


#tvec=utils.h.tvec.to_python()
#gidvec=utils.h.gidvec.to_python()
#print type(tvec)
#print type(gidvec)
def plot_raster(tvec,gidvec):
    pallete=[[0.42,0.67,0.84],[0.50,0.80,1.00],[0.90,0.32,0.00],[0.34,0.67,0.67],[0.42,0.82,0.83],[0.90,0.59,0.00], 
                [0.33,0.67,0.47],[0.42,0.83,0.59],[0.90,0.76,0.00],[1.00,0.85,0.00],[0.71,0.82,0.41],[0.57,0.67,0.33]]

    fig = plt.figure()
    fig.clf()
    color=[1.00,0.38,0.60] # Choose differe Colors for each cell population
    plt.title("Raster Plot")
    plt.hold(True)
    plt.plot(tvec,gidvec,'.',c=color, markeredgecolor = 'none')
    plt.savefig('raster'+str(utils.COMM.rank)+'.png')
    
if utils.COMM.rank==0:
    plot_raster(tvec,gidvec)
    

#Probably just get the spike distance.

#import http_server as hs
#hs.load_url('force.json')

'''

def mkjson(): #Only ascii as in dictionary contents
    from allensdk.core.json_utilities import write

    #can be serialised into dictionary contents.
    for m in morphs:
        write(str(m)+'.json',m.root)
    return 0


'''
