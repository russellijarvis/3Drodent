# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 17:01:22 2015

@author: russell
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
utils = Utils(config,NCELL=40,readin=0)
info_swc=utils.gcs(utils.NCELL)
utils.wirecells()#wire cells on different hosts.
utils.matrix_reduce()
utils.h('forall{ for(x,0){ uninsert xtra}}')   #mechanism only needed for wiring cells not for simulating them. 
from rigp import NetStructure
if utils.COMM.rank==0:
    hubs=NetStructure(utils,utils.my_ecm,utils.my_icm,utils.visited,utils.celldict)
    print 'experimental rig'
    #utils.plotgraph()
    hubs.save_matrix()
    hubs.hubs()    
    hubs.insert_cclamp(hubs.outdegree,hubs.indegree)
    utils.dumpjsongraph()

hubs=NetStructure(utils,utils.ecm,utils.icm,utils.visited,utils.celldict)
hubs.insert_cclamp(hubs.outdegree,hubs.indegree)
hubs.insert_cclamp(0,1)
#utils.graph_reduce()
utils.spikerecord()
vec = utils.record_values()
print 'setup recording'
tstop = 1150
utils.COMM.barrier()
utils.prun(tstop)
import matplotlib 
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
fig = plt.figure()
fig.clf()
plt.hold(True) #seems to be unecessary function call.
#TODO outsource management of membrane traces to neo/elephant.
#TODO use allreduce to reduce python dictionary to rank0
#This is different to Allreduce.
for gid,v in vec['v'].iteritems():
    #print v.to_python()
    plt.plot(vec['t'].to_python(),v.to_python())
fig.savefig('membrane_traces'+str(utils.COMM.rank)+'.png')    
plt.hold(False) #seems to be unecessary function call.
plt.xlabel('time (ms)')
plt.ylabel('Voltage (mV)')
plt.title('traces')
plt.grid(True)

print type(utils.tvec.to_python())
print type(utils.gidvec.to_python())
def plot_raster(tvec,gidvec):
    fig = plt.figure()
    fig.clf()
    colors=np.array([[0.42,0.67,0.84],[0.50,0.80,1.00],[0.90,0.32,0.00],[0.34,0.67,0.67],[0.42,0.82,0.83],[0.90,0.59,0.00],[0.33,0.67,0.47],[0.42,0.83,0.59],[0.90,0.76,0.00],[1.00,0.85,0.00],[0.71,0.82,0.41],[0.57,0.67,0.33],[1.00,0.38,0.60]]) # Colors for each cell population
    plt.title("Raster Plot")
    plt.hold(True)
    j=len(colors)-1
    plt.plot(tvec,gidvec,'.',c=colors[j], markeredgecolor = 'none')
    plt.savefig('raster'+str(utils.COMM.rank)+'.png')

plot_raster(utils.tvec.to_python(),utils.gidvec.to_python())
    

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