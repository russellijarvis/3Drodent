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
import pdb
config = Config().load('config.json')
utils = Utils(config,NCELL=20,readin=1)
info_swc=utils.gcs(utils.NCELL)

utils.wirecells_test()#wire cells on different hosts.
utils.matrix_reduce()
#utils.graph_reduce()
utils.h('forall{ for(x,0){ uninsert xtra}}')    
from rigp import NetStructure
if utils.COMM.rank==0:
    hubs=NetStructure(utils,utils.my_ecm,utils.my_icm,utils.visited,utils.celldict)
    print 'experimental rig'
    utils.plotgraph()
    hubs.save_matrix()
    hubs.hubs()
    print '\n', utils.COMM.rank, hubs.outdegree, " outdegree"
#Does the insertion of an IClamp work 
hubs=NetStructure(utils,utils.ecm,utils.icm,utils.visited,utils.celldict)
hubs.hubs()
print '\n', utils.COMM.rank, hubs.outdegree, " outdegree"
#In addition to stimulating the out degree hub, stimulate the first cell on each host,
#To make activity more likely.
utils.setup_iclamp_step(utils.cells[0], 0.27, 1020.0, 750.0)
# configure recording
utils.spikerecord()
vec = utils.record_values()
print 'setup recording'
tstop = 1500
utils.COMM.barrier()
utils.prun(tstop)
tvec=utils.tvec.to_python()
idvec=utils.idvec.to_python()
#utils.vec_reduce()
print tvec
print idvec

import matplotlib 
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
fig = plt.figure()
fig.clf()

for gid,v in vec['v'].iteritems():
    print v.to_python()
    plt.plot(vec['t'].to_python(),v.to_python())
fig.savefig('membrane_traces'+str(utils.COMM.rank)+'.png')    

plt.xlabel('time (ms)')
plt.ylabel('Voltage (mV)')
plt.title('traces')
plt.grid(True)

utils.vec_reduce()
if utils.RANK==0:
    vecs=zip(utils.my_tvec,utils.my_idvec)
with open(str(utils.COMM.rank)+'vectors.p', 'wb') as handle:
    pickle.dump(utils.my_tvec, handle)    
#system.
chtodir=os.getcwd()+"../../tigramite_1.3"
os.chdir("/home/russell/tigramite_1.3")
#if COMM.rank==0:
#execfile('tigramite_gui.py')
#utils.prun(10)
#a=tigramite.Tigramite()

#import http_server as hs
#hs.load_url('force.json')

'''

def mkjson(): #Only ascii as in dictionary contents
    from allensdk.core.json_utilities import write

    #can be serialised into dictionary contents.
    for m in morphs:
        write(str(m)+'.json',m.root)
    return 0


#morphs,swclist,cells1=read_local_swc()        
#cells1


#bp.cache_data(395310469, working_directory='neuronal_model')
#for i in information['msg']:    
#    bp.cache_data(i['id'], working_directory='neuronal_model')

#f1 = open('/home/russell/git/allen/neuron_models_from_query_builder.json')
#f1= open('neuron_models_from_query_builder.json')
#information = api.load(f1)
#return information
'''
#def read_local_json():
#    from allensdk.api import api
#    api.json_utilities.read('/home/russell/git/allen/neuron_models_from_query_builder.json')
#read_local_json()
#pdb.set_trace()

