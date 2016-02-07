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
import d3py
import pickle



bp = BiophysicalPerisomaticApi('http://api.brain-map.org')


#Above and below are the same. Above is online version, below is offline version that acts on local files.
#Use the other methods in the biophysical_perisomatic_api to one by one do the things needed to cache whole models.

'''
def read_local_json():
    from allensdk.api import api
    f1= open('neuron_models_from_query_builder.json')
    information = api.load(f1)
    return information
information=read_local_json()
'''

#bp.cache_data(395310469, working_directory='neuronal_model')
#for i in information['msg']:    
#    bp.cache_data(i['id'], working_directory='neuronal_model')

#==============================================================================
#Excitatory and inhibitory cells are categorised according to the suffix's listed below:
#Pvalb:  Inhibitory aspiny cells (i.e. it clustered the fast-spiking Pvalb cells)
#Scnn1a: Layer 4 excitatory pyramidal neurons.
#CRE virus connectivity?
from utils import Utils
config = Config().load('config.json')
utils = Utils(config,NCELL=4000)
info_swc=utils.gcs(utils.NCELL)
'''
utils.wirecells()#wire cells on different hosts.
utils.matrix_reduce()
utils.h('forall{ for(x,0){ uninsert xtra}}')    


from rigp import NetStructure
hubs=NetStructure(utils.ecm,utils.icm,utils.celldict)
print 'experimental rig'
hubs.hubs()

if utils.COMM.rank==0:
    utils.plotgraph()
    hubs.save_matrix()

#In addition to stimulating the out degree hub, stimulate the first cell on each host,
#To make activity more likely.
utils.setup_iclamp_step(utils.cells[0], 0.27, 1020.0, 750.0)

# configure recording
utils.spikerecord()

vec = utils.record_values()

#utils.h.dt = 0.025
#util.h.xopen("rigp.hoc")
tstop = 3000
utils.COMM.barrier()


utils.prun(tstop)


#with open(str(utils.COMM.rank)+'vectors.p', 'wb') as handle:
#       pickle.dump(vec, handle)    

#system.
chtodir=os.getcwd()+"../../tigramite_1.3"
os.chdir("/home/russell/tigramite_1.3")
#if COMM.rank==0:
#execfile('tigramite_gui.py')
#utils.prun(10)
#a=tigramite.Tigramite()

#import http_server as hs
#hs.load_url('force.json')

#morphs,swclist,cells1=read_local_swc()        
#cells1
def mkjson(): #Only ascii as in dictionary contents
    from allensdk.core.json_utilities import write

    #can be serialised into dictionary contents.
    for m in morphs:
        write(str(m)+'.json',m.root)
    return 0

'''