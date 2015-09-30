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
from mpi4py import MPI

# initialize the MPI interface
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

PROJECTROOT=os.getcwd()

bp = BiophysicalPerisomaticApi('http://api.brain-map.org')


#Above and below are the same. Above is online version, below is offline version that acts on local files.
#Use the other methods in the biophysical_perisomatic_api to one by one do the things needed to cache whole models.
def read_local_json():
    from allensdk.api import api
    f1= open('neuron_models_from_query_builder.json')
    information = api.load(f1)
    return information
information=read_local_json()

#bp.cache_data(395310469, working_directory='neuronal_model')
#for i in information['msg']:    
#    bp.cache_data(i['id'], working_directory='neuronal_model')

#==============================================================================
#Excitatory and inhibitory cells are categorised according to the suffix's listed below:
#Pvalb:  Inhibitory aspiny cells (i.e. it clustered the fast-spiking Pvalb cells)
#Scnn1a: Layer 4 excitatory pyramidal neurons.


NCELL=18    
import glob
from allensdk.model.biophysical_perisomatic.utils import Utils
from allensdk.model.biophys_sim.config import Config


#import brain_functions as bf
from neuron import h

h('objref py, pc')
h('py = new PythonObject()')
h('pc = new ParallelContext()')



def prep_list():                
    import pickle
    allrows = pickle.load(open('allrows.p', 'rb'))
    allrows.remove(allrows[0])#The first list element are the column titles. 
    allrows2 = [i for i in allrows if int(len(i))>9 ]
    return allrows2 


allrows2 = prep_list()#NFILE)
gidvec = []
s1=''

os.chdir(os.getcwd() + '/main')        

#h.cells, allrows2, ie0, ie1 =bf.mb(RANK, NCELL, SIZE, allrows2, gidvec,h,s1)

#info_swc=utils.gcs(utils.NCELL)
#utils.seclists()
from utils import Utils

config = Config().load('config.json')
utils = Utils(config)
NCELL=utils.NCELL=18

from neuron import h

h('{load_file("stdgui.hoc")}')
h('{load_file("nrngui.hoc")}')
h('load_file("import3d.hoc")')

info_swc=utils.gcs(utils.NCELL)
nclist, ecm, icm=utils.wirecells_s()#Wire cells on same host.
nclist, ecm, icm=utils.wirecells4()#wire cells on different hosts.
h('forall{ for(x,0){ uninsert xtra }}') 

def progress(pinvl, swlast):
#From http://senselab.med.yale.edu/ModelDB/ShowModel.asp?model=151681&file=\bulb3d\util.py
  sw = h.startsw()
  print "t=%g wall interval %g"% (h.t, sw-swlast)
  h.cvode.event(h.t+pinvl, (progress, (pinvl , sw)))

def show_progress(invl):
#From http://senselab.med.yale.edu/ModelDB/ShowModel.asp?model=151681&file=\bulb3d\util.py

  global fih
  if RANK == 0:
    fih = h.FInitializeHandler(2, (progress, (invl, h.startsw())))

show_progress(200)
prun(tstop)
h.xopen('tools.hoc')
h.xopen('analysis2.hoc')

execfile('analysisp3.py')
h('time_finish=pc.time()')
h('print time_finish')

import numpy as np
print 'sums of ecm and icm'
print np.sum(ecm), np.sum(icm)
print SIZE, RANK
import mpi4py as mpi4py
mpi4py.get_config(), 'config'
pc = utils.h.ParallelContext()
s = "mpi4py thinks I am %d of %d,\
# NEURON thinks I am %d of %d\n"
cw = MPI.COMM_WORLD
print s % (cw.rank, cw.size, pc.id(), pc.nhost())
time_start=pc.time()
print time_start




def read_local_swc():
    morphs=[]
    cells1=[]
    Utils.initialize_hoc()
    swclist=glob.glob('*.swc')
    for swcf in swclist:
        #morphology = swc.read_swc(swcf)
        morphology=Utils.generate_morphology(swcf)
        cell1=Utils.load_cell_parameters()       
        cells1.append(cell1)
        print type(cells1)
        print type(cell1)
        morphology.root
        morphs.append(morphology)
    return morphs,swclist,cells1

#morphs,swclist,cells1=read_local_swc()        
#cells1
def mkjson(): #Only ascii as in dictionary contents
    from allensdk.core.json_utilities import write

    #can be serialised into dictionary contents.
    for m in morphs:
        write(str(m)+'.json',m.root)
    return 0

