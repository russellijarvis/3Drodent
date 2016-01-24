#import mpi4py  
#from mpi4py import MPI
#import pdb
#pdb.set_trace()

import neuron
from neuron import h


import os
import neuron 
from numpy import *
import scipy
import numpy as np
from neuron import h #//Can call Python from HOC, and this line allows, to call HOC from python, allowing nesting of calls. 
#//Now the below imports are duplicated once here and once in import_data2.py
#//will this destroy functionality?


import matplotlib
matplotlib.use('Agg') 


#import neuron
#import hoc
#h=hoc.HocObject()
# There is problem with passing tstop in as a variable. IT has to be hard coded, in its most native env. 
# There is no evidence that other variables are not successfuly being passed in however.

# The problem seems to be that tstop needs to be defined in python not in hoc.

h.dt = 0.025
tstop = 1000
h("dc=-1")#delete cell. This retains the cell, but makes all of its section diameter 0us, and it makes all projection weights to it 0.
h("plastic=1")
h("ncell=1182/2")#minimum 17
#can take 100 now.
h("minr=5.3")
#h("minr=6.3")
#was 7.3
h("prunenet=25")
h("offset=0")
h("parallels=0")
#h("tstop=140//1000//1 second =1000 ms")
h("delay=3000")
h("iw=1000 // divisor/ or scaler, so if its bigger weight should be smaller. Was 1000 but this should be stronger. delay equals internal weight.")

h("iwp=1")
h("ew=0.0185// External weight.This value very appropriate. But perhaps smaller how to automate the tuning of synaptic weights?")
h("run_iter=1")

h("msdp=0 //basically prune some connections distal from the soma.")

h("fast_wire=0")
h("ff=0")
h("if(ff==1){ ncell=40 }")
h("if(ff==1){ prunenet=0 }")
h("numcell=ncell")#should really use sed. Such as to reduce duplicate variable to one.
h("no_b=0")
h("lpy=0")



#Although run is ultimately called from inside NEURON, I have created NEURON from inside Python and not the other way around



# run the simulation





lfp=1
h.xopen("init3.hoc")



def runp(tstop):
 if lfp==1:
  h.init()

 h.finitialize() #set initial variables.
 while h.t<tstop: #Perform integration over the cable equation.
     if lfp==1:
       h.advance() #this method is designed to work with field and may take longer.
     if lfp!=1:
       h.fadvance()



runp(tstop)

"""
if lfp!=1:
 #This method does not generate any tvec or idvec. Possibly only for 
 # feed forward.
 h('dt=0.0025')
 h('tstop=1000')
 h('run()')
"""

h.xopen("post_analysis4.hoc")
h.summary()
#execfile('fft.py')
execfile('plotting7.py')
h.xopen('fft.hoc')
execfile('fft.py')

h.summary()
#h('if(strcmp(machinename,"zaza3")==0){nrnpython("execfile('fmri.py')")}')
if h.machinename=='zaza3':
 execfile('fmri.py')

#execfile('fmri.py')
pe=nx.pagerank(dirg)
ps=nx.pagerank(dirg2)
prs=highest_centrality(ps)
pre=highest_centrality(pe)
print 'summary 3:'
print 'pagerank, sgc'
print prs
print 'pagerank, ent'
print pre

cells=h.cells
allsectionsinpython=h.allsec()
#for sec in h.allsec(): print sec.name()                                                                                                                                  

