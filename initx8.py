#!/usr/bin/python
# -*- coding: utf-8 -*-

import neuron
from neuron import h

import os
import neuron
from numpy import *
import scipy
import numpy as np

import matplotlib
matplotlib.use('Agg')
import pylab as p1

import networkx as nx

# I fixed a bug in the spatial wiring procedure which should allow the program to be satisfied. more easily.

h.dt = 0.025
tstop = 2500
h('dc=-1')  # delete cell. This retains the cell, but makes all of its section diameter 0us, and it makes all projection weights to it 0.
h('plastic=1')
h('get_dist=0')
h('dlc=0')  # dilate centres.
structural_plotting=1
h('ncell=200')
#I normaly use alphabet size of two symbols (0,1). There is no reason why this could not be increased. Using multilevel quantisation.
prunenet=0
h('qd=0')#quick and dirty for debugging
h('lpy=1')#declare lpy here, the value is not important, its only important that the declaration exists.

   
h('if(qd==1){prunenet=900}')
h('objref py')
#h('objref py')
#h('py = new PythonObject()')
#h('if(qd==1){ py.structural_plotting=0')
h('minr=10')#minumum r large, for some reason which I do not yet 
h('minrs=10')#minumum r small, understand
#the number of connections made, depends somewhat on the last run, this may be 
# caused by morph3.hoc, deleting inappropriate cells occasionally.

# h("offset=0")
##OLD

h('offset=0')

##

h('parallels=0')

# h("tstop=1000//1 second =1000 ms")

h('delay=3000')
h('iw=2500 // divisor/ or scaler, so if its bigger weight should be smaller. Was 1000 but this should be stronger. delay equals internal weight.'
  )

h('iwp=1')
h('ew=0.0185')
#// External weight.This value very appropriate. But perhaps smaller how to automate the t#uning of synaptic weights?'
#  )
h('run_iter=1')

h('msdp=0 //basically prune some connections distal from the soma.')
h('fast_wire=0')
h('ff=0')
h('if(ff==1){ get_dist=0 }')
h('if(ff==1){ ncell=40 }')
h('if(ff==1){ prunenet=0 }')
h('numcell=ncell')  # should really use sed. Such as to reduce duplicate variable to one.
ncell=h.ncell
h('no_b=0')

h('if (numcell>99){large=1}')
h('if (numcell>99){large_scale=1}')

# Although run is ultimately called from inside NEURON, I have created NEURON from inside Python and not the other way around.

# run the simulation


#I think the reason my simulation takes up more than 60 cells worth of memory is also because of all the variables I have added.
lfp = 1
h.xopen('init3.hoc')


# execfile('structure.py')
# This file contents is now in find_in_degree()
#execfile('structure2.py')
def runp(tstop):
    if lfp == 1:
        h.init()
    if lfp != 1:
    	h.finitialize()  # set initial variables.
    while h.t < tstop:  # Perform integration over the cable equation.
        if lfp == 1:
            h.advance()  # this method is designed to work with field and may take longer.
        if lfp != 1:
            h.fadvance()


runp(tstop)



fin=0
execfile('structure2.py')
sfin0 = str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
    + str(int(h.ff)) + str(ncell) + str(int(fin)) + str(int(outdegree)) \
    + str(int(outdegreeg))

if h.machinename != 'zaza3':
 
    h.xopen("post_analysis6.hoc")
    #execfile('plotting9.py')
    execfile('spike_distance.py')# somehow the 'times[j].x[i]' object ie the array of vectors, gets destroyed. 
    # Probably by altering a reference to it, so its not available.
    h.xopen("isi.hoc")

    execfile('isi_distance.py')
    h.xopen('fft.hoc')
    execfile('fft2.py')
    h.xopen("post_analysis6.hoc")
    execfile('plotting9.py')
    h.isi_summary()
    h.summary()
    hoc_string = 'hoc_stdout("summary_txt' + sfin0 + '")'
    h(hoc_string)
    h.summary()
    h.isi_summary()
	#summary()
    #print 'fourier, sum of power', sum(lfp1s), ' ', sum(vfr)

    h('hoc_stdout()')

    h.xopen("supplament.hoc")
    execfile('supplament.py')    
    quit()
    # os.system('cd -')

if h.machinename == 'zaza3':
    h.xopen("post_analysis6.hoc")
    execfile('plotting9.py')
    execfile('spike_distance.py')
    h.xopen("isi.hoc")

    execfile('isi_distance.py')
    h.xopen('fft.hoc')
    execfile('fft2.py')

    hoc_string = 'hoc_stdout("summary_txt' + sfin0 + '")'
    h(hoc_string)
    h.summary()
    h.ent_table_sq2()
    execfile('plotting9.py')
    print 'fourier, sum of power', sum(lfp1s), ' ', sum(vfr)
    h('hoc_stdout()')
    ##Cannot load coherance on trifid because of unsupported packages.
    h.xopen("supplament.hoc")
    execfile('coherance.py')

    
    

for sec in h.allsec(): print sec.name()

if h.machinename != 'zaza3':
    print 'The program has executed cleanly'
cells = h.cells
allsectionsinpython = h.allsec()
help(h.Vector)
#quit()

