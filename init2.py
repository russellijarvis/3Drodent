#import mpi4py  
#from mpi4py import MPI
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
h("plastic=0")
h("ncell=200")#minimum 17
h("minr=30")
h("prunenet=50")
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
h("numcell=ncell")#should really use sed. Such as to reduce duplicate variable to one.
h("no_b=0")
h("lpy=0")

"""

h("entpy=1") # I have worked around the problem that I am launching from both Python and NEURON 

h.entpy=1

h("delay=iw=1000") #//delay equals internal weight.
h("iwp=1") #delay should not include the somatic distance. only r, and if fast wire its actual 
#Distance. Propogation delay is radius dependent even in single neurons.
h("ew=0.0185// External weight.This value very appropriate. But perhaps smaller how to automate the tuning of synaptic weights?")
h("run_iter=1")
h("prunenet=9000 //basically prune some connections distal from the soma.")
#//data=1
h("fast_wire=0")
h("ff=0")
h("ncell=60")
h("ninhib=3*ncell/4")
h("nexec=ncell-ninhib")
h("minr=3")
h("plot_syn=0")
h("get_dist=0")
# //to scale the delay
h("plastic=1")
"""


#Although run is ultimately called from inside NEURON, I have created NEURON from inside Python and not the other way around



# run the simulation






h.xopen("init3.hoc")
"""
def runp(tstop):
 h.finitialize() #set initial variables.
 while h.t<tstop: #Perform integration over the cable equation.
     h.fadvance()
runp(tstop)
h.xopen("post_analysis4.hoc")
p1.show()
"""
h('tstop=1000')
h('run()')


"""
stim_one(outdegree)
stim_one(indegree)
stim_one_inter(outdegreeg)
stim_one_inter(indegreeg)
runp(tstop)
h.xopen("post_analysis4.hoc")

"""
#execfile('plotting4.py')

execfile('plotting7.py')
execfile('fft.py')
#execfile('fmri.py')
cells=h.cells
allsectionsinpython=h.allsec()
for sec in h.allsec(): print sec.name()

