#import mpi4py  
#from mpi4py import MPI
import neuron
from neuron import h
#import neuron
#import hoc
#h=hoc.HocObject()
# There is problem with passing tstop in as a variable. IT has to be hard coded, in its most native env. 
# There is no evidence that other variables are not successfuly being passed in however.

# The problem seems to be that tstop needs to be defined in python not in hoc.


h.dt = 0.025
tstop = 1000.0
h("ncell=100")
h("minr=30")
h("prunenet=0")
h("parallels=0")
#h("tstop=140//1000//1 second =1000 ms")
h("delay=3000")
h("iw=1000 // divisor/ or scaler, so if its bigger weight should be smaller. Was 1000 but this should be stronger. delay equals internal weight.")

h("iwp=1")
h("ew=0.0185// External weight.This value very appropriate. But perhaps smaller how to automate the tuning of synaptic weights?")
h("run_iter=1	 ")

h("msdp=0 //basically prune some connections distal from the soma.")

h("fast_wire=0")
h("ff=0	")
h("no_b=0")

h("plastic=0")
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
"""


#Although run is ultimately called from inside NEURON, I have created NEURON from inside Python and not the other way around



# run the simulation






h.xopen("init3.hoc")

def runp():
 h.finitialize() #set initial variables.
 while h.t<tstop: #Perform integration over the cable equation.
     h.fadvance()
runp()
h.xopen("post_analysis4.hoc")
execfile('fmri.py')
execfile('plotting3.py')

cells=h.cells
allsectionsinpython=h.allsec()
for sec in h.allsec(): print sec.name()

