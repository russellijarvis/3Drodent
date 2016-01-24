from mpi4py import MPI
import os
import neuron
from neuron import h
#from nrn import h


pc = h.ParallelContext()

s = "mpi4py thinks I am %d of %d, NEURON thinks I am %d of %d\n"
cw = MPI.COMM_WORLD
print s % (cw.rank, cw.size, pc.id(), pc.nhost())
h("parallels=1")
h.xopen('initp.hoc')
pc.done()
