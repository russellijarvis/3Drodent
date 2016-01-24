


#!/usr/bin/env python

# All original code, that appeals to idioms described by Hines and Carnevale.

# Copyright (C) 2012, 2013 Russell Jarvis
# This file is part of This.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License

import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
import os
import urllib2
import zipfile
#import mpi4py
from mpi4py import MPI
#GPIO pins on microcontoller are both TX and RX
#import neuron
from neuron import h
#import LFPy#Importing LFPy imports neuron. NEURON has to be imported after mpi4py


#initialize the MPI interface
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


h('objref pc')
pc = h.ParallelContext()
h.pc=pc
#h('pc = new ParallelContext()')
#h('pc = py.pc')
s = "mpi4py thinks I am %d of %d,\
 NEURON thinks I am %d of %d\n"
cw = MPI.COMM_WORLD
print s % (cw.rank, cw.size, \
           pc.id(),pc.nhost())
#pc.done()

#h.dt = 0.5
#tstop = 1025
#h('print dt')
h('pc.set_maxstep(10)')
h('load_file("nrngui.hoc")')
h('stdinit()')
h('pc.psolve(1000)')# Interestingly this does cause an error its just non fatal.
