import os
import glob

import pickle
import pdb

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from mpi4py import MPI
from neuron import h
from bsmart import granger  # Load the Granger calculation tool
from matplotlib.colors import LogNorm
#import natsort
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

pc = h.ParallelContext()
h('objref pc')
h.pc=pc
s = "mpi4py thinks I am %d of %d,\
 NEURON thinks I am %d of %d\n"
cw = MPI.COMM_WORLD
print s % (cw.rank, cw.size, \
           pc.id(),pc.nhost())
h('time_start=pc.time()')
h('objref py')
h('py = new PythonObject()')




import re

#Natural sort is different to conventional comp alg sort.
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

#This code will add time. I can comment it out later,
#once it has executed on trifid properly.


#Run first time after simulation. This should be on the end of the simulation code.

sortedlist=glob.glob('*.npy')#key=int
sortedlist=natural_sort(sortedlist)
allmemb=[]
#size_variable=[]
#size_variable=pickle.load(open('allmemb.p', 'rb' ) )
sizev=np.load(sortedlist[0])

#fsize=len(size_variable)

result_matrix=np.zeros_like(sizev)      


allmemb=[]
for fname in sortedlist:#sorted(glob.glob('*soma*.dat_')):
    z = np.load(fname)
    np.add(result_matrix,z,result_matrix)
    allmemb.append(z)
    print fname
    print np.sum(z)
    print np.sum(result_matrix)

      #cnt+=1


fig = plt.figure()
fig.clf()
im = plt.imshow(result_matrix, interpolation='nearest')
plt.autoscale(True)
plt.colorbar(im)
plt.xlabel('columns = targets')
plt.ylabel('rows = sources')
plt.title('Transfer Entropy Matrix')
plt.grid(True)
#sfin = str(starti)+str(stopi)+str(startj)+str(stopj)+str(flag)+'Transfer_Entropy_matrix2.png'
fig.savefig('final_nte.png')
#np.save(result_matrix)
