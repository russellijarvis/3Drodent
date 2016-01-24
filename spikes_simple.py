# This buffer is for notes you don't want to save, and for Lisp evaluation.
# If you want to create a file, visit that file with C-x C-f,
# then enter the text in that file's own buffer.

import os
import glob
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from mpi4py import MPI
#GPIO pins on microcontoller are both TX and RX
#import neuron
from neuron import h
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


#allmemb=[]
'''
Run first time after simulation. This should be on the end of the simulation code.
for fname in glob.glob('*soma*.dat_'):
    z = np.loadtxt(fname)
    allmemb.append(z)
    print fname
pickle.dump(allmemb, open("allmemb.p", "wb" ) )
'''
#allmemb=pickle.load(open('allmemb.p','rb'))

x=np.loadtxt('_spt.dat')
'''
h('objref cell')
h('print cell')

h('{load_file("stdgui.hoc")}')
h('{load_file("nrngui.hoc")}')

h('load_file("import3d.hoc")')

h.xopen("morph4.hoc")
os.chdir(os.getcwd()+'/main')


for i in range(RANK, int(np.shape(allmemb)[0]), SIZE): #20 was int(len(allrows))

    
    s=allmemb[i]    
    h.cell=h.mkcell('cell-93-trace.CNG.swc')
    
    h('cell.gid1=py.i')
    h('cell.recvec=cell.recvec.from_python(py.s)')
    print i, RANK
    h('print cell.recvec.sum()')
    for j in x:
'''
