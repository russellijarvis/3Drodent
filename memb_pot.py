# This buffer is for notes you don't want to save, and for Lisp evaluation.
# If you want to create a file, visit that file with C-x C-f,
# then enter the text in that file's own buffer.

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

class acell:
   'Common base class for analysis of cells'
   cellcnt = 0

   def __init__(self,gid, membv,spike):
      self.spike = spike
      self.membv= membv#np.zeros(0,int(tstop/mindelay))
      self.gid = gid
      
      acell.cellcnt += 1

allmemb=[]



import re

#Natural sort is different to conventional comp alg sort.
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

#This code will add time. I can comment it out later,
#once it has executed on trifid properly.


#Run first time after simulation. This should be on the end of the simulation code.

sortedlist=glob.glob('*soma*.dat_')#key=int
#sortedlist=natsort.natsorted(sortedlist, key=lambda y: y.lower())
sortedlist=natural_sort(sortedlist)
#sortedlist.sort()
cnt=0
if RANK==0:
   for fname in sortedlist:#sorted(glob.glob('*soma*.dat_')):
      if cnt<1001:
         z = np.loadtxt(fname)
         allmemb.append(z)
         print fname
         cnt+=1
   pickle.dump(allmemb, open("allmemb.p", "wb" ) )
   
COMM.Barrier()#Make sure all hosts are finished writing, before starting reading.

allmemb=pickle.load(open('allmemb.p','rb'))
allmemb#key=int
#allmemb=sorted(allmemb)
x=np.loadtxt('_spt.dat')
'''
tvec=[]
idvec=[]


for i in xrange(0,int(np.max(x[:,1]))):
   print int(x[i,1])
   indexs=np.where(x[:,1]==i)
   print x[indexs,0]
   print 'times'
'''
h('objref cell')
h('print cell')

h('{load_file("stdgui.hoc")}')
h('{load_file("nrngui.hoc")}')

h('load_file("import3d.hoc")')

h.xopen("morph4.hoc")
os.chdir(os.getcwd()+'/main')

acelllist=[]
indexs=[]
times=[]
#global NCELL
#NCELL= int(np.shape(allmemb)[0])
NCELL=1000
gidvec=[]
h('objref nc')
#
# Must be stored and saved in a file.
#
global tstop
tstop=1000
#
h('objref cells')
h('cells = new List()')
for i in range(RANK, NCELL-1, SIZE-1): #20 was int(len(allrows))
    #'s' needs to be a private variable in this method because later on 
    #it gets used differently in a different method, and only the s from 
    #here is understood because of the conflict.
    

    #This will be the reason for the skipping pattern. The simulation was skipping every SIZEth cell. Because it was incrementing by SIZE, when the true number was SIZE-1
    s=allmemb[i]    
    gidvec.append(i)
    indexs=np.where(x[:,1]==i)
    #print x[indexs,0]
    times=x[indexs,0]
    acelllist.append(acell(i,s,times))
    print getattr(acelllist[acell.cellcnt-1],'gid')
    print getattr(acelllist[acell.cellcnt-1],'membv')
    print getattr(acelllist[acell.cellcnt-1],'spike')

    
    h.cell=h.mkcell('cell-93-trace.CNG.swc')

    h('cell.gid1=py.i')
    h('cell.recvec=cell.recvec.from_python(py.s)')
    h('cell.spk_train.from_python(py.times)')
	
    h('pc.set_gid2node(py.i, pc.id)')# // associate gid i with this host
    h('cell.soma[0] nc =  new NetCon(&v(0.5), nil)') 
    h('pc.cell(py.i, nc)')#//                   

    #This is very important for the persistance of individual HOC cell objects.
    #Its not that I actually need a list of HOC, cell Objects I just need to keep all the references to old cell objects, such that these references don't get replaced/overwritten by new HOC cell objects
    h('cells.append(cell)')


    #for j in xrange(0,np.shape(x)[0]):
    #    if int(x[j,1])==i:
    #        spikes=x[j,0]
    #        h('cell.spk_train=cell.spk_train.append(py.spikes)')
    #        h('cell.spk_train.printf')

COMM.Barrier()
    
gidn=0
os.chdir('../')
h.xopen("nqs.hoc")

h('xopen("tools.hoc")')
COMM.Barrier()

#execfile('analysisp3.py')
import analysisp4 as ap4
trentm = np.zeros((NCELL, NCELL))
trentm=ap4.t_e(SIZE, NCELL, RANK, trentm,np)
# Do some plotting here!
#
trentm2 = np.zeros((NCELL, NCELL))
trentm2=ap4.t_e2(SIZE, NCELL, RANK, trentm2,np)
# Do some plotting here!
#


if RANK == 0:

    fig = plt.figure()
    fig.clf()
    im = plt.imshow(my_trentm, interpolation='nearest')
    plt.autoscale(True)
    plt.colorbar(im)
    plt.xlabel('columns = targets')
    plt.ylabel('rows = sources')
    plt.title('Transfer Entropy Matrix')
    plt.grid(True)
    sfin = str(SIZE) + str(NCELL) + 'Transfer_Entropy_matrix.png'
    fig.savefig(sfin)


    fig = plt.figure()
    fig.clf()
    im = plt.imshow(my_trentm2, interpolation='nearest')
    plt.autoscale(True)
    plt.colorbar(im)
    plt.xlabel('columns = targets')
    plt.ylabel('rows = sources')
    plt.title('Transfer Entropy Matrix')
    plt.grid(True)
    sfin = str(SIZE) + str(NCELL) + 'Transfer_Entropy_matrix.png'
    fig.savefig(sfin)
