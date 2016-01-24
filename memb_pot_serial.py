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

samp=len(sortedlist)-1
cnt=0
for fname in sortedlist:#sorted(glob.glob('*soma*.dat_')):
   if cnt<samp:
      z = np.loadtxt(fname)
      allmemb.append(z)
      print fname
      cnt+=1
pickle.dump(allmemb, open("allmemb.p", "wb" ) )
   

np.array(allmemb)
h.xopen("nqs.hoc")
h.xopen('analysis2.hoc')

h('objref storval')
h('objref vec1, vec2, t1, t2')
h('vec1=new Vector()')
h('vec2=new Vector()')




def loopplot(starti,stopi,startj, stopj,flag):
   resultm=np.zeros((len(allmemb),len(allmemb)))
   resultm2=np.zeros((len(allmemb),len(allmemb)))
   i=0
   j=0
   for i in xrange(starti,stopi):
      for j in xrange(startj,stopj):
         h('objref vec1, vec2, vec1c, vec2c')
         h('vec1=new Vector()')
         h('vec2=new Vector()')
         h('vec1c=new Vector()')
         h('vec2c=new Vector()')
         h('t1=new Vector()')
         h('t2=new Vector()')


         h('vec1=vec1.from_python((py.allmemb[i]))')
         h('vec2=vec2.from_python((py.allmemb[j]))')
         h('vec1.resample(vec1,200/40000)')
         h('vec2.resample(vec2,200/40000)')
         h('vec1.scale(0,130)')
         h('vec2.scale(0,130)')
         h('storval=normte(vec1,vec2,20)')
      
         h.storval.printf 
         resultm[i][j]=h.storval.x[2]
         print i, j, resultm[i][j]
         print np.sum(resultm)
         '''
         if flag==1:
            h('vec1c.copy(vec1)')
            h('vec2c.copy(vec2)')
            h('vec2c=vec2c.spikebin(vec2c,0)')
            h('vec1c=vec1c.spikebin(vec1c,0)')
            tt = np.linspace(0.0, 0.025 * 10000, 10000)
            t1=[]
            t2=[]
            t2=h.vec2c.to_python()
            t1=h.vec1c.to_python()
            h('py.t2=vec2c.to_python()')
            h('py.t1=vec1c.to_python()')
            #print t1
            #print t2
            tspike1=[]
            tspike1=tt[np.nonzero(t1)]
            print tspike1
            tspike2=[]
            tspike2=tt[np.nonzero(t2)]
            print tspike2
            h.t2=h.t2.from_python(tspike2)
            h.t1=h.t1.from_python(tspike1)
         
            #h('t2.printf')
            #h('t1.printf')
         
            h('storval=normte(t1,t2,20)')
            print np.shape(resultm2), i, j
            print np.shape(resultm), i, j
            resultm2[i][j]=h.storval.x[2]
            print i, j, resultm[i][j]
            print np.sum(resultm)
         if flag==0:
         '''
   '''
   np.save(str(starti)+str(stopi)+str(startj)+str(stopj)+'resultm',resultm)
 
   fig = plt.figure()
   fig.clf()
   im = plt.imshow(resultm, interpolation='nearest')
   plt.autoscale(True)
   plt.colorbar(im)
   plt.xlabel('columns = targets')
   plt.ylabel('rows = sources')
   plt.title('Transfer Entropy Matrix')
   plt.grid(True)
   sfin = str(starti)+str(stopi)+str(startj)+str(stopj)+'Transfer_Entropy_matrix2.png'   
   fig.savefig(sfin)
   '''
   fig = plt.figure()
  
   fig.clf()
   im = plt.imshow(resultm, interpolation='nearest')
   plt.autoscale(True)
   plt.colorbar(im)
   plt.xlabel('columns = targets')
   plt.ylabel('rows = sources')
   plt.title('Transfer Entropy Matrix')
   plt.grid(True)
   sfin = str(starti)+str(stopi)+str(startj)+str(stopj)+str(flag)+'Transfer_Entropy_matrix2.png'
   fig.savefig(sfin)
   
   np.save(str(starti)+str(stopi)+str(startj)+str(stopj)+str(flag)+'resultm2',resultm)



#loopplot(starti,stopi,startj, stopj)



loopplot(0,len(allmemb),0,len(allmemb),0)

#for j in xrange (0,1500,10): #make two lists of numbers.
#   for i in xrange (500,2000,10):
#      loopplot(i,j,i,j,0)

# Takes int arguments not list. 
# 
# loopplot(range(0,1500,150),range(500,2000,150),range(0,1500,150),range(500,2000,150),0)
'''
for j in xrange (0,1500,10): #make two lists of numbers.
   for i in xrange (500,2000,10):
      loopplot(i,j,i+500,j+500,0)


for j in xrange (0,1500,10): #make two lists of numbers.
   for i in xrange (500,2000,10):
      loopplot(i,j,i+1000,j+1000,0)
for j in xrange (0,1500,10): #make two lists of numbers.
   for i in xrange (500,2000,10):
      loopplot(i,j,i+1500,j+1500,0)
'''
'''
loopplot(0,500,0,500,0)
loopplot(500,1000,500,1000,0)
loopplot(1000,1500,1000,1500,0)
loopplot(1500,2000,1500,2000,0)

loopplot(0,500,500,1000,0)
loopplot(500,1000,1000,1500,0)
loopplot(1000,1500,1500,2000,0)
loopplot(1500,2000,0,500,0)

loopplot(0,500,1000,1500,0)
loopplot(500,1000,1500,2000,0)
loopplot(1000,1500,0,500,0)
loopplot(1500,2000,500,1000,0)

loopplot(0,500,1500,2000,0)
loopplot(500,1000,0,500,0)
loopplot(1000,1500,500,1000,0)
loopplot(1500,2000,1000,1500,0)

###FLAG 2
loopplot(0,500,0,500,1)
loopplot(500,1000,500,1000,1)
loopplot(1000,1500,1000,1500,1)
loopplot(1500,2000,1500,2000,1)

loopplot(0,500,500,1000,1)
loopplot(500,1000,1000,1500,1)
loopplot(1000,1500,1500,2000,1)
loopplot(1500,2000,0,500,1)

loopplot(0,500,1000,1500,1)
loopplot(500,1000,1500,2000,1)
loopplot(1000,1500,0,500,1)
loopplot(1500,2000,500,1000,1)

loopplot(0,500,1000,1500,1)
loopplot(500,1000,1500,2000,1)
loopplot(1000,1500,0,500,1)
loopplot(1500,2000,500,1000,1)

loopplot(0,500,1500,2000,1)
loopplot(500,1000,0,500,1)
loopplot(1000,1500,500,1000,1)
loopplot(1500,2000,1000,1500,1)
'''



'''
i=0,500  j=0,500
i=500,1000  j=500,1000
i=1000,1500  j=1000,1500
i=1500,2000  j=1500,2000

'''
