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

#Note needs dat files as inputs.
sortedlist=glob.glob('*soma*.dat_')#key=int
#sortedlist=natsort.natsorted(sortedlist, key=lambda y: y.lower())
sortedlist=natural_sort(sortedlist)
#sortedlist.sort()

samp=len(sortedlist)-1
cnt=0

'''
for fname in sortedlist:#sorted(glob.glob('*soma*.dat_')):
   if cnt<samp:
      z = np.loadtxt(fname)
      allmemb.append(z)
      print fname
      cnt+=1
pickle.dump(allmemb, open("allmemb.p", "wb" ) )
'''
x=np.loadtxt('_spt.dat')


   

np.array(allmemb)
h.xopen("nqs.hoc")
h.xopen('analysis2.hoc')

h('objref storval')
h('objref vec1, vec2, t1, t2')
h('vec1=new Vector()')
h('vec2=new Vector()')

indexs0=[]
times0=[]
indexs1=[]
times1=[]
h('tstop1=0.0')
h.tstop1=np.max(x)
h('binsz=10')
     
def loopplot2(starti,stopi,startj, stopj,x):
   
   resultm=np.zeros((len(allmemb),len(allmemb)))
   resultm2=np.zeros((len(allmemb),len(allmemb)))
   i=0
   j=0
   for i in xrange(starti,stopi):
      for j in xrange(startj,stopj):
         #print np.shape(x)
         #print np.shape(indexs1)
         #print np.shape(times0)
         #print np.shape(indexs0)
         
         indexs0=np.where(x[:,1]==i)
         #if np.sum(indexs0)!=0:
         times0=x[indexs0,0][0]
         indexs1=np.where(x[:,1]==j)
         #  if np.sum(indexs1)!=0:
         times1=x[indexs1,0][0]
         print np.shape(times0)
         print np.shape(times1)
         print type(times1)
         h('objref vec0, vec1')
         h('vec0=new Vector()')
         h('vec1=new Vector()')

         h('vec0=vec0.from_python(py.times0)')
         h.vec0=h.vec0.from_python(times0)

         h('vec0.printf')
         #h.vec0.hist(h.vec0,0,(1000+10-1)/10,10)	  

         h('vec1=vec1.from_python(py.times1)')
         h.vec1=h.vec1.from_python(times1)

         #h.vec1.hist(h.vec1,0,(1000+10-1)/10,10)	  
         #h.vec0.printf
         #h.vec1.printf
        
         h('storval=normte(vec1,vec2,20)')
      
         h.storval.printf 
         #resultm[i][j]=h.storval.x[2]
         #print i, j, resultm[i][j]
         #print np.sum(resultm)
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
            
loopplot2(0,249,0,249,x)


def loopplot(starti,stopi,startj, stopj):
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
         #import rpy2.robjects.numpy2ri
      

         h('vec1.resample(vec1,200/40000)')
         h('vec2.resample(vec2,200/40000)')
         h('vec1.scale(0,130)')
         h('vec2.scale(0,130)')
         h('storval=normte(vec1,vec2,20)')
      
         h.storval.printf 
         resultm[i][j]=h.storval.x[2]
         print i, j, resultm[i][j]
         print np.sum(resultm)
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
import sys
pf=sys.argv


if(pf==0):
   loopplot(0,249,0,249)
if(pf==1):
   loopplot(250,499,250,499)
if(pf==2):
   loopplot(500,749,500,749)
if(pf==3):
   loopplot(750,999,750,999)
if(pf==4):
   loopplot(1000,1249,1000,1249)
if(pf==5):
   loopplot(1250,1500,1250,1500)
if(pf==6):
   loopplot(1500,1749,1500,1749)
if(pf==7):
   loopplot(1749,2000,1749,2000)

if(pf==8):
   loopplot(0,249,249,500)
if(pf==9):
   loopplot(249,500,500,749)
if(pf==10):
   loopplot(500,749,749,1000)
if(pf==11):
   loopplot(749,1000,1000,1249)
if(pf==12):
   loopplot(1000,1249,1249,1500)
if(pf==13):
   loopplot(1000,1249,1249,1500)
if(pf==14):
   loopplot(1249,1500,1500,1749)
if(pf==15):
   loopplot(1500,1749,1749,2000)



#loopplot(0,500,500,1000)
#loopplot(500,1000,1000,1500)
#loopplot(1000,1500,1500,2000)
#loopplot(1500,2000,0,500)
if(pf==16):
   looplot(0,249,500,749)
if(pf==17):
   looplot(249,500,749,1249)
if(pf==18):
   looplot(500,749,1000,1249)
if(pf==19):
   looplot(749,1000,1249,1500)
if(pf==20):
   looplot(1000,1249,1500,1749)
if(pf==21):
   looplot(1249,1500,1749,2000)
if(pf==22):
   looplot(1500,1749,0,249)
if(pf==23):
   looplot(1749,2000,249,500)


if(pf==24):
   loopplot(0,249,749,1000)
if(pf==25):
   loopplot(249,500,1000,1249)
if(pf==26):
   loopplot(500,749,1249,1500)
if(pf==27):
   loopplot(749,1000,1500,1749)
if(pf==28):
   loopplot(1000,1249,1749,2000)
if(pf==29):
   loopplot(1249,1500,0,249)
if(pf==30):
   loopplot(1500,1749,249,500)
if(pf==31):
   loopplot(1749,2000,500,749)

#loopplot(0,500,1000,1500)
#loopplot(500,1000,1500,2000)
#loopplot(1000,1500,0,500)
#loopplot(1500,2000,500,1000)

if(pf==32):
   loopplot(0,249,1500,1749)
if(pf==33):
   loopplot(249,500,1749,2000)
if(pf==34):
   loopplot(500,749,0,249)
if(pf==35):
   loopplot(749,1000,249,500)
if(pf==36):
   loopplot(1000,1249,500,749)
if(pf==37):
   loopplot(1249,1500,749,1000)
if(pf==38):
   loopplot(1500,1749,1000,1249)
if(pf==39):
   loopplot(1749,2000,1249,1500)


#loopplot(0,500,1500,2000)
#loopplot(500,1000,0,500)
#loopplot(1000,1500,500,1000)
#loopplot(1500,2000,1000,1500)







###FLAG 2
'''
loopplot(0,500,0,500)
loopplot(500,1000,500,1000)
loopplot(1000,1500,1000,1500)
loopplot(1500,2000,1500,2000)

loopplot(0,500,500,1000)
loopplot(500,1000,1000,1500)
loopplot(1000,1500,1500,2000)
loopplot(1500,2000,0,500)

loopplot(0,500,1000,1500)
loopplot(500,1000,1500,2000)
loopplot(1000,1500,0,500)
loopplot(1500,2000,500,1000)

loopplot(0,500,1000,1500)
loopplot(500,1000,1500,2000)
loopplot(1000,1500,0,500)
loopplot(1500,2000,500,1000)

loopplot(0,500,1500,2000)
loopplot(500,1000,0,500)
loopplot(1000,1500,500,1000)
loopplot(1500,2000,1000,1500)
'''


#loopplot(0,750,0,750)#len(allmemb),0,len(allmemb),0)
#loopplot(750,1500,750,1500)#len(allmemb),0,len(allmemb),0)

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
