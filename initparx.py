
import pdb
import glob
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from matplotlib.collections import PolyCollection, LineCollection
import os
import urllib2
import zipfile
from mpi4py import MPI

#import random1 as r
#rand=r.random1


from neuron import h
import brain_functions as bf



COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

pc = h.ParallelContext()
h('objref pc')
h.pc = pc


s = "mpi4py thinks I am %d of %d,\
 NEURON thinks I am %d of %d\n"
cw = MPI.COMM_WORLD
print s % (cw.rank, cw.size, pc.id(), pc.nhost())
h('time_start=pc.time()')

NFILE = 3175# Maximum on this machine.
#Number of files to consider. This always ends up being a number more than the number of cells used. However NCELL and NFILE are of a similar magnitude. Some cell files have to filtered out, and cannot be used.


h('iwp=1')
h('ew=0.0185')
h('delay=3000')
h('iw=1000 // divisor/ or scaler, so if its bigger weight should be smaller. Was 1000 but this should be stronger. delay equals internal weight.'
  )

h('dc=-1')  # delete cell. This retains the cell, but makes all of its section diameter 0us, and it makes all projection weights to it 0.
h('default_channels=0')
h('plastic=1')
h('get_dist=0')
h('dlc=0')  # dilate centres.
structural_plotting = 1

# I normaly use alphabet size of two symbols (0,1). There is no reason why this could not be increased. Using multilevel quantisation.

prunenet = 0
h('qd=0')  # quick and dirty for debugging
h('lpy=1')  # declare lpy here, the value is not important, its only important that the declaration exists.

# h('if(qd==1){prunenet=900}')

h('objref py')
h('py = new PythonObject()')

###

tstop = 1250  # 1500#550
h('tstop1=1250')

f = open('tstopfile', 'w')
f.write(str(tstop))
f.close()
###

h('minr=10')  # minumum r large, for some reason which I do not yet
h('minrs=10')  # minumum r small, understand
h('delay=3000')
h('iw=1000 // divisor/ or scaler, so if its bigger weight should be smaller. Was 1000 but this should be stronger. delay equals internal weight.'
  )
h('ew=0.0185// External weight.This value very appropriate. But perhaps smaller how to automate the tuning of synaptic weights?'
  )

h('{load_file("stdgui.hoc")}')
h('{load_file("nrngui.hoc")}')

h('load_file("import3d.hoc")')
h('strdef workingdir')
h('workingdir = getcwd()')
h('strdef neurons_home')
h('neurons_home = neuronhome()')
h('strdef machinename')
h('machine_name(machinename)')
h('strdef worm')
h('strdef morphrat')
h('strdef morphhuman')
h('strdef graph_dir')
h('strdef results_dir')
h('strdef matrix_dir')
h('sprint(matrix_dir,"%s%s",workingdir,"matrix_dir")')

h('sprint(results_dir,"%s%s",workingdir,"results")')
h('sprint(graph_dir,"%s%s",workingdir,"graphs")')
h('sprint(morphrat,"%s%s",workingdir,"main")')
h('sprint(morphhuman,"%s%s",workingdir,"SWC-2013human")')

# h('sprint(worm,"%s%s",workingdir,"/openworm/CNG version/")')

checkpoint_interval = 50000.

allrows = pickle.load(open('allrows.p', 'rb'))
h.xopen('morph4.hoc')
h.xopen('nqs.hoc')

os.system('pwd')
import re

# os.chdir('/home/zaza3/Downloads/trunk/examples/expericomp20140421/')

cnt1 = pc.id

# I do not know how to refer to relative paths in Python,
# the below emulates a call to a relative path.

os.chdir(os.getcwd() + '/main')


allrows2 = bf.prep_list(NFILE,allrows)

NCELL=int(np.shape(allrows2)[0])


cnt = 0
i = 0
gidvec = []
h('objref nc')
h('objref cell')
h('objref cells')
h('cells = new List()')
h('objref gidvec')
h('gidvec =new Vector()')
h('objref inter_list')
h('inter_list=new List()')
h('objref pyr_list')
h('pyr_list = new List()')

h('objref aspiny_list, hipp, neoc, basalf, basalg')
h('aspiny_list=new List()')
h('hipp=new List()')
h('neoc=new List()')
h('basalf=new List()')
h('basalg=new List()')




ie = np.zeros((2, NCELL / SIZE + 1))

# I think its because of the previous manipulation of NCELL
iesize = NCELL# A bit confused about why this works, and not the former
ie = np.zeros((2, NCELL / SIZE + 1))


ie0 = np.zeros((iesize, 2))
ie1 = np.zeros((iesize, 2))


for (i, row) in enumerate(allrows2):
    if 'interneuron' in row:
        print row
    if 'pyramid' in row:
        print row




#h.cells=bf.mb(RANK, NCELL, SIZE, allrows2, gidvec,h)
h.cells, allrows,ie0,ie1=bf.mb(RANK, NCELL, SIZE, allrows2, gidvec,h,1150,1250)

print np.shape(allrows2)[0], NCELL, h('cells.count')
h('forall{ for(x,0){ insert xtra }}')
h('forall{ for(x,0){ insert extracellular}}')

h('chdir("../")')
h('xopen("interpxyz.hoc")')
h('grindaway()')

# h('forall{ for(x,0){ print x_xtra }}')

h('system("pwd")')
h('xopen("seclists.hoc")')

h('objref coords')
h('coords = new Vector(5)')

h('objref coords2')
h('coords2 = new Vector(3)')

# print coords, RANK

h('objref syn_')
h('objref synlist')
h('synlist= new List()')
h('objref nc')
h('objref nclist')
h('nclist=new List()')
COMM.Barrier()  # wait until all hosts get to this point

h('objref strobj')
h('strobj = new StringFunctions()')
h('strdef sec_string2')
h('strdef cell_here')
h('strdef tail_target')
h('objref source, target')


icm = np.zeros((NCELL, NCELL))
ecm = np.zeros((NCELL, NCELL))

i = 0
j = 0
cnt = 0
cnt_sec = 0

for cell in h.cells:
    h('print cell.gid1')


h.nclist,ecm,icm=bf.wirecells(RANK,NCELL,SIZE,h,icm,ecm)

#h.nclist,ie0,ie1=ip.wirecells(RANK,NCELL,SIZE,h,icm,ecm)


#wirecells(RANK,NCELL,SIZE,h)
