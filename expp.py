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
#GPIO pins on microcontoller are both TX and RX. As in A CPU can trivially transmit and recieve.
#import neuron
from neuron import h
import LFPy#Importing LFPy imports neuron. NEURON has to be imported after mpi4py


#initialize the MPI interface
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()



pc = h.ParallelContext()
h('objref pc')
h.pc=pc
#h('pc = new ParallelContext()')
#h('pc = py.pc')
s = "mpi4py thinks I am %d of %d,\
 NEURON thinks I am %d of %d\n"
cw = MPI.COMM_WORLD
print s % (cw.rank, cw.size, \
           pc.id(),pc.nhost())

#h.dt = 0.05
#tstop = 1025
#h('print dt')

#h('pc.psolve(1000)')# Interestingly this does cause an error its just non fatal. Raising the question. Is the installation broken? What parallel python NEURON simulation would run at this stage?

#

h('dc=-1')  # delete cell. This retains the cell, but makes all of its section diameter 0us, and it makes all projection weights to it 0.
h('default_channels=0')
h('plastic=1')
h('get_dist=0')
h('dlc=0')  # dilate centres.
structural_plotting=1
h('n_cell=1')
n_cell=h.n_cell
h('numcell=n_cell')  # should really use sed. Such as to reduce duplicate variable to one.
n_cell=h.n_cell
NCELL=60

#I normaly use alphabet size of two symbols (0,1). There is no reason why this could not be increased. Using multilevel quantisation.
prunenet=0
h('qd=0')#quick and dirty for debugging
h('lpy=1')#declare lpy here, the value is not important, its only important that the declaration exists.
h('if(qd==1){prunenet=900}')
h('objref py')
#h('objref py')
#h('py = new PythonObject()')
#h('if(qd==1){ py.structural_plotting=0')
h('minr=10')#minumum r large, for some reason which I do not yet 
h('minrs=10')#minumum r small, understand
#the number of connections made, depends somewhat on the last run, this may be 
# caused by morph3.hoc, deleting inappropriate cells occasionally.
# h("offset=0")
##OLD
h('offset=0')
##
h('parallels=0')
# h("tstop=1000//1 second =1000 ms")
h('delay=3000')
h('iw=1000 // divisor/ or scaler, so if its bigger weight should be smaller. Was 1000 but this should be stronger. delay equals internal weight.'
  )
h('iwp=1')
h('ew=0.0185// External weight.This value very appropriate. But perhaps smaller how to automate the tuning of synaptic weights?'
  )
h('run_iter=1')
h('msdp=0 //basically prune some connections distal from the soma.')
h('fast_wire=0')
h('ff=0')
h('if(ff==1){ get_dist=0 }')
h('if(ff==1){ n_cell=40 }')
h('if(ff==1){ prunenet=0 }')
h('no_b=0')
h('if (numcell>99){large=1}')
h('if (numcell>99){large_scale=1}')
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
h('sprint(worm,"%s%s",workingdir,"/openworm/CNG version/")')

checkpoint_interval = 50000.

def prun(tstop):
  cvode = h.CVode()
  cvode.cache_efficient(1)
  #pc.spike_compress(0,0,1)
  pc.setup_transfer()
  mindelay = pc.set_maxstep(10)
  if RANK == 0: print 'mindelay = %g'%mindelay
  runtime = h.startsw()
  exchtime = pc.wait_time()

  inittime = h.startsw()
  h.stdinit()
  inittime = h.startsw() - inittime
  if RANK == 0: print 'init time = %g'%inittime
  
  while h.t < tstop:
    told = h.t
    tnext = h.t + checkpoint_interval
    if tnext > tstop:
      tnext = tstop
    pc.psolve(tnext)
    if h.t == told:
      if RANK == 0:
        print "psolve did not advance time from t=%.20g to tnext=%.20g\n"%(h.t, tnext)
      break

  runtime = h.startsw() - runtime
  comptime = pc.step_time()
  splittime = pc.vtransfer_time(1)
  gaptime = pc.vtransfer_time()
  exchtime = pc.wait_time() - exchtime
  if RANK == 0: print 'runtime = %g'% runtime
  print comptime, exchtime, splittime, gaptime
  #printperf([comptime, exchtime, splittime, gaptime])/

def stationary_poisson(nsyn,lambd,tstart,tstop):
    ''' Generates nsyn stationary possion processes with rate lambda between tstart and tstop'''
    interval_s = (tstop-tstart)*.001
    spiketimes = []
    for i in xrange(nsyn):
        spikecount = np.random.poisson(interval_s*lambd)
        spikevec = np.empty(spikecount)
        if spikecount==0:
            spiketimes.append(spikevec)
        else:
            spikevec = tstart + (tstop-tstart)*np.random.random(spikecount)
            spiketimes.append(np.sort(spikevec)) #sort them too!

    return spiketimes


#Fetch Mainen&Sejnowski 1996 model files
if not os.path.isfile('/home/zaza3/trunk/examples/patdemo/cells/j4a.hoc') and RANK==0:
    #get the model files:
    u = urllib2.urlopen('http://senselab.med.yale.edu/ModelDB/eavBinDown.asp?o=2488&a=23&mime=application/zip')
    localFile = open('patdemo.zip', 'w')
    localFile.write(u.read())
    localFile.close()
    #unzip:
    myzip = zipfile.ZipFile('patdemo.zip', 'r')
    myzip.extractall('.')
    myzip.close()

#resync MPI threads
COMM.Barrier()

# Define cell parameters

# Define synapse parameters
synapse_parameters = {
    'idx' : 0, # to be set later
    'e' : 0.,                   # reversal potential
    'syntype' : 'AMPA',       # synapse type
    'tau' : 5.,                 # syn. time constant
    'weight' : .001,            # syn. weight
    'record_current' : True,
}

# Define electrode parameters
point_electrode_parameters = {
    'sigma' : 0.3,      # extracellular conductivity
    'x' : 0.,  # electrode requires 1d vector of positions
    'y' : 0.,
    'z' : 0.,
}


#number of units



n_cells = int(n_cell)#SIZE
cell_id = RANK
#execfile('/home/zaza3/Downloads/trunk/examples/expericomp20140421/atlas_reader3.py')

allrows=pickle.load(open('allrows.p','rb'))
h.xopen("/home/zaza3/trunk/examples/expericomp20140421/morph4.hoc")
h.xopen("/home/zaza3/trunk/examples/expericomp20140421/nqs.hoc")

os.system('pwd')
import re


os.chdir('/home/zaza3/Downloads/trunk/examples/expericomp20140421/')
#h.xopen("find_distances_.hoc")

#set the numpy random seeds
global_seed = 1234
np.random.seed(global_seed)


#assign cell positions
x_cell_pos = np.linspace(-250., 250., n_cells)

z_rotation = np.random.permutation(np.arange(0., np.pi, np.pi / n_cells))

#synaptic spike times
n_pre_syn = 1000
pre_syn_sptimes = stationary_poisson(nsyn=n_pre_syn, lambd=5., tstart=0, tstop=1000)

#re-seed the random number generator
cell_seed = global_seed + cell_id
np.random.seed(cell_seed)
cnt1=pc.id
#fn = re.compile('*.swc')
h('load_file("pyramidal_cell_14Vb.hoc")')
#h('objref cell')
#h('cell = new PyramidalCell()')
#h('print cell')
#quit()
os.chdir('/home/zaza3/trunk/examples/expericomp20140421/main/')
h('objref cell')
h('objref cells')
h('cells = new List()')
h('objref gidvec')
h('gidvec =new Vector()')


#cells=0
#cells=[cells]*40
#cells=[ LFPy.Cell(**cell_parameters) for i in range(int(n_cell+1))];
cells=[ 0 for i in range(int(n_cell+1))];


color_vec = [plt.cm.rainbow(int(x*256./n_cells)) for x in xrange(int(len(allrows)))]

#cells=[for x in xrange(0,50)]
#This is to make every cell a different color.

#prun(0.1)

h('objref py')
h('py = new PythonObject()')
cnt=0 #This value cnt, will be different on different hosts and thats what is required.
#prun(0.1)

#def make_cells(NCELL):
allrows2=[]
i=0
for i in range(RANK, int(len(allrows)), SIZE): #20 was int(len(allrows))
    #if((cnt<n_cell) and (0<i)):#While there are less than 100 cells on every node.
    s=allrows[i]    

    if(cnt>1):
        if(int(len(s))>9):#//This condition is counter intuitive many cells 
            storename= str(s[3]) #//simply being inside a loop, may be the main pro
            #print storename, s[3]
            allrows2.append(allrows[i])
    cnt+=1
print np.shape(allrows2), np.shape(allrows)

cnt=0
i=0
gidvec=[]
h('objref nc')

for i in range(RANK, NCELL, SIZE): #20 was int(len(allrows))
    s=allrows2[i]    
    storename= str(s[3]) #//simply being inside a loop, may be the main pro
    if(re.search('.swc',storename)):                                    
        h.cell=h.mkcell(storename)
        #h('cell.init()')
        #h('cell = new PyramidalCell()')
        #h('cell.synapses()')
        
        h('cell.geom_nseg()')
        h('cell.pyr()')
        h('cell.gid1=py.i')
        h('cell.gvpre.printf')
        h('cell.gvpost.printf')
        h('cell.nametype=py.str(py.s[5])') 
        h('cell.num_type=py.int(py.s[6])')
        h('cell.population=py.str(py.s[7])')
        h('cell.reponame=py.str(py.storename)')


        h('strdef cellposition')
        h('sprint(cellposition,"%s%d%s%d%s%d%s","cell.position(",py.int(py.s[0]),",",py.int(py.s[1]),",",py.int(py.s[2]),")")')
        #h('execute(cellposition)')
	
        h('pc.set_gid2node(py.i, pc.id)')# // associate gid i with this host
        #h('nc = cell.connect2target(nil)')# // attach spike detector to cell
        h('cell.soma nc =  new NetCon(&v(0.5), nil)') 
        h('pc.cell(py.i, nc)')#//                   
	#cell=h.cell
        #pc.cell(i, h.NetCon(cell.soma(.5)._ref_v, None, sec=cell.soma))
 
        h('cells.append(cell)')
        h('gidvec.append(py.i)')
        gidvec.append(i)
        h('print py.i')
        cnt+=1
                
             
#h('pc.psolve(1000)')


#prun(0.1)
#quit()
print 'past... the test'
#h('forall{ for(x,0){ insert xtra }}') 
#h('forall{ for(x,0){ insert extracellular}}') 

h('chdir("../")')
print 'death by extracellular? will it pass the basic test'
#prun(10)
print 'past basic test'

#h('xopen("interpxyz.hoc")')
#h('grindaway()')
#h('forall{ for(x,0){ print x_xtra }}')
print 'will it pass the basic test'
prun(10)
print 'past basic test'

h('system("pwd")')
h('xopen("seclists.hoc")')

h('objref coords')
h('coords = new Vector(5)')


h('objref coords2')
h('coords2 = new Vector(3)')



#print coords, RANK

h('objref syn_')
h('objref synlist')
h('synlist= new List()')
h('objref nc')
h('objref nclist')
h('nclist=new List()')
COMM.Barrier()#wait until all hosts get to this point

#h('objref cell2')
#h('cell2=pc.gid2cell(gidvec.x[0])')
'''
1
'''

h('objref strobj')
h('strobj = new StringFunctions()')
h('strdef sec_string2')
h('strdef cell_here')
h('strdef tail_target')

h('objref source, target')

'''
h('objref synx')
if pc.gid_exists(0):
  h('cell=pc.gid2cell(0)')
  h('print cell')
  h('cell.soma[0] synx = new ExpSyn(0.5)')
  h('cell.synlist.append(synx)')
  h('print synx')
  h('print cell.synlist.o(0)')
  h('nc = pc.gid_connect (0, synx)')
  h('print nc')

  h('nc.delay = 100')

  h('nc.weight = 0.01')
  h('nclist.append (nc)')
  h('print nclist')
print 'will it pass the basic test'
prun(0.1)
print 'past basic test'
'''
'''
h('objref target, nc, syn')
for i in xrange(0, NCELL-1): 
  targid = (i+1)%NCELL

  if (not pc.gid_exists(targid)): continue
  h('target = pc.gid2cell(targid)')
  h('print target')
  h('syn = target.pre_list.object(0)')
  h('print syn')
  h('nc = pc.gid_connect(py.i, syn)')
  h('nclist.append (nc)')
  h('nc.delay = 1')
  h('nc.weight = 0.01')

'''
'''
h('objref synx')
if pc.gid_exists(0):
  cell=pc.gid2cell(0)
  #print 
  h('print cell')
  syn = h.ExpSyn(0.5,sec=cell.soma[0])
  h('help(syn)')
  #print h.syn
  nc = pc.gid_connect (10, syn)
  #print h.nc
  h('print nc')
  h('nclist.append (nc)')

  h('nc.delay = 10')

  h('nc.weight = 0.01')

  #nc = h.NetCon(pre(0.5)._ref_v, syn)
  #nc.weight[0] = 2.0
'''
h('hoc_stdout("debug.txt")')





i=0
j=0
cnt=0
cnt_sec=0
for s in xrange(0,SIZE):
    for j in xrange(0,NCELL):
        #for i, cells.count-1{
        #h.cells.o(j)
        if (RANK==s):
            coordlist=[[0 for x in xrange(0,2)] for x in xrange(0,2)]
     
            test=0
            print j, ' j ', s, ' s '
            h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
            if (test!=0):
                h('source=pc.gid2cell(py.int(py.j))')
                for sec in h.source.all:
                    for seg in sec:
                        '''
                        get_cox=str('coords.x[0]=x_xtra('+str(seg.x)+')')
                        h(get_cox)
                        h('print py.get_cox')
                        get_coy=str('coords.x[1]=y_xtra('+str(seg.x)+')')
                        h(get_coy)
                        get_coz=str('coords.x[2]=z_xtra('+str(seg.x)+')')
                        h(get_coz)
                        '''
                        h.coords.x[3]=int(j)
                        h('coords.x[4]=py.seg.x')#ie the gidvec.  
                        
                        coords=np.array(h.coords.to_python(),dtype=np.float32)
                    #print 's: ', s , 'coords', coords, '\n'
                        coordlist[0][1]=coords
                        #h('strdef secnames')
                        #h('cells.o(0).soma[0]{ secnames=secname() }')
                        secnames=sec.name() #h.secnames
                        coordlist[1][0]=str(secnames)
                        print coordlist
            if (test==0):
                coordlist=None
        else: 
            coordlist=None
        #Only gid 6 is only every being broadcast.
        data = COMM.bcast(coordlist, root=s)#ie root = rank

        if (data!=None):#and(RANK!=s):
            for i in xrange(0,NCELL):
                gidn=int(data[0][1][3])
                if(gidn==6):
                    print coordlist
                    print 'yes, not equal to gid 6'
                    break
                if(i!=int(gidn)):#if the gids are not the same
                    test2=0
                    print j, ' j ', i, ' i ', s, ' s '
                    h('py.test2=py.int(pc.gid_exists(py.int(py.i)))')
                    if (test2!=0):
                        h('target=pc.gid2cell(py.int(py.i))')
                        for sec in h.target.all:
                            for seg in sec:
            
            #for sec in h.spk_rx_ls:#iterate through post synaptic locations.
                #print "bcast finished and data on RANK %d is: "%RANK, data[0][1], data[1][0]#, secnames,' data ', data
                #for seg in sec:
                                '''    
                                print sec.name(), h(str('z_xtra(')+str(seg.x)+')'), 'rank= ', RANK, '\n'


                                h(str('coords2.x[2]=')+str('z_xtra(')+str(seg.x)+')')
                                h(str('coords2.x[1]=')+str('y_xtra(')+str(seg.x)+')')
                                h(str('coords2.x[0]=')+str('x_xtra(')+str(seg.x)+')')
                  

                                h('coordsx=0.0')
                                h.coordsx=data[0][1][0]
                                h('coordsy=0.0')
                                h.coordsy=data[0][1][1]
                                h('coordsz=0.0')
                                h.coordsz=data[0][1][2]
                                r=0.
                #print sec.name(), data[1][0][:]
                    
                #print 'code hangs here'
                                h('py.r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)')
                                r=float(r)
                                print ' r ', r
                                if(r<1):#I still need to debug the code, to find out why it won't search the finest levels of granularity.
                                    gidcompare=''
'''
                                secnames=sec.name()
                                cellind=int(secnames[secnames.find("Cell[")+5:secnames.find("].")])#This is the index of the post synaptic cell.
                        #h('py.pol=strcmp(Cell[py.int(py.cellin))]
                                post_syn=secnames+' '+'syn_ = new ExpSid('+str(seg.x)+')'
                                h('gidn=0')
                                h.gidn=int(data[0][1][3])
                        
                                h(post_syn)
                                h('syn_.cid=py.int(py.i)') #This makes the gid a property of the synapse. 
                                h('Cell[py.cellind].ampalist.append(syn_)')
                        
                        #cellind=(re.search(r'\[(.*)\]', secnames).group(1))
                                print cellind, 'cellind'
                #place the post synapse.
                                h('synlist.append (syn_)')
                                h('gidn=0')
                                h.gidn=int(data[0][1][3])
                                h('Cell[py.cellind].gvpre.append(gidn)')
                        
                        #print True(gidcompare==int(data[0][1][3])), ' is this the same cell?'
                                h('nc = pc.gid_connect(gidn, syn_)')
                        #connect the pre synapse to the post synapse.
                        #h('print coords.x[3], "this is the global identifier"')
                        #h('Cell[py.int(tail_target)].gvpre.append(gidn)')
                        #h('Cell[py.int(tail_target)].gvpost.append(py.int(py.gidcompare))')
                               
                                h('print nc')
                                print 'the gids of the connected cells are: ', i, ' ', gidn, '\n'
                                h('nclist.append(nc)')

                                h('nc.delay = (1/dt)+1000')
                                h('nc.weight = 1')
                                    

h('hoc_stdout("")')

COMM.Barrier()
h('objref tgt')
h('srcid =0')
h('system("rm net_schematic.txt")')
h('hoc_stdout("net_schematic.txt")')
COMM.Barrier() 
if RANK==0: print 'source\ttarget\tsynapse\n'
for i in xrange(0,SIZE):
   if i==RANK:
      nclist=0
      h('py.nclist=py.int(nclist.count-1)')
      for j in xrange(0,int(nclist)):
         
         h('srcid=nclist.o(py.j).srcgid()')
         h('tgt = nclist.o(py.j).syn')
         h('print srcid, " ", tgt.cid')
         #h('printf("%d\t%d\n", srcid, tgt.cid)')
h('hoc_stdout()')

prun(0.11)

#The fact that code executres above and not here, suggests that mindelay is still a problem
#h('pc.psolve(1000)')# This executes okay. Suggesting that the problem really is 
h('print nclist.count, " size nclist"')

#h('for i=0,cells.count-1{ cells.o(i).gvpre.printf } ') 
#h('for i=0,cells.count-1{ print cells.o(i).gid1 } ')

'''
COMM.Barrier() 
if RANK==0: print '\n gidvecs on the various hosts \n'
for i in xrange(0,SIZE): 
    if (i==RANK):
        print 'host', RANK, gidvec[:], '\n'
        h('gidvec.printf()')
    
COMM.Barrier() 



'''


COMM.Barrier() 
if RANK==0: print '\n cell numbers on the various hosts \n'
for i in xrange(0,SIZE): 
    if (i==RANK):
        print 'host', RANK, '\n'
        h('print "cells.count ", cells.count')
COMM.Barrier() 

'''
#h('xopen("rigp.hoc")')
tstop=1000

#h.pc.psolve(tstop)

h('print dt')
h('pc.set_maxstep(10)')
#h('load_file("nrngui.hoc")')#Already been called before above.
h('stdinit()')

h('pc.psolve(1000)')
'''
'''

def runp(tstop):
    lfp=1
    if lfp == 1:
        h.init()
    if lfp != 1:
    	h.finitialize()  # set initial variables.
    while h.t < tstop:  # Perform integration over the cable equation.
        if lfp == 1:
            h.advance()  # this method is designed to work with field and may take longer.  
            h.pc.psolve(tstop)
        if lfp != 1:
            h.fadvance()
            h.pc.psolve(tstop)



runp(1000)
'''
                
        #if(test!=0):

'''            
    for sec in h.spk_trig_ls:#iterate through pre synaptic locations.
        #if not(cnt<79):
        #if not(cnt<cnt_sec-1):
        #    print 'broke here'
        #    break
        #Its not a memory capacity problem that limits us to 80.
        if (RANK==s):
            coordlist=[[0 for x in xrange(0,2)] for x in xrange(0,2)]
            h('coords.x[0]=cells.o(0).soma[0].x_xtra()')
            h('coords.x[1]=cells.o(0).soma[0].y_xtra()')
            h('coords.x[2]=cells.o(0).soma[0].z_xtra()')
            h('coords.x[3]=py.int(gidvec.x[0])')
            coords=np.array(h.coords.to_python(),dtype=np.float32)
            #print 's: ', s , 'coords', coords, '\n'
            coordlist[0][1]=coords
            h('strdef secnames')
            h('cells.o(0).soma[0]{ secnames=secname() }')
            secnames=sec.name() #h.secnames
            coordlist[1][0]=str(secnames)
        else: 
            coordlist=None
        #broadcast post synaptic site is wrong. It should be presynaptic cell that is broadcast.
        #print 's RANK, SIZE coordlist, ',s,' ',' ', RANK,' ', SIZE,' ',coordlist, '\n'
        #if coordlist!=None:
        #COMM.Barrier()    
        
        cnt+=1
        #print cnt, ' ', bool(cnt<cnt_sec), sec.name(), ' ' , cnt_sec

            #if(cnt<80):
        
        data = COMM.bcast(coordlist, root=s)#ie root = rank
            #print ' rank ', RANK, ' data', data

        
        if (data!=None)and(RANK!=s):
            for sec in h.spk_rx_ls:#iterate through post synaptic locations.
                #print "bcast finished and data on RANK %d is: "%RANK, data[0][1], data[1][0]#, secnames,' data ', data
                h('coords2.x[0]=cells.o(0).soma[0].x_xtra()')
                h('coords2.x[1]=cells.o(0).soma[0].y_xtra()')
                h('coords2.x[2]=cells.o(0).soma[0].z_xtra()')
                h('coordsx=0.0')
                h.coordsx=data[0][1][0]
                h('coordsy=0.0')
                h.coordsy=data[0][1][1]
                h('coordsz=0.0')
                h.coordsz=data[0][1][2]
                r=0.
                print sec.name(), data[1][0][:]
                    
                #print 'code hangs here'
                h('py.r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)')
                r=float(r)
                if(r<55):
                    #print r, RANK, ' r, RANK '
                    #h('print r, "this is r"')
                #Assign a post synapse.
                #secnames=data[1][0][:]
                    secnames=sec.name()
                    post_syn=secnames+' '+'syn_ = new AMPA(0.5)'
                    h(post_syn)
                #place the post synapse.
                    h('synlist.append (syn_)')
                    h('gidn=0')
                    h.gidn=int(data[0][1][3])
                    h('nc = pc.gid_connect(gidn, syn_)')#connect the pre synapse to the post synapse.
                #h('print coords.x[3], "this is the global identifier"')
                    h('print nc')
                    h('nc.delay = 1')
                    h('nc.weight = 0.01')
                    #print 'code hangs here'
'''
#COMM.Barrier()    

'''2'''



'''
i=0
cnt=0
cnt_sec=0
for s in xrange(0, SIZE) : 
    #for sec in h.allsec():#iterate through pre synaptic locations.
    #   cnt_sec+=1
    for sec in h.spk_trig_ls:#iterate through pre synaptic locations.
        #if not(cnt<79):
        #if not(cnt<cnt_sec-1):
        #    print 'broke here'
        #    break
        #Its not a memory capacity problem that limits us to 80.
        if (RANK==s):
            coordlist=[[0 for x in xrange(0,2)] for x in xrange(0,2)]
            h('coords.x[0]=cells.o(0).soma[0].x_xtra()')
            h('coords.x[1]=cells.o(0).soma[0].y_xtra()')
            h('coords.x[2]=cells.o(0).soma[0].z_xtra()')
            h('coords.x[3]=py.int(gidvec.x[0])')
            coords=np.array(h.coords.to_python(),dtype=np.float32)
            #print 's: ', s , 'coords', coords, '\n'
            coordlist[0][1]=coords
            h('strdef secnames')
            h('cells.o(0).soma[0]{ secnames=secname() }')
            secnames=sec.name() #h.secnames
            coordlist[1][0]=str(secnames)
        else: 
            coordlist=None
        #broadcast post synaptic site is wrong. It should be presynaptic cell that is broadcast.
        #print 's RANK, SIZE coordlist, ',s,' ',' ', RANK,' ', SIZE,' ',coordlist, '\n'
        #if coordlist!=None:
        #COMM.Barrier()    
        
        cnt+=1
        #print cnt, ' ', bool(cnt<cnt_sec), sec.name(), ' ' , cnt_sec

            #if(cnt<80):
        
        data = COMM.bcast(coordlist, root=s)#ie root = rank
            #print ' rank ', RANK, ' data', data

        
        if (data!=None)and(RANK!=s):
            for sec in h.spk_rx_ls:#iterate through post synaptic locations.
                #print "bcast finished and data on RANK %d is: "%RANK, data[0][1], data[1][0]#, secnames,' data ', data
                h('coords2.x[0]=cells.o(0).soma[0].x_xtra()')
                h('coords2.x[1]=cells.o(0).soma[0].y_xtra()')
                h('coords2.x[2]=cells.o(0).soma[0].z_xtra()')
                h('coordsx=0.0')
                h.coordsx=data[0][1][0]
                h('coordsy=0.0')
                h.coordsy=data[0][1][1]
                h('coordsz=0.0')
                h.coordsz=data[0][1][2]
                r=0.
                print sec.name(), data[1][0][:]
                    
                #print 'code hangs here'
                h('py.r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)')
                r=float(r)
                if(r<55):
                    #print r, RANK, ' r, RANK '
                    #h('print r, "this is r"')
                #Assign a post synapse.
                #secnames=data[1][0][:]
                    secnames=sec.name()
                    post_syn=secnames+' '+'syn_ = new AMPA(0.5)'
                    h(post_syn)
                #place the post synapse.
                    h('synlist.append (syn_)')
                    h('gidn=0')
                    h.gidn=int(data[0][1][3])
                    h('nc = pc.gid_connect(gidn, syn_)')#connect the pre synapse to the post synapse.
                #h('print coords.x[3], "this is the global identifier"')
                    h('print nc')
                    h('nc.delay = 1')
                    h('nc.weight = 0.01')
                    #print 'code hangs here'

#COMM.Barrier()    

'''    
"""

        #for sec in h.allsec():
coordlist=[[0 for x in xrange(0,2)] for x in xrange(0,2)]
if RANK==1:
    h('coords.x[0]=cells.o(0).soma[0].x_xtra()')
    h('coords.x[1]=cells.o(0).soma[0].y_xtra()')
    h('coords.x[2]=cells.o(0).soma[0].z_xtra()')
            #h('coords.x[4]=x')
    h('coords.x[3]=gidvec.x[0]')

    coords=np.array(h.coords.to_python(),dtype=np.float32)
    coordlist[0][1]=coords
    h('strdef secnames')
    h('cells.o(0).soma[0]{ secnames=secname() }')
    secnames=h.secnames
    coordlist[1][0]=str(secnames)
            #coordlist[1][0]=str(sec.name())
            #comm.Send([coords, MPI.float], dest=j, tag=77)
    #print coordlist
else: 
    coords= None
    secnames= None
    coordlist=None
            #data = COMM.bcast(secnames, root=RANK)
            #COMM.Barrier()
            #COMM.Recv([coords, MPI.FLOAT], source=0, tag=77)
data = COMM.bcast(coordlist, root=1)

COMM.Barrier()    
        #print secnames, 'secnames', RANK, data
#print coordlist, ' coords ',' rank ', RANK
print ' rank ', RANK, ' data', data

        #The secnames of the post synapses, and coordinates have been broadcast.
        #Because its a broadcast. One statement sends and recieves. Its implied that all hosts will recieve.
#if coordlist!=None:
if (data!=None)and(RANK!=1):

        #if secnames!=None:
                #print "bcast finished and data on RANK %d is: "%RANK, coords, secnames,' data ', data
    h('coords2.x[0]=cells.o(0).soma[0].x_xtra()')
    h('coords2.x[1]=cells.o(0).soma[0].y_xtra()')
    h('coords2.x[2]=cells.o(0).soma[0].z_xtra()')
                #h('print coords2.x[2]')
    h('coordsx=0.0')
    coordlist=data

    h.coordsx=coordlist[0][1]
    h('coordsy=0.0')
    h.coordsy=coordlist[1][1]
    h('coordsz=0.0')
    h.coordsz=coordlist[2][1]
    r=0.
    h('r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)')
    r=float(h.r)
    print r
                #h('print r, "this is r"')
                #Assign a post synapse.
    post_syn=secnames+' '+'syn_ = new AMPA(0.5)'
    h(post_syn)
                #place the post synapse.
    h('synlist.append (syn_)')
    h('nc = pc.gid_connect(coords.x[3], syn_)')#connect the pre synapse to the post synapse.
                #h('print coords.x[3], "this is the global identifier"')
    h('nc.delay = 1')
    h('nc.weight = 0.01')
            
    #RANK+=1
#COMM.Barrier()
"""
        

        
