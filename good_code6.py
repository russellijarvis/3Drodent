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
import matplotlib
matplotlib.use('Agg')
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

h('iwp=1')
h('ew=0.0185')
h('delay=3000')
h('iw=1000 // divisor/ or scaler, so if its bigger weight should be smaller. Was 1000 but this should be stronger. delay equals internal weight.')


h('dc=-1')  # delete cell. This retains the cell, but makes all of its section diameter 0us, and it makes all projection weights to it 0.
h('default_channels=0')
h('plastic=1')
h('get_dist=0')
h('dlc=0')  # dilate centres.
structural_plotting=1

NCELL=35
tstop=50
h('tstop=0')
h.tstop=tstop
h('print tstop')
quit()
#I normaly use alphabet size of two symbols (0,1). There is no reason why this could not be increased. Using multilevel quantisation.
prunenet=0
h('qd=0')#quick and dirty for debugging
h('lpy=1')#declare lpy here, the value is not important, its only important that the declaration exists.
h('if(qd==1){prunenet=900}')
h('objref py')
h('minr=10')#minumum r large, for some reason which I do not yet 
h('minrs=10')#minumum r small, understand
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
#h('sprint(worm,"%s%s",workingdir,"/openworm/CNG version/")')

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
    if h.t%100==0:
      print 'working', h.t
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


allrows=pickle.load(open('allrows.p','rb'))
h.xopen("/home/zaza3/trunk/examples/expericomp20140421/morph4.hoc")
h.xopen("/home/zaza3/trunk/examples/expericomp20140421/nqs.hoc")

os.system('pwd')
import re


os.chdir('/home/zaza3/Downloads/trunk/examples/expericomp20140421/')

cnt1=pc.id
os.chdir('/home/zaza3/trunk/examples/expericomp20140421/main/')


h('objref py')
h('py = new PythonObject()')

def prep_list():
  cnt=0
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
  return allrows2
  print np.shape(allrows2), np.shape(allrows)


allrows2=prep_list()

#def mbs(allrows2): Need to declare global variables to make the below algorithm a definition.
cnt=0
i=0
gidvec=[]
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


for i in range(RANK, NCELL, SIZE): #20 was int(len(allrows))
  s=allrows2[i]    
  storename= str(s[3]) #//simply being inside a loop, may be the main pro
  if(re.search('.swc',storename)):                                    
    h.cell=h.mkcell(storename)
    #h('cell.init()')
    h('cell.geom_nseg()')
    #h('cell.pyr()')
    h('cell.gid1=py.i')
    h('cell.gvpre.printf')
    h('cell.gvpost.printf')
    h('cell.nametype=py.str(py.s[5])') 
    h('cell.num_type=py.int(py.s[6])')
    h('cell.population=py.str(py.s[7])')
    h('cell.reponame=py.str(py.storename)')
    h('cell.div.resize(py.int(py.NCELL))')
    h('cell.conv.resize(py.int(py.NCELL))')
    h('cell.nametype=py.str(py.s[5])') 
    h('cell.num_type=py.int(py.s[6])')
    h('cell.population=py.str(py.s[7])')
    h('cell.reponame=py.str(py.storename)')
  


    h('if(strcmp("pyramid",py.s[5])==0){pyr_list.append(cell)}')
    h('if(strcmp("pyramid",py.s[5])==0){ cell.pyr()}')

    h('if(strcmp("interneuron",py.s[5])==0){ inter_list.append(cell) }')
    h('if(strcmp("interneuron",py.s[5])==0){ cell.basket()}')
    h('if(strcmp("interneuron",py.s[5])==0){ print "inter neuron"}')
  
    h('if(strcmp("aspiny",py.s[5])==0){ aspiny_list.append(cell) }')  
    h('if(strcmp(cell.population,"hippocampus")==0){ hipp.append(cell)  }') 
    h('if(strcmp(cell.population,"neocortex")==0){ neoc.append(cell) }') 
    #h('if(strcmp(cell.population,"basalforebrain")==0){ basalf.append(cell) }') 
   # h('if(strcmp(cell.population,"basalganglia")==0){ basalg.append(cell) }')
    h('strdef cellposition')
    h('sprint(cellposition,"%s%d%s%d%s%d%s","cell.position(",py.int(py.s[0]),",",py.int(py.s[1]),",",py.int(py.s[2]),")")')
      #h('execute(cellposition)')

	
    h('pc.set_gid2node(py.i, pc.id)')# // associate gid i with this host
        #h('nc = cell.connect2target(nil)')# // attach spike detector to cell
    h('cell.soma[0] nc =  new NetCon(&v(0.5), nil)') 
    h('pc.cell(py.i, nc)')#//                   
	#cell=h.cell
        #pc.cell(i, h.NetCon(cell.soma(.5)._ref_v, None, sec=cell.soma))
 
    h('cells.append(cell)')
    h('gidvec.append(py.i)')
    gidvec.append(i)
    h('print py.i')
    cnt+=1

#mbs(allrows2)             
#h('pc.psolve(1000)')



h('forall{ for(x,0){ insert xtra }}') 
h('forall{ for(x,0){ insert extracellular}}') 

h('chdir("../")')
h('xopen("interpxyz.hoc")')
h('grindaway()')
#h('forall{ for(x,0){ print x_xtra }}')

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



#h('hoc_stdout("debug.txt")')





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
            #print j, ' j ', s, ' s '
            h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
            if (test!=0):
                h('source=pc.gid2cell(py.int(py.j))')
                for sec in h.source.spk_trig_ls:
                    for seg in sec:
                        get_cox=str('coords.x[0]=x_xtra('+str(seg.x)+')')
                        h(get_cox)
                        #h('print py.get_cox')
                        get_coy=str('coords.x[1]=y_xtra('+str(seg.x)+')')
                        h(get_coy)
                        get_coz=str('coords.x[2]=z_xtra('+str(seg.x)+')')
                        h(get_coz)
                        h.coords.x[3]=int(j)
                        h('coords.x[4]=py.seg.x')#ie the gidvec.  
 
                        coords=np.array(h.coords.to_python(),dtype=np.float32)
                    #print 's: ', s , 'coords', coords, '\n'
                        coordlist[0][1]=coords
                        #h('strdef secnames')
                        #h('cells.o(0).soma[0]{ secnames=secname() }')
                        secnames=sec.name() #h.secnames
                        coordlist[1][0]=str(secnames)
                        #print coordlist
            if (test==0):
                coordlist=None
        else: 
            coordlist=None
        #Only gid 6 is only every being broadcast.
        data = COMM.bcast(coordlist, root=s)#ie root = rank

        if (data!=None):#and(RANK!=s):
            for i in xrange(0,NCELL):
                gidn=int(data[0][1][3])
                #if(gidn==6):
                #    print coordlist
                #    print 'yes, not equal to gid 6'
                #    break
                if(i!=int(gidn)):#if the gids are not the same
                    test2=0
                    #print j, ' j ', i, ' i ', s, ' s '
                    h('py.test2=py.int(pc.gid_exists(py.int(py.i)))')
                    if (test2!=0):
                        h('target=pc.gid2cell(py.int(py.i))')
                        #if int(h.target.conv.x[int(gidn)]>3): 
                        #    break 
                        for sec in h.target.spk_rx_ls:
                            for seg in sec:
                                #print sec.name(), h(str('z_xtra(')+str(seg.x)+')'), 'rank= ', RANK, '\n'
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
                                #print ' r ', r
                                if(r<15):#I still need to debug the code, to find out why it won't search the finest levels of granularity.
                                    gidcompare='' 

                                    secnames=sec.name()
                                    cellind=int(secnames[secnames.find("Cell[")+5:secnames.find("].")])#This is the index of the post synaptic cell.
                                    #h('py.pol=strcmp(Cell[py.int(py.cellin))]
                                    # post_syn=secnames+' '+'syn_ = new ExpSid('+str(seg.x)+')'
                                    polarity=0
                                    h('py.polarity=Cell[py.int(py.cellind)].num_type==3') 
                                    if(polarity==1):
                                        post_syn=secnames+' '+'syn_ = new GABAa('+str(seg.x)+')'
                                    else:  
                                        post_syn=secnames+' '+'syn_ = new AMPA('+str(seg.x)+')'
                                    h('gidn=0')
                                    h.gidn=int(data[0][1][3])
                        
                                    h(post_syn)
                                    h('print syn_')

                                    h('syn_.cid=py.int(py.i)') #This makes the gid a property of the synapse. 
                                    h('Cell[py.cellind].ampalist.append(syn_)')
                                
                        #cellind=(re.search(r'\[(.*)\]', secnames).group(1))
                                    #print cellind, 'cellind'
                #place the post synapse.
                                    h('synlist.append(syn_)')
                                    h('gidn=0')
                                    h.gidn=int(data[0][1][3])
                                    h('Cell[py.cellind].gvpre.append(gidn)')
                                    h('Cell[py.cellind].div.append(gidn)')
                                    
                        #print True(gidcompare==int(data[0][1][3])), ' is this the same cell?'
                                    h('nc = pc.gid_connect(gidn, syn_)')
                        #connect the pre synapse to the post synapse.
                        #h('print coords.x[3], "this is the global identifier"')
                        #h('Cell[py.int(tail_target)].gvpre.append(gidn)')
                        #h('Cell[py.int(tail_target)].gvpost.append(py.int(py.gidcompare))')
                               
                                    h('print nc," ", nc.srcgid()," ", gidn')
                                    print 'the gids of the connected cells are: ', i, ' ', gidn, '\n'
                                    #quit()
                                    h('nclist.append(nc)')

                                    h('nc.delay = 1+((200+py.r)/delay)')
                                    #h('nc.weight = r')
                                    h('nc.weight = py.r/iw')
#h('hoc_stdout()')

COMM.Barrier()
h('objref tgt')
h('srcid =0')
#h('system("rm net_schematic.txt")')
#h('hoc_stdout("net_schematic.txt")')
COMM.Barrier() 
if RANK==0: print 'source\ttarget\tsynapse\n'
#conn_m=[]
for i in xrange(0,SIZE):
   if i==RANK:
      precell=[]
      postcell=[]
      nclist=0
      h('py.nclist=py.int(nclist.count-1)')
      for j in xrange(0,int(nclist)):
         
         h('srcid=nclist.o(py.j).srcgid()')
         h('tgt = nclist.o(py.j).syn')
         srcid=0
         tgt=0
         
         h('print srcid, " ", tgt.cid')
         h('py.srcid=srcid')
         h('py.tgt=tgt.cid')

         precell.append(int(srcid))
         postcell.append(int(tgt))
      sconm=zip(precell,postcell)
      if np.sum(np.array(sconm))==0: 
        print 'no connections, quiting!'
        exit()
        quit()
      
      #dim1=np.shape(sconm)[0]
      #dim2=np.shape(sconm)[1]
      #print dim1, dim2, 'dimensions'
      #Cl = np.zeros([np.shape(sconm)[0], np.shape(sconm)[1]], dtype='i')

      sconm=COMM.gather(sconm,root=0)

      #if RANK==0:
      #C = np.zeros([dim1, dim2], dtype='i')
      #rowtype = MPI.INT.Create_contiguous(dim2)
      #rowtype.Commit()
  
      #COMM.Gatherv(sendbuf=[sconm, MPI.INT], recvbuf=[C, (sconm, None), rowtype], root=0)

COMM.Barrier()

  


sconm2=[]
if RANK==0:
  
    for k in xrange(0,SIZE):
        sconm2=sconm2+sconm[k]#Concatonate the lists so that they are all the size.

        print sconm2
        conm= [[0 for x in xrange(0,int(NCELL))] for x in xrange(0,int(NCELL))]
        print len(sconm2[:][:])

        for j in xrange(0,len(sconm2[:][:])):
           print sconm2[j][0],sconm2[j][1], sconm2[j] 
           conm[sconm2[j][0]][sconm2[j][1]]=conm[sconm2[j][0]][sconm2[j][1]]+1   
    print conm
    fin =0
    #fin = fin + 1
    fig = plt.figure()
    fig.clf()


    im = plt.imshow(conm, interpolation='nearest')

    plt.autoscale(True)
    plt.colorbar(im)
    plt.xlabel('columns = targets')
    plt.ylabel('rows = sources')
    plt.title('Parallel Conection Matrix')
    plt.grid(True)
        
    sfin = 'Parallel_connection_matrix.png'
    fig.savefig(sfin)
    print sfin

    #h('system("eog Parallel_connection_matrix.png")')
else:
  conm=None
conm = COMM.bcast(conm, root=0)
print RANK, np.shape(conm)


#It may not be necessary to remove this mechanism if pointers are probably set, as they are in rig7.hoc in the serial version. I remember this being a problem.


if RANK==0: 
  def get_in(ma):
    old_row = 0
    row_index = 0
    old_j = 0
    for j in xrange(0, int(ma.shape[0])):
      #print ma[j,:]

      if sum(ma[j, :]) > old_row:  # old :
        old_row = sum(ma[j, :])
        row_index = j
            #print row_index, 'j= ', j, old_row, old_j
        old_j = j
    return row_index
  def get_out(ma):
    old_column = 0
    column_index = 0
    old_i = 0
    for i in xrange(0, int(ma.shape[1])):
      if sum(ma[:, i]) > old_column:  # old :
        old_column = sum(ma[:, i])
        column_index = i
            #print column_index, 'i= ', i, old_column, old_i
        old_i = i
            #print ma
    return column_index

  outdegree=0
  indegree=0
  outdegreeg=0
  indegreeg=0
  def obtain_values(q,r):
    global outdegree, indegree, indegreeg, outdegreeg
    ma=q
    mg=r
    #structural=s
    outdegree = get_out(ma)
    indegree = get_in(ma)
    outdegreeg = get_out(mg)
    indegreeg = get_in(mg)

    #print outdegreeg
    #print indegreeg

    #print outdegree
    #print indegree
    return

  obtain_values(np.array(conm),np.array(conm))
  degrees=[]
  degrees.append(outdegreeg)
  degrees.append(indegreeg)  
  degrees.append(outdegree)  
  degrees.append(indegree)
else:
  degrees=None
degrees = COMM.bcast(degrees, root=0)
print degrees, RANK
h('xopen("instrumentation.hoc")')
#This procedure will be crippled by the addition of NetStim netconns.
#Such that tracenet has to be executed before the netstims.


h('xopen("rigp.hoc")')#Check for strength of external stimulation.
h('forall{ for(x,0){ uninsert xtra }}') #Parallel Simulation causes seg fault wi

COMM.Barrier() 
if RANK==0: print '\n cell numbers on the various hosts \n'
for i in xrange(0,SIZE): 
    if (i==RANK):
        print 'host', RANK, '\n'
        h('print "cells.count ", cells.count')
COMM.Barrier() 
h('spikerecordg()')#has to be called before the run only, otherwise it will delet
COMM.Barrier()

prun(tstop)

COMM.Barrier()
h('spikeout()')
COMM.Barrier()

timec=[]
ncell=0.
h('py.ncell=cells.count-1')
for i in xrange(0,int(ncell)):
   timec.append(h.recvectors[int(i)].to_python())
   #print h.recvectors[int(i)].to_python()

timec=COMM.gather(timec,root=0)
COMM.Barrier()
if RANK==0:
   timec2=[]
   for k in xrange(0,SIZE):
        timec2=timec2+timec[k]#Concatonate the lists so that they are all the size.

        print timec2
   print timec2
   tc = np.array(timec2[int(0)])
   N = len(tc)
   t = np.linspace(0.0, 0.025 * N, N)

   fin += 1
   fig = plt.figure()
   fig.clf()
   plt.hold(True)
   for i in xrange(0, int(len(timec2))):
     tc = np.array(timec2[int(i)])
     plt.plot(t, tc, linewidth=3)

     plt.title('cells index membrane potential')
     plt.xlabel('ms')
     plt.ylabel('mV')
    
   # out1=downsample(time_courses[int(i)],oldrate=40000,newrate=200)
   # downsampled.append(out1)
   plt.hold(False)
   sfin = 'traces_all'+str(RANK)+'.png'
   fig.savefig(sfin)
   #h('system("eog traces_all0.png")')


COMM.Barrier()
h('spikeout()')
COMM.Barrier()
h('xopen("post_analysisp.hoc")')
h('objref vb[2]')
for s in xrange(0,SIZE):
    for j in xrange(0,NCELL):
        #for i, cells.count-1{
        #h.cells.o(j)
        if (RANK==s):
            #coordlist=[[0 for x in xrange(0,2)] for x in xrange(0,2)]
          histogram=[]  
          test=0
            #print j, ' j ', s, ' s '
          h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
          if (test!=0):
              h('vb[0] = new Vector() //vb stands for vector binned.')
              h('vb[0].hist(times[i],0,(maxt+binsz-1)/binsz,binsz)')
              h('vb[0].x[0]=j')
              h('py.histogram=vb[0].hist(times[i],0,(maxt+binsz-1)/binsz,binsz).to_python()')
              print histogram, 'not indexed'
              #print histogram[0], 'index1'
          else :
             histogram=None
        else: 
             histogram=None
        #Only gid 6 is only every being broadcast.
        data = COMM.bcast(histogram, root=s)#ie root = rank

        if (data!=None):#and(RANK!=s):
            for i in xrange(0,NCELL):
                gidn=int(data[0])
                #if(gidn==6):
                #    print coordlist
                #    print 'yes, not equal to gid 6'
                #    break
                if(i!=int(gidn)):#if the gids are not the same
                    test2=0
                    #print j, ' j ', i, ' i ', s, ' s '
                    h('py.test2=py.int(pc.gid_exists(py.int(py.i)))')
                    if (test2!=0):

                      h('vb[1] = new Vector() //vb stands for vector binned.')
                      #h('vb[0].hist(times[i],0,(maxt+binsz-1)/binsz,binsz)               ')
                      h('vb[1].hist(times[j],0,(maxt+binsz-1)/binsz,binsz)')
                      
                      h('vb[1].x[0]=0')    
                      h('vo1=normte(py.data,vb[1],20)')
                      h('d=vo1.x(2)')
                      h('py.indexi=py.int(i)')
                      h('py.indexj=py.int(j)')
                      h('py.cellgidi=cells.o(i).gid1')
                      h('py.cellgidj=cells.o(j).gid1')
                      h(' storval=GetTENQ(vb[0],vb[1],20,int(py.i),int(py.data[0]))')
                      h(' nqte3.append(storval)')
                      h('//storval.v[5]')
                      h('// from(0) to(1) TE(2) NTE(3) HX2|X2P(4) prefdir(5) TEshufavg(6) TEshufstd(7) sig(8)')	  
                      h('py.dirin=py.float(storval.v[5].x[0])')
                      h('py.ent=py.float(storval.v[3].x[0])')
                      h('target=pc.gid2cell(py.int(py.i))')
               
