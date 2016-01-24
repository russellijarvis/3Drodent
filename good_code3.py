#!/usr/bin/env python
'''
LFPs from a population of cells relying on MPI
'''
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
#pc.done()

h.dt = 0.025
tstop = 1025
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
NCELL=16

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
h.xopen("/home/zaza3/Downloads/trunk/examples/expericomp20140421/morph4.hoc")
h.xopen("/home/zaza3/Downloads/trunk/examples/expericomp20140421/nqs.hoc")

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


h('objref py')
h('py = new PythonObject()')
cnt=0 #This value cnt, will be different on different hosts and thats what is required.

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
for i in range(RANK, NCELL, SIZE): #20 was int(len(allrows))
    s=allrows2[i]    
    storename= str(s[3]) #//simply being inside a loop, may be the main pro
    if(re.search('.swc',storename)):                                    
        h.cell=h.mkcell(storename)
        h('cell.geom_nseg()')
        h('cell.pyr()')
        h('strdef cellposition')
        h('sprint(cellposition,"%s%d%s%d%s%d%s","cell.position(",py.int(py.s[0]),",",py.int(py.s[1]),",",py.int(py.s[2]),")")')
        h('execute(cellposition)')
	
        h('pc.set_gid2node(py.i, pc.id)')# // associate gid i with this host
        h('objref nc')
        h('nc = cell.connect2target(nil)')# // attach spike detector to cell
        h('pc.cell(py.i, nc)')#//                    
        h('cells.append(cell)')
        h('gidvec.append(py.i)')
        gidvec.append(i)
        h('print py.i')
        cnt+=1
                
             




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

i=0
cnt=0
cnt_sec=0
h('objref target')
for s in xrange(0,SIZE):
    for j in xrange(0,NCELL-1):
        #for i, cells.count-1{
        #h.cells.o(j)
        if (RANK==s):
            coordlist=[[0 for x in xrange(0,2)] for x in xrange(0,2)]
     
            test=0
            print j, ' j ', s, ' s '
            h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
            if (test!=0):
                h('target=pc.gid2cell(py.int(py.j))')
                for sec in h.target.spk_trig_ls:
                    for seg in sec:
                        #print seg.x
                        #for(x,0){
                        #print sec.name(), 'got here', 'test ', test
                        get_cox=str('x_xtra('+str(seg.x)+')')
                        #get_cox=str('coords.x[0]='+sec.name()+'.x_xtra()')
                        h(get_cox)
                        get_coy=str('y_xtra('+str(seg.x)+')')
                        #get_coy=str('coords.x[1]='+sec.name()+'.y_xtra()')
                        h(get_coy)
                        get_coz=str('z_xtra('+str(seg.x)+')')
                        #get_coz=str('coords.x[2]='+sec.name()+'.z_xtra()')
                        h(get_coz)
                        #h('print x_xtra('+str(seg.x)+')')
                        h('coords.x[3]=py.int(j)')#ie the gidvec.  
                        h('coords.x[4]=py.seg.x')#ie the gidvec.  
 
                        coords=np.array(h.coords.to_python(),dtype=np.float32)
                    #print 's: ', s , 'coords', coords, '\n'
                        coordlist[0][1]=coords
                        h('strdef secnames')
                        h('cells.o(0).soma[0]{ secnames=secname() }')
                        secnames=sec.name() #h.secnames
                        coordlist[1][0]=str(secnames)
            if (test==0):
                coordlist=None
        else: 
            coordlist=None
        
        data = COMM.bcast(coordlist, root=s)#ie root = rank

        if (data!=None)and(RANK!=s):
            for sec in h.spk_rx_ls:#iterate through post synaptic locations.
                #print "bcast finished and data on RANK %d is: "%RANK, data[0][1], data[1][0]#, secnames,' data ', data
                for seg in sec:
                    #print seg.x
                    #h('coords2.x[0]=cells.o(0).soma[0].x_xtra()')
                    #h('coords2.x[1]=cells.o(0).soma[0].y_xtra()')
                    h(str('coords2.x[2]='+sec.name()+'.z_xtra()'))
                    h(str('coords2.x[1]='+sec.name()+'.y_xtra()'))
                    h(str('coords2.x[0]='+sec.name()+'.x_xtra()'))



                    #h('coords2.x[0]=cells.o(0).soma[0].x_xtra()')
                    #h('coords2.x[1]=cells.o(0).soma[0].y_xtra()')
                    #h('coords2.x[2]=cells.o(0).soma[0].z_xtra()')
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
                    if(r<10):
                    #print r, RANK, ' r, RANK '
                    #h('print r, "this is r"')
                #Assign a post synapse.
                #secnames=data[1][0][:]
                        secnames=sec.name()
                        post_syn=secnames+' '+'syn_ = new AMPA('+str(seg.x)+')'
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
                        h('nclist.append(nc)')
                        h('print nclist.count, " size nclist"')
COMM.Barrier()


COMM.Barrier() 
if RANK==0: print '\n gidvecs on the various hosts \n'
for i in xrange(0,SIZE): 
    if (i==RANK):
        print 'host', RANK, gidvec[:], '\n'
        h('gidvec.printf()')
    
COMM.Barrier() 



COMM.Barrier() 
if RANK==0: print '\n cell numbers on the various hosts \n'
for i in xrange(0,SIZE): 
    if (i==RANK):
        print 'host', RANK, '\n'
        h('print "cells.count ", cells.count')
COMM.Barrier() 

h('indegree=4')
h('outdegree=6')
indegree=4
outdegree=6
h('xopen("rigp.hoc")')
h
def runp(tstop):
    lfp=1
    if lfp == 1:
        h.init()
    if lfp != 1:
    	h.finitialize()  # set initial variables.
    while h.t < tstop:  # Perform integration over the cable equation.
        if lfp == 1:
            h.advance()  # this method is designed to work with field and may take longer.  
            pc.psolve(1000)
        if lfp != 1:
            h.fadvance()
            pc.psolve(1000)



runp(1000)
h('pc.psolve(1000)')

                
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
        

        
