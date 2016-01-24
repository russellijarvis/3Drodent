#!/usr/bin/python
# -*- coding: utf-8 -*-

# All original code, that appeals to idioms described by Hines and Carnevale.

# Copyright (C) 2012, 2013 Russell Jarvis
# This file is part of Open Blue Brain.
#
# Open Blue Brain is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Open Blue Brain is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License

import pdb
import glob
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


import socket
from matplotlib.collections import PolyCollection, LineCollection
import os
import urllib2
import zipfile
from mpi4py import MPI

#import random1 as r
#rand=r.random1


from neuron import h

# Dont import LFPy
#Importing LFPy imports neuron. NEURON has to be imported after mpi4py
hostname=socket.gethostname()

# initialize the MPI interface

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

NFILE = 3175#1175 # Maximum on this machine.
#Number of files to consider. This always ends up being a number more than the number of cells used. However NCELL and NFILE are of a similar magnitude. Some cell files have to filtered out, and cannot be used.

h('ew=0.0185')
h('delay=3000')
h('iw=1000 // divisor/ or scaler, so if its bigger weight should be smaller. Was 1000 but this should be stronger. delay equals internal weight.'
  )

'''
h('iwp=1')

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
'''
# h('if(qd==1){prunenet=900}')

h('objref py')
h('py = new PythonObject()')

###

tstop = 1500#550
h('tstop1=1500')

f = open('tstopfile', 'w')
f.write(str(tstop))
f.close()
###

#h('minr=10')  # minumum r large, for some reason which I do not yet
#h('minrs=10')  # minumum r small, understand
#h('delay=3000')

#internal weight. Average scaler of weights internal to the network.
h('iw=1000 // divisor/ or scaler, so if its bigger weight should be smaller. Was 1000 but this should be stronger. delay equals internal weight.'
  )
#external weight. Wheight applied to external inputs.
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
def gather(data):

    # https://github.com/NeuralEnsemble/PyNN/blob/master/examples/distrib_example.py
    # Andrew Davison

    assert isinstance(data, np.ndarray)

    # first we pass the data size

    size = data.size
    sizes = COMM.gather(size, root=0) or []

    # now we pass the data

    displacements = [sum(sizes[:i]) for i in range(len(sizes))]
    print (COMM.rank, 'sizes=', sizes, 'displacements=', displacements)
    gdata = np.empty(sum(sizes))
    COMM.Gatherv([data, size, MPI.DOUBLE], [gdata, (sizes,
                 displacements), MPI.DOUBLE], root=0)
    return gdata
'''
def broadcast(RANK,NCELL,SIZE,h,):

    global s,j,i,test
    s=0
    j=0
    i=0
    gidcompare = ''    

    for s in xrange(0, SIZE):#Was full SIZE, not SIZE-1
        for j in xrange(0,NCELL):
            coordlist=[]
            if RANK == s:
                coordlist = [[0 for x in xrange(0, 2)] for x in xrange(0,
                         2)]
                test = 0
                h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
                if test != 0:
                    h('source=pc.gid2cell(py.int(py.j))')
                   
                    #Here execute a Lambda function. A function that I don't know what it is yet but it
                    #Defines what coorlist is
                    #coordlist = 
                if test == 0:
                    coordlist = None
            else:
                coordlist = None
            data = COMM.bcast(coordlist, root=s)  # ie root = rank

            if data != None:  
                for i in xrange(0, NCELL-1):
                    gidn = int(data[0][1][3])
                    test2 = 0
                    h('py.test2=py.int(pc.gid_exists(py.int(py.i)))')
                    if test2 != 0:
                        h('target=pc.gid2cell(py.int(py.i))')
                        #Here execute a Lambda function. A function that I don't know what it is yet but it
                        # But it connects up neurons, or calculates the transfer entropy.

    COMM.Barrier()

    return (h.nclist, ecm, icm)
'''

def pprint(str='', end='\n', comm=MPI.COMM_WORLD):
    """Print for MPI parallel programs: Only rank 0 prints *str*."""

    if comm.rank == 0:
        print str + end,


h('proc pprint(){ if (pc.id()==0){ print $s1 }}')


def prun(tstop):
    #This code is from:
    #http://senselab.med.yale.edu/ModelDB/ShowModel.asp?model=151681
    cvode = h.CVode()
    cvode.cache_efficient(1)

  # pc.spike_compress(0,0,1)

    pc.setup_transfer()
    mindelay = pc.set_maxstep(10)
    if RANK == 0:
        print 'mindelay = %g' % mindelay
    runtime = h.startsw()
    exchtime = pc.wait_time()

    inittime = h.startsw()
    h.stdinit()
    inittime = h.startsw() - inittime
    if RANK == 0:
        print 'init time = %g' % inittime

    while h.t < tstop:
        told = h.t
        tnext = h.t + checkpoint_interval
        if tnext > tstop:
            tnext = tstop
        pc.psolve(tnext)
        if h.t == told:
            if RANK == 0:
                print 'psolve did not advance time from t=%.20g to tnext=%.20g\n' \
                    % (h.t, tnext)
            break

    # if h.t%2==0: The problem is h.t is float multiple not integer multiple

        print 'working', h.t
    runtime = h.startsw() - runtime
    comptime = pc.step_time()
    splittime = pc.vtransfer_time(1)
    gaptime = pc.vtransfer_time()
    exchtime = pc.wait_time() - exchtime
    if RANK == 0:
        print 'runtime = %g' % runtime
    print comptime, exchtime, splittime, gaptime


  # printperf([comptime, exchtime, splittime, gaptime])/

def stationary_poisson(
    #This code is from LFPy.
    nsyn,
    lambd,
    tstart,
    tstop,
    ):
    ''' Generates nsyn stationary possion processes with rate lambda between tstart and tstop'''

    interval_s = (tstop - tstart) * .001
    spiketimes = []
    for i in xrange(nsyn):
        spikecount = np.random.poisson(interval_s * lambd)
        spikevec = np.empty(spikecount)
        if spikecount == 0:
            spiketimes.append(spikevec)
        else:
            spikevec = tstart + (tstop - tstart) \
                * np.random.random(spikecount)
            spiketimes.append(np.sort(spikevec))  # sort them too!

    return spiketimes


# Test of analysis functions.

allrows = pickle.load(open('allrows.p', 'rb'))
h.xopen('morph4.hoc')
h.xopen('nqs.hoc')

os.system('pwd')
import re

cnt1 = pc.id

# I do not know how to refer to relative paths in Python,
# the below emulates a call to a relative path.

os.chdir(os.getcwd() + '/main')


def prep_list(NFILE):
    cnt = 0
    allrows2 = []
    i = 0
    for i in xrange(0, NFILE - 1):
        s = allrows[i]
        if cnt > 1:
            if int(len(s)) > 9:  # //This condition is counter intuitive many cells
                storename = str(s[3])  # //simply being inside a loop, may be the main pro

                allrows2.append(allrows[i])
        cnt += 1
    return allrows2
    print np.shape(allrows2), np.shape(allrows)
allrows2 = prep_list(NFILE)


###
#
NCELL=int(np.shape(allrows2)[0])
# Should be max_NCELL=int(np.shape(allrows2)[0])
#
###

   

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




#GABA neuro transmitters are made by this code.
# start, end, step size.
#If size =1, then step size is 0.
#If size =2, then step size is 1.
def mb(RANK, NCELL, SIZE, allrows2, gidvec, h,soff,boff):
#make soff and boff optional parameters in python function.

    #h('objref py')
    #h('py = new PythonObject()')
    global i, cnt, cnti,s
    i=0
    cnti = 0
    cnt=0
    s=''
    storename=''
    for i in range(RANK+soff, boff,SIZE):#NCELL-1, SIZE):  # 20 was int(len
    
        s = allrows2[i]
        storename = str(s[3])  # //simply being inside a loop, may be the main problem
        if re.search('.swc', storename):
            h.cell = h.mkcell(storename)


            h('cell.geom_nseg()')
            h('cell.gid1=py.i')
            h('cell.gvpre.printf')
            h('cell.gvpost.printf')




            h.cell.nametype=str(s[5])
            #h('cell.nametype=py.str(py.s[5])')
            h.cell.num_type=int(s[6])
            #h('cell.num_type=py.int(py.s[6])')
            h('cell.population=py.str(py.s[7])')
            #h('cell.reponame=py.str(py.storename)')
            h.cell.reponame=str(storename)
            h('cell.div.resize(py.int(py.NCELL))')
            h('cell.conv.resize(py.int(py.NCELL))')
            #h('cell.nametype=py.str(py.s[5])')

            h('if(strcmp("pyramid",py.s[5])==0){pyr_list.append(cell)}')
            h('if(strcmp(cell.population,"neocortex")==0){ if(strcmp("pyramid",py.s[5])==0){ cell.pyr() }}')


            h('if(strcmp(cell.population,"hippocampus")==0){ if(strcmp("pyramid",py.s[5])==0){ cell.pyr2() }}')

            if 'hippocampus' in s:
                h('cell.pyr2()')


            if 'interneuron' in s:
                ie0[cnti] = i
                ie1[cnti] = 1
                cnti += 1
                h.cell.polarity=int(1)
                h('cell.polarity=1')

            if 'pyramid' in s:
                ie0[cnti] = i
                ie1[cnti] = 0
                cnti += 1
                h.cell.polarity=int(0)
                h('cell.polarity=0')
            print s[5]
            h('print py.s[5]')
            h('if(strcmp("interneuron",py.s[5])==0){ cell.polarity=1 }')
            h('if(strcmp("pyramid",py.s[5])==0){ cell.polarity=0 }')
            
                

            #h('cell.polarity=1')

            h('if(strcmp("interneuron",py.s[5])==0){ inter_list.append(cell) }')
            h('if(strcmp("interneuron",py.s[5])==0){ cell.basket()}')
            h('if(strcmp("interneuron",py.s[5])==0){ print "inter neuron"}')

            h('if(strcmp("aspiny",py.s[5])==0){ aspiny_list.append(cell) }')

            h('if(strcmp(cell.population,"hippocampus")==0){ hipp.append(cell)  }')

            h('if(strcmp(cell.population,"neocortex")==0){ neoc.append(cell) }')

            #h('strdef cellposition')
            #h('sprint(cellposition,"%s%d%s%d%s%d%s","cell.position(",py.float(py.s[0]),",",py.float(py.s[1]),",",py.float(py.s[2]),")")')
            #h('print cellposition')
            #h('execute(cellposition)')

            h('pc.set_gid2node(py.i, pc.id)')  # // associate gid i with this host

            h('cell.soma[0] nc =  new NetCon(&v(0.5), nil)')
            h('pc.cell(py.i, nc)')  # //


            h('cells.append(cell)')
            h('gidvec.append(py.i)')
            gidvec.append(i)
            h('print py.i')
            cnt += 1
    #destroy the allrows list to free up memory, then return the empty list.
    #I can't seem to destroy allrows here without wrecking something in rigp later.
    #allrows=[]
    return (h.cells, allrows2)
#import brain_functions as bf

#h.cells, allrows2, ie0, ie1=bf.mb(RANK, NCELL, SIZE, allrows2, gidvec,h,1150,1250)
h.cells, allrows=mb(RANK, NCELL, SIZE, allrows2, gidvec,h,1250,1350)
h('for i=0,cells.count-1{ print cells.o(i).polarity }')
#raise SystemExit(0)
allrows2=[]

#for i in range(RANK+1150, 1350,SIZE):#NCELL-1, SIZE):  # 20 was int(len(allrows))


'''
for i in range(RANK, NCELL-1, SIZE):  # 20 was int(len(allrows))
    s = allrows2[i]
    storename = str(s[3])  # //simply being inside a loop, may be the main problem
    if re.search('.swc', storename):
        h.cell = h.mkcell(storename)


        h('cell.geom_nseg()')
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

        h('if(strcmp("pyramid",py.s[5])==0){pyr_list.append(cell)}')
        h('if(strcmp(cell.population,"neocortex")==0){ if(strcmp("pyramid",py.s[5])==0){ cell.pyr() }}')


        h('if(strcmp(cell.population,"hippocampus")==0){ if(strcmp("pyramid",py.s[5])==0){ cell.pyr2() }}')

        if 'hippocampus' in s:
            h('cell.pyr2()')


        if 'interneuron' in s:
            ie0[cnti] = i
            ie1[cnti] = 1
            cnti += 1
            h('cell.polarity=1')
        if 'pyramid' in s:
            ie0[cnti] = i
            ie1[cnti] = 0
            cnti += 1
            h('cell.polarity=0')

        # h('if(strcmp("interneuron",py.s[5])==0){py.ie[py.i][0]=py.i}')
        # h('if(strcmp("interneuron",py.s[5])==0){py.ie[py.i][1]=0}')

        h('if(strcmp("interneuron",py.s[5])==0){ inter_list.append(cell) }')
        h('if(strcmp("interneuron",py.s[5])==0){ cell.basket()}')
        h('if(strcmp("interneuron",py.s[5])==0){ print "inter neuron"}')

        h('if(strcmp("aspiny",py.s[5])==0){ aspiny_list.append(cell) }')

        h('if(strcmp(cell.population,"hippocampus")==0){ hipp.append(cell)  }')

        h('if(strcmp(cell.population,"neocortex")==0){ neoc.append(cell) }')

        h('strdef cellposition')
        h('sprint(cellposition,"%s%d%s%d%s%d%s","cell.position(",py.int(py.s[0]),",",py.int(py.s[1]),",",py.int(py.s[2]),")")')

        # h('execute(cellposition)')

        h('pc.set_gid2node(py.i, pc.id)')  # // associate gid i with this host

        # h('nc = cell.connect2target(nil)')# // attach spike detector to cell

        h('cell.soma[0] nc =  new NetCon(&v(0.5), nil)')
        h('pc.cell(py.i, nc)')  # //

        # cell=h.cell
        # pc.cell(i, h.NetCon(cell.soma(.5)._ref_v, None, sec=cell.soma))

        h('cells.append(cell)')
        h('gidvec.append(py.i)')
        gidvec.append(i)
        h('print py.i')
        cnt += 1
''' 

print np.shape(allrows2)[0], NCELL, h('cells.count')

# mbs(allrows2)
# h('pc.psolve(1000)')

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

#NCELL=boff-soff
#To make this code work.
icm = np.zeros((NCELL, NCELL))
ecm = np.zeros((NCELL, NCELL))

i = 0
j = 0
cnt = 0
cnt_sec = 0

#for cell in h.cells:
#    h('print cell.gid1')



def wirecells(RANK,NCELL,SIZE,h,icm,ecm):

    global s,j,i,test,test2,r
    s=0
    j=0
    i=0
    gidcompare = ''

    secnames = ''# sec.name()
    cellind =0 #int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.

    polarity = 0
    

    for s in xrange(0, SIZE):#Was full SIZE, not SIZE-1
        for j in xrange(0,NCELL):
            print s, j
            if RANK == s:
                coordlist = [[0 for x in xrange(0, 2)] for x in xrange(0,
                         2)]

                test = 0

            
                h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
                if test != 0:
                    h('source=pc.gid2cell(py.int(py.j))')
                    for sec in h.source.spk_trig_ls:
                        for seg in sec:
                            get_cox = str('coords.x[0]=x_xtra('
                                          + str(seg.x) + ')')
                            h(get_cox)


                            get_coy = str('coords.x[1]=y_xtra('
                                          + str(seg.x) + ')')
                            h(get_coy)
                            get_coz = str('coords.x[2]=z_xtra('
                                          + str(seg.x) + ')')
                            h(get_coz)
                            h.coords.x[3] = int(j)
                            h.coords.x[4] = seg.x # ie the gidvec
                            #h('coords.x[4]=py.seg.x')  # ie the gidvec.

                            coords = np.array(h.coords.to_python(),
                                              dtype=np.float64)

                            coordlist[0][1] = coords


                            secnames = sec.name()  # h.secnames
                            coordlist[1][0] = str(secnames)


                if test == 0:
                    coordlist = None
            else:
                coordlist = None
            data = COMM.bcast(coordlist, root=s)  # ie root = rank

            if data != None:  
                for i in xrange(0, NCELL-1):
                    gidn = int(data[0][1][3])
                    
                    if i != int(gidn):  # if the gids are not the same
                        test2 = 0
                        h('py.test2=py.int(pc.gid_exists(py.int(py.i)))')
                        if test2 != 0:
                            h('target=pc.gid2cell(py.int(py.i))')
                            for sec in h.target.spk_rx_ls:
                                for seg in sec:
                                    h(str('coords2.x[2]=') + str('z_xtra(')
                                      + str(seg.x) + ')')
                                    h(str('coords2.x[1]=') + str('y_xtra(')
                                      + str(seg.x) + ')')
                                    h(str('coords2.x[0]=') + str('x_xtra(')
                                    + str(seg.x) + ')')

                                    h('coordsx=0.0')
                                    h.coordsx = data[0][1][0]
                                    h('coordsy=0.0')
                                    h.coordsy = data[0][1][1]
                                    h('coordsz=0.0')
                                    h.coordsz = data[0][1][2]
                                    r = 0.
                                    h('py.r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)')
                                    r = float(r)
                                    if r < 10:  
                                        gidcompare = ''

                                        secnames = sec.name()
                                        cellind = int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.

                                        polarity = 0
                                   

                                        #cellind is a cell index, that is relative to the host. So the identifier repeats on different hosts.
                                        #gidn is a global identifier. These numbers are not repeated on different hosts.
                                        polarity=int(h.Cell[int(cellind)].polarity)
                                        
                                        print polarity

                                        if int(polarity) == int(1):
                                            post_syn = secnames + ' ' + 'syn_ = new GABAa(' + str(seg.x) + ')'
                                            icm[i][gidn] = icm[i][gidn] + 1
                                        else:

                                            post_syn = secnames + ' ' + 'syn_ = new AMPA(' + str(seg.x) + ')'
                                            ecm[i][gidn] = ecm[i][gidn] + 1

                                        h('gidn=0')
                                        h.gidn = int(data[0][1][3])

                                        h(post_syn)
                                        h('print syn_')

                                        h('syn_.cid=py.int(py.i)')  # This makes the gid a property of the synapse.
                                        h.Cell[cellind].ampalist.append(h.syn_)



                                        h('synlist.append(syn_)')
                                        h('gidn=0')
                                        h.gidn = int(data[0][1][3])
                                        h.Cell[cellind].div.append(h.gidn)
                                        h.Cell[cellind].gvpre.append(h.gidn)


                                        h('nc = pc.gid_connect(gidn, syn_)')

             

                                        h('print nc," ", nc.srcgid()," ", gidn')
                                        print 'the gids of the connected cells are: ', i, ' ', gidn, '\n'


                                        h('nclist.append(nc)')

                                        h('nc.delay = 1+((py.r)/delay)')
                                        h('nc.weight = py.r/iw')

    COMM.Barrier()

    return (h.nclist, ecm, icm)





def wirecells2(RANK,NCELL,SIZE,h,icm,ecm):

    global s,j,i,test,test2,r
    s=0
    j=0
    i=0
    gidcompare = ''

    secnames = ''# sec.name()
    cellind =0 #int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.

    polarity = 0
    

    for s in xrange(0, SIZE):#Was full SIZE, not SIZE-1
        for j in xrange(0,NCELL):
            print s, j
            if RANK == s:
                coordlist = [[0 for x in xrange(0, 2)] for x in xrange(0,
                         2)]

                test = 0

            
                h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
                if test != 0:
                    h('source=pc.gid2cell(py.int(py.j))')
                    for sec in h.source.spk_trig_ls:
                        for seg in sec:
                            get_cox = str('coords.x[0]=x_xtra('
                                          + str(seg.x) + ')')
                            h(get_cox)


                            get_coy = str('coords.x[1]=y_xtra('
                                          + str(seg.x) + ')')
                            h(get_coy)
                            get_coz = str('coords.x[2]=z_xtra('
                                          + str(seg.x) + ')')
                            h(get_coz)
                            h.coords.x[3] = int(j)
                            h.coords.x[4] = seg.x # ie the gidvec
                            #h('coords.x[4]=py.seg.x')  # ie the gidvec.

                            coords = np.array(h.coords.to_python(),
                                              dtype=np.float64)

                            coordlist[0][1] = coords


                            secnames = sec.name()  # h.secnames
                            coordlist[1][0] = str(secnames)


                if test == 0:
                    coordlist = None
            else:
                coordlist = None
            data = COMM.bcast(coordlist, root=s)  # ie root = rank

            if data != None:  
                for i in xrange(0, NCELL-1):
                    gidn = int(data[0][1][3])
                    
                    #if i != int(gidn):  # if the gids are not the same
                    test2 = 0
                    h('py.test2=py.int(pc.gid_exists(py.int(py.i)))')
                    if test2 != 0:
                        h('target=pc.gid2cell(py.int(py.i))')

                        for sec in h.target.spk_rx_ls:
                            for seg in sec:
                                #Get the coordinates of locations on this machine.
                                h(str('coords2.x[2]=') + str('z_xtra(')
                                  + str(seg.x) + ')')
                                h(str('coords2.x[1]=') + str('y_xtra(')
                                  + str(seg.x) + ')')
                                h(str('coords2.x[0]=') + str('x_xtra(')
                                  + str(seg.x) + ')')
                                #Unpack the broadcast coordinates.
                                #Coordinates of locations possibly on a different machine.
                                h('coordsx=0.0')
                                h.coordsx = data[0][1][0]
                                h('coordsy=0.0')
                                h.coordsy = data[0][1][1]
                                h('coordsz=0.0')
                                h.coordsz = data[0][1][2]
                                secnamesg=data[1][0]
                                #The idea is that the string secnamesg cannot be the same string as sec.name() however this needs to be done with a string comparison.
                                #Allow self wiring of neurons, but not on the same segment/section. Because that is a short circuit.
                                #if secnamesg!=sec.name():

                                #if re.search('.swc', storename):

                                if re.search(str(secnamesg), str(sec.name())):
                                    print secnamesg, sec.name()

                                    r = 0.
                                    h('py.r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)')
                                    r = float(r)
                                    if r < 10:  
                                        gidcompare = ''

                                        secnames = sec.name()
                                        cellind = int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.

                                        polarity = 0
                                   

                                        #cellind is a cell index, that is relative to the host. So the identifier repeats on different hosts.
                                        #gidn is a global identifier. These numbers are not repeated on different hosts.
                                        polarity=int(h.Cell[int(cellind)].polarity)
                                        
                                        print polarity

                                        if int(polarity) == int(1):
                                            post_syn = secnames + ' ' + 'syn_ = new GABAa(' + str(seg.x) + ')'
                                            icm[i][gidn] = icm[i][gidn] + 1
                                        else:

                                            post_syn = secnames + ' ' + 'syn_ = new AMPA(' + str(seg.x) + ')'
                                            ecm[i][gidn] = ecm[i][gidn] + 1

                                        h('gidn=0')
                                        h.gidn = int(data[0][1][3])

                                        h(post_syn)
                                        h('print syn_')

                                        h('syn_.cid=py.int(py.i)')  # This makes the gid a property of the synapse.
                                        h.Cell[cellind].ampalist.append(h.syn_)



                                        h('synlist.append(syn_)')
                                        h('gidn=0')
                                        h.gidn = int(data[0][1][3])
                                        h.Cell[cellind].div.append(h.gidn)
                                        h.Cell[cellind].gvpre.append(h.gidn)


                                        h('nc = pc.gid_connect(gidn, syn_)')

             

                                        h('print nc," ", nc.srcgid()," ", gidn')
                                        print 'the gids of the connected cells are: ', i, ' ', gidn, '\n'


                                        h('nclist.append(nc)')

                                        h('nc.delay = 1+((py.r)/delay)')
                                        h('nc.weight = py.r/iw')

    COMM.Barrier()

    return (h.nclist, ecm, icm)
h.nclist=wirecells2(RANK,NCELL,SIZE,h,icm,ecm)
#h.nclist=wirecells2(RANK,NCELL,SIZE,h,icm,ecm)


#wirecells(RANK,NCELL,SIZE,h)
'''
#if SIZE=nhosts=4, to iterate all over them, must stop at 3.
#This would mean ranks, 0,1,2,3 were visited, 4 ranks.
for s in xrange(0, SIZE-1):#Was full SIZE, not SIZE-1
    for j in xrange(0, NCELL-1):
        # for i, cells.count-1{
        # h.cells.o(j)

        if RANK == s:
            coordlist = [[0 for x in xrange(0, 2)] for x in xrange(0,
                         2)]

            test = 0

            # print j, ' j ', s, ' s '

            h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
            if test != 0:
                h('source=pc.gid2cell(py.int(py.j))')
                for sec in h.source.spk_trig_ls:
                    for seg in sec:
                        get_cox = str('coords.x[0]=x_xtra('
                                + str(seg.x) + ')')
                        h(get_cox)

                        # h('print py.get_cox')

                        get_coy = str('coords.x[1]=y_xtra('
                                + str(seg.x) + ')')
                        h(get_coy)
                        get_coz = str('coords.x[2]=z_xtra('
                                + str(seg.x) + ')')
                        h(get_coz)
                        h.coords.x[3] = int(j)
                        h('coords.x[4]=py.seg.x')  # ie the gidvec.

                        coords = np.array(h.coords.to_python(),
                                dtype=np.float32)

                    # print 's: ', s , 'coords', coords, '\n'

                        coordlist[0][1] = coords

                        # h('strdef secnames')
                        # h('cells.o(0).soma[0]{ secnames=secname() }')

                        secnames = sec.name()  # h.secnames
                        coordlist[1][0] = str(secnames)

                        # print coordlist

            if test == 0:
                coordlist = None
        else:
            coordlist = None

        # Only gid 6 is only every being broadcast.

        data = COMM.bcast(coordlist, root=s)  # ie root = rank

        if data != None:  # and(RANK!=s):
            for i in xrange(0, NCELL-1):
                gidn = int(data[0][1][3])
                if i != int(gidn):  # if the gids are not the same
                    test2 = 0
                    h('py.test2=py.int(pc.gid_exists(py.int(py.i)))')
                    if test2 != 0:
                        h('target=pc.gid2cell(py.int(py.i))')
                        for sec in h.target.spk_rx_ls:
                            for seg in sec:
                                h(str('coords2.x[2]=') + str('z_xtra(')
                                  + str(seg.x) + ')')
                                h(str('coords2.x[1]=') + str('y_xtra(')
                                  + str(seg.x) + ')')
                                h(str('coords2.x[0]=') + str('x_xtra(')
                                  + str(seg.x) + ')')

                                h('coordsx=0.0')
                                h.coordsx = data[0][1][0]
                                h('coordsy=0.0')
                                h.coordsy = data[0][1][1]
                                h('coordsz=0.0')
                                h.coordsz = data[0][1][2]
                                r = 0.
                                h('py.r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)'
                                  )
                                r = float(r)
                                if r < 15:  # I still need to debug the code, to find out why it won't search the finest levels of granularity.
                                    gidcompare = ''

                                    secnames = sec.name()
                                    cellind = \
    int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.

                                    polarity = 0
                                    #h('py.polarity=py.int(Cell[py.int(py.cellind)].num_type==3)'
        #)

                                    h('py.polarity=py.int(Cell[py.int(py.cellind)].polarity)')

                                    if int(polarity) == int(1):
                                        post_syn = secnames + ' ' \
    + 'syn_ = new GABAa(' + str(seg.x) + ')'
                                        icm[i][gidn] = icm[i][gidn] + 1
                                    else:

                                        # break;
                                        # if num_type, is 4, or 5 then the cell type is excitatory

                                        post_syn = secnames + ' ' \
    + 'syn_ = new AMPA(' + str(seg.x) + ')'
                                        ecm[i][gidn] = ecm[i][gidn] + 1

                                    h('gidn=0')
                                    h.gidn = int(data[0][1][3])

                                    h(post_syn)
                                    h('print syn_')

                                    h('syn_.cid=py.int(py.i)')  # This makes the gid a property of the synapse.
                                    h('Cell[py.cellind].ampalist.append(syn_)'
        )

                        # cellind=(re.search(r'\[(.*)\]', secnames).group(1))
                                    # print cellind, 'cellind'
                # place the post synapse.

                                    h('synlist.append(syn_)')
                                    h('gidn=0')
                                    h.gidn = int(data[0][1][3])
                                    h('Cell[py.cellind].gvpre.append(gidn)'
        )
                                    h('Cell[py.cellind].div.append(gidn)'
        )


                                    h('nc = pc.gid_connect(gidn, syn_)')

                        # connect the pre synapse to the post synapse.
                        # h('print coords.x[3], "this is the global identifier"')
                        # h('Cell[py.int(tail_target)].gvpre.append(gidn)')
                        # h('Cell[py.int(tail_target)].gvpost.append(py.int(py.gidcompare))')

                                    h('print nc," ", nc.srcgid()," ", gidn'
        )
                                    print 'the gids of the connected cells are: ', \
    i, ' ', gidn, '\n'

                                    # quit()

                                    h('nclist.append(nc)')

                                    h('nc.delay = 1+((py.r)/delay)')
                                    #Was this value
                                    #h('nc.delay = 1+((200+py.r)/delay)')
                                    #
                                    # h('nc.weight = r')

                                    h('nc.weight = py.r/iw')


COMM.Barrier()
'''
n_cell = 0
h('py.n_cell=cells.count')

# cell_ax = np.linspace(0.0, len(ie), len(ie))#sampling frequency

ie0 = gather(ie0)
ie1 = gather(ie1)
#ie1, and ie2 are already 2 column. Zipping them togethor would create 4 columns.
ie = zip(ie0, ie1)
COMM.Barrier()


my_icm = np.zeros_like(icm)

COMM.Reduce([icm, MPI.DOUBLE], [my_icm, MPI.DOUBLE], op=MPI.SUM,
            root=0)

my_ecm = np.zeros_like(ecm)

COMM.Reduce([ecm, MPI.DOUBLE], [my_ecm, MPI.DOUBLE], op=MPI.SUM,
            root=0)


COMM.Barrier()


h('for i=0,cells.count-1{ print cells.o(i).num_type }')

h('objref tgt')
h('srcid =0')


COMM.Barrier()
if RANK == 0:
    print 'source\ttarget\tsynapse\n'

# conn_m=[]



'''
for i in xrange(0, SIZE-1):
    if i == RANK:
        precell = []
        postcell = []
        nclist = 0
        h('py.nclist=py.int(nclist.count-1)')
        for j in xrange(0, int(nclist)):

            h('srcid=nclist.o(py.j).srcgid()')
            h('tgt = nclist.o(py.j).syn')
            srcid = 0
            tgt = 0

            h('print srcid, " ", tgt.cid')
            h('py.srcid=srcid')
            h('py.tgt=tgt.cid')

            precell.append(int(srcid))
            postcell.append(int(tgt))
        sconm = zip(precell, postcell)
        #sconm = np.transpose(np.matrix(sconm))
        if np.sum(np.array(sconm)) == 0:
            print 'no connections, quiting!'
            raise SystemExit(0)

sconm = COMM.gather(sconm, root=0)


COMM.Barrier()
'''
#sconm2 = []
#

'''
    for k in xrange(0, SIZE):
        sconm2 = sconm2 + sconm[k]  # Concatonate the lists so that they are all the size.

        print sconm2
        conm = [[0 for x in xrange(0, int(NCELL))] for x in xrange(0,
                int(NCELL))]
        print len(sconm2[:][:])

        for j in xrange(0, len(sconm2[:][:])):
            print sconm2[j][0], sconm2[j][1], sconm2[j]
            conm[sconm2[j][0]][sconm2[j][1]] = \
                conm[sconm2[j][0]][sconm2[j][1]] + 1
    print conm
    fin = 0

    # fin = fin + 1

    fig = plt.figure()
    fig.clf()

    im = plt.imshow(conm, interpolation='nearest')

    plt.autoscale(True)
    plt.colorbar(im)
    plt.xlabel('columns = targets')
    plt.ylabel('rows = sources')
    plt.title('Parallel Conection Matrix')
    plt.grid(True)

    sfin = str(NCELL) + ' ' + str(SIZE) \
        + 'Parallel_connection_matrix.png'
    fig.savefig(sfin)
    print sfin
    
'''
if RANK == 0:
#These methods are only failing when there is either no ecm or icm matrix. 
    fig = plt.figure()
    fig.clf()
    im = plt.imshow(my_ecm, interpolation='nearest')
    plt.autoscale(True)
    plt.colorbar(im)
    plt.xlabel('columns = targets')
    plt.ylabel('rows = sources')
    plt.title('Excitatory Conection Matrix')
    plt.grid(True)
    sfin = str(NCELL) + ' ' + str(SIZE) \
        + 'Excitatory_connection_matrix.png'
    fig.savefig(sfin)
    print sfin

    fig = plt.figure()
    fig.clf()
    im = plt.imshow(my_icm, interpolation='nearest')
    plt.autoscale(True)
    plt.colorbar(im)
    plt.xlabel('columns = targets')
    plt.ylabel('rows = sources')
    plt.title('Inhibitory Conection Matrix')
    plt.grid(True)
    sfin = str(NCELL) + ' ' + str(SIZE) \
        + 'Inhibitory_connection_matrix.png'
    fig.savefig(sfin)
    print sfin
#else:



if RANK == 0:
# only define and execute these functions on RANK 0


    def get_in(ma):
        old_row = 0
        row_index = 0
        old_j = 0
        for j in xrange(0, int(ma.shape[0])):

      # print ma[j,:]

            if sum(ma[j, :]) > old_row:  # old :
                old_row = sum(ma[j, :])
                row_index = j

            # print row_index, 'j= ', j, old_row, old_j

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

            # print column_index, 'i= ', i, old_column, old_i

                old_i = i

            # print ma

        return column_index


    outdegree = 0
    indegree = 0
    outdegreeg = 0
    indegreeg = 0


    def obtain_values(q, r):
        global outdegree, indegree, indegreeg, outdegreeg
        ma = q
        mg = r


        outdegree = get_out(ma)
        indegree = get_in(ma)
        outdegreeg = get_out(mg)
        indegreeg = get_in(mg)


        return

    #Call the function.
    obtain_values(np.array(ecm), np.array(ecm))
    #
    degrees = []
    degrees.append(outdegreeg)
    degrees.append(indegreeg)
    degrees.append(outdegree)
    degrees.append(indegree)
else:
    degrees = None
degrees = COMM.bcast(degrees, root=0)
print degrees, RANK

# It may not be necessary to remove this mechanism if pointers are probably set, as they are in rig7.hoc in the serial version. I remember this being a problem.

COMM.Barrier()

h('xopen("rigp.hoc")')  # Check for strength of external stimulation.
h('forall{ for(x,0){ uninsert xtra }}')  # Parallel Simulation causes seg fault wi

COMM.Barrier()
if RANK == 0:
    print 'cell numbers on the various hosts' 
    for i in xrange(0, SIZE):
        if i == RANK:
            print 'host', RANK, '\n'
            h('print "cells.count ", cells.count')
COMM.Barrier()
h('spikerecordg()')  # has to be called before the run only, otherwise it will delet
#I think I called this twice once in rigp. And that may have been the cause of the problem.
COMM.Barrier()


prun(tstop)

###After Run

#idveci=gather(idveci)
#tveci=gather(tveci)
#idvece=gather(idvece)
#tvece=gather(tvece)

h('spikeout()')

COMM.Barrier()

timec = []
ncell = 0.
h('py.ncell=cells.count-1')
for i in xrange(0, int(ncell)):
    timec.append(h.recvectors[int(i)].to_python())

   # print h.recvectors[int(i)].to_python()

timec = COMM.gather(timec, root=0)
COMM.Barrier()
if RANK == 0:
    timec2 = []
    for k in xrange(0, SIZE):
        timec2 = timec2 + timec[k]  # Concatonate the lists so that they are all the size.


    tc = np.array(timec2[int(0)])

   # tc = np.array(timec2[int(0)])

    N = len(tc)
    t = np.linspace(0., 0.025 * N, N)  # sampling frequency

    #fin += 1
    fig = plt.figure()
    fig.clf()
    plt.hold(True)
    for i in xrange(0, int(len(timec2))):
        tc = np.array(timec2[int(i)])
        tc = tc + i * 80
        plt.plot(t, tc, linewidth=1)

        plt.title('cells index membrane potential')
        plt.xlabel('ms')
        plt.ylabel('mV')

   # out1=downsample(time_courses[int(i)],oldrate=40000,newrate=200)
   # downsampled.append(out1)

    plt.hold(False)
    sfin = str(RANK) + ' ' + str(SIZE) + 'traces_all.png'
    fig.savefig(sfin)  # ,dpi='1000')

   # h('system("eog traces_all0.png")')

COMM.Barrier()
h('spikeout()')
h('vout()')

# Run first time after simulation. This should be on the end of the simulation code.


COMM.Barrier()
h('xopen("post_analysisp.hoc")')
h('objref vb[2]')
h('binsz=10//ms')

COMM.Barrier()

n_cell = 0
h('py.n_cell=idvec.size')
print 'completed'

COMM.Barrier()


tvec = gather(np.array(h.tvec.to_python()))
idvec = gather(np.array(h.idvec.to_python()))
COMM.Barrier()
#import brain_functions as bf

#bf.spk_plt(RANK, NCELL, SIZE)


#h('polvector = new Vector()')
#h('for i=0,cells.count-1{ polvector.x[i]=cells.o(i).polarity}')

#polvector = np.array(h.polvector.to_python(), dtype=np.float64)
#polvector = gather(polvector)





if RANK==0:




    #idvece=idvec[np.where(ie[1]==0)] 
    #tvece=tvec[np.where(ie[1]==0)]
    #idveci=idvec[np.where(ie[1]==1)] 
    #tveci=tvec[np.where(ie[1]==1)]
    '''
    fig = plt.figure()
    fig.clf()
    plt.ylabel('Cell number')
    plt.xlabel('spike time (ms)')
    plt.plot(tveci, idveci, 'b.')
    plt.plot(tvece, idvece, 'g.')
    sfin = 'spikes_colour_coded' + str(NCELL) + str(SIZE) + '.png'
    fig.savefig(sfin)
    '''
    fig = plt.figure()
    fig.clf()
    plt.plot(tvec, idvec, 'b.')
    #idvec and tvec are not the same length for some reason.
    sfin = 'spikesbw' + str(NCELL) + str(SIZE) + '.png'
    fig.savefig(sfin)


    '''
    tvece = []
    tveci = []
    isivec = []
    idveci = []
    isiveci = []
    idvece = []
    isivece = []

    
    i = 0
    for row in ie:
        i = int(i)
        for i in xrange(0, len(idvec)):
            if int(row[0]) == int(idvec[i]):  # if we have selected the idvec from the list

            # of cell polarities.
            # if the polarity is + add the corresponding
            # identity and time to the appropriate lists.

                if int(row[1]) == int(0):
                    print row[1], 'excitatory'
                    idvece.append(idvec[i])
                    tvece.append(tvec[i])
                if int(row[1]) == int(1):
                    print row[1], 'inhibitory'

                    idveci.append(idvec[i])
                    tveci.append(tvec[i])

                    # All the polarities are 0 Excitatory

                    # ie[i][0]=gidn
                    # ie[i][1]=polarity

    vecs1=zip(idveci,tveci)
    vecs2=zip(idvece,tvece)
    #Zipping an empty vector with a full one creates another empty vector.
    #vecs=zip(vecs1,vecs2)
    '''

    #print np.shape(tvece), np.shape(idvece)
    #print np.shape(tveci), np.shape(idveci)
    #print vecs[:]
    #print vecs2[:]
#    return vecs

COMM.Barrier()

#if RANK == 0:
#    print tvec
#    print idvec
#    vecs=spkplt2(tvec,idvec,ie)

h.xopen("tools.hoc")
execfile('analysisp3.py')
h('time_finish=pc.time()')
h('print time_finish')



