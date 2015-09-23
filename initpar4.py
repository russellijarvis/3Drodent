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
import unittest

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
PROJECTROOT=os.getcwd()

NFILE = 3175#1175 # Maximum on this machine.
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

#prunenet = 0
h('qd=0')  # quick and dirty for debugging
h('lpy=1')  # declare lpy here, the value is not important, its only important that the declaration exists.

# h('if(qd==1){prunenet=900}')

h('objref py')
h('py = new PythonObject()')

###

tstop = 100  # 1500#550
h('tstop1=500')

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
soff=1315
boff=1450
#NCELL=boff-soff
#1350,1450)
#NFILE = 3175#1175 # Maximum on this machine.
#NCELL=int(np.shape(allrows2)[0])
NCELL=185
# Should be max_NCELL=int(np.shape(allrows2)[0])
# Keep NFILE the same, change 
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

#import brain_functions as bf
os.chdir(PROJECTROOT)
#execfile('brain_functions.py')

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
os.chdir(os.getcwd() + '/main')
#os.chdir(PROJECTROOT)

import brain_functions as bf
#For a variable to be accessible to the HOC object it has to be declared out here. 
s1=''

#h.cells, allrows2, ie0, ie1 =bf.mbd(RANK, NCELL, SIZE, allrows2, gidvec,h,1250,1350,s1)
h.cells, allrows2, ie0, ie1 =bf.mb(RANK, NCELL, SIZE, allrows2, gidvec,h,s1)
allrows2=[]
allrows=[]
h('for i=0,cells.count-1{print cells.o(i).polarity }')
#raise SystemExit(0)
allrows2=[]
#h.cells, allrows2, ie0, ie1=mb(RANK, NCELL, SIZE, allrows2, gidvec,h,soff,boff)#1350,1450)
#allrows2=[]#destroy the memory expensive allrows2 list.
#raise SystemExit(0)

h('forall{ for(x,0){ insert xtra }}')
h('forall{ for(x,0){ insert extracellular}}')

os.chdir(PROJECTROOT)

h('xopen("interpxyz.hoc")')
h('grindaway()')

# h('forall{ for(x,0){ print x_xtra }}')

#h('system("pwd")')


#just changed 2015

#h('xopen("seclists.hoc")')
#raise SystemExit(0)

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



#import brain_functions as bf
#z=0
#NCELL=boff

from utils import Utils
nclist, ecm, icm=utils.wirecells_s()#Wire cells on same host.
nclist, ecm, icm=utils.wirecells3()#wire cells on different hosts.

#h.nclist,ecm,icm=bf.wirecells3(RANK,NCELL,SIZE,h,icm,ecm)

#h.nclist,ecm,icm=bf.wirecells2(RANK,NCELL,SIZE,h,icm,ecm)
h('print nclist')
#raise SystemExit(0)
print np.sum(ecm)
print np.sum(icm)

os.chdir(PROJECTROOT)
#execfile('brain_functions.py')

#h('chdir(workingdir)')
COMM.Barrier()
n_cell = 0
#h('py.n_cell=cells.count')

# cell_ax = np.linspace(0.0, len(ie), len(ie))#sampling frequency

#ie0 = gather(ie0)
#ie1 = gather(ie1)
#ie1, and ie2 are already 2 column. Zipping them togethor would create 4 columns.
#ie = zip(ie0, ie1)
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

            #precell.append(int(srcid))
            #postcell.append(int(tgt))
        #sconm = zip(precell, postcell)
        #sconm = np.transpose(np.matrix(sconm))
        #if np.sum(np.array(sconm)) == 0:
        if nclist==0:
            print 'no connections, quiting!'
            raise SystemExit(0)



COMM.Barrier()
if RANK == 0:
 
    fig = plt.figure()
    fig.clf()
    im = plt.imshow(my_ecm, interpolation='nearest')
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
    plt.colorbar(im)
    plt.xlabel('columns = targets')
    plt.ylabel('rows = sources')
    plt.title('Inhibitory Conection Matrix')
    plt.grid(True)
    sfin = str(NCELL) + ' ' + str(SIZE) \
        + 'Inhibitory_connection_matrix.png'
    fig.savefig(sfin)
    print sfin



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
    obtain_values(np.array(my_ecm), np.array(my_icm))
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
#raise SystemExit(0)

h('xopen("rigp.hoc")')  # Check for strength of external stimulation.
h('forall{ for(x,0){ uninsert xtra }}')  # Parallel Simulation causes seg fault wi
h.spikerecordg()
COMM.Barrier()
if RANK == 0:
    print 'cell numbers on the various hosts' 
    for i in xrange(0, SIZE):
        if i == RANK:
            print 'host', RANK, '\n'
            h('print "cells.count ", cells.count')
COMM.Barrier()
h.spikerecordg()  # has to be called before the run only, otherwise it will delet
#I think I called this twice once in rigp. And that may have been the cause of the problem.
COMM.Barrier()


def progress(pinvl, swlast):
#From http://senselab.med.yale.edu/ModelDB/ShowModel.asp?model=151681&file=\bulb3d\util.py
  sw = h.startsw()
  print "t=%g wall interval %g"% (h.t, sw-swlast)
  h.cvode.event(h.t+pinvl, (progress, (pinvl , sw)))

def show_progress(invl):
#From http://senselab.med.yale.edu/ModelDB/ShowModel.asp?model=151681&file=\bulb3d\util.py

  global fih
  if RANK == 0:
    fih = h.FInitializeHandler(2, (progress, (invl, h.startsw())))

show_progress(200)
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







if RANK==0:
    fig = plt.figure()
    fig.clf()
    plt.plot(tvec, idvec, 'b.')
    sfin = 'spikesbw' + str(NCELL) + str(SIZE) + '.png'
    fig.savefig(sfin)

COMM.Barrier()


h.xopen('tools.hoc')
h.xopen('analysis2.hoc')

execfile('analysisp3.py')
h('time_finish=pc.time()')
h('print time_finish')


#import analysisp3 as ap

'''
trentm = np.zeros((NCELL, NCELL))
     
trentm=bf.t_e(SIZE, NCELL, RANK, trentm,np)

my_trentm = np.zeros_like(trentm)

COMM.Reduce([trentm, MPI.DOUBLE], [my_trentm, MPI.DOUBLE], op=MPI.SUM,
            root=0)

COMM.Barrier()
print 'np.sum(trentm)= ', np.sum(trentm)

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



COMM.Barrier()


trentm2 = np.zeros((NCELL, NCELL))
visited = np.zeros((NCELL, NCELL))


trentm2,visited=bf.t_e2(SIZE, NCELL, RANK, histogram, trentm2, visited, h)


#trentm2,visited=t_e2(SIZE, NCELL, RANK, histogram, trentm2, visited, h)
my_trentm2 = np.zeros_like(trentm2)
COMM.Reduce([trentm2, MPI.DOUBLE], [my_trentm2, MPI.DOUBLE],
            op=MPI.SUM, root=0)

my_visited = np.zeros_like(visited)
COMM.Reduce([visited, MPI.DOUBLE], [my_visited, MPI.DOUBLE],
            op=MPI.SUM, root=0)


COMM.Barrier()

if RANK == 0:
    
    fig = plt.figure()
    fig.clf()

    im = plt.imshow(my_visited, interpolation='nearest')

    plt.autoscale(True)
    plt.colorbar(im)
    plt.xlabel('columns = targets')
    plt.ylabel('rows = sources')
    plt.title('visited')
    plt.grid(True)

    sfin = str(SIZE) + str(NCELL) + 'Visited.png'
    fig.savefig(sfin)
    fig = plt.figure()
    fig.clf()
    
    im = plt.imshow(my_trentm, interpolation='nearest')

    plt.autoscale(True)
    plt.colorbar(im)
    plt.xlabel('columns = targets')
    plt.ylabel('rows = sources')
    plt.title('Transfer Entropy Matrix Two')
    plt.grid(True)

    sfin = str(SIZE) + ' ' + str(NCELL) + 'Transfer_Entropy_matrix2.png'
    fig.savefig(sfin)

COMM.Barrier()

'''


