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

import re

#import random1 as r
#rand=r.random1


#from neuron import h

# Dont import LFPy
#Importing LFPy imports neuron. NEURON has to be imported after mpi4py
hostname=socket.gethostname()

# initialize the MPI interface

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

'''
pc = h.ParallelContext()
h('objref pc')
h.pc = pc


s = "mpi4py thinks I am %d of %d,\
 NEURON thinks I am %d of %d\n"
cw = MPI.COMM_WORLD
print s % (cw.rank, cw.size, pc.id(), pc.nhost())
h('time_start=pc.time()')
'''
NFILE = 3175#1175 # Maximum on this machine.
#Number of files to consider. This always ends up being a number more than the number of cells used. However NCELL and NFILE are of a similar magnitude. Some cell files have to filtered out, and cannot be used.


#h('iwp=1')
#h('ew=0.0185')
#h('delay=3000')
#h('iw=1000 // divisor/ or scaler, so if its bigger weight should be smaller. Was 1000 but this should be stronger. delay equals internal weight.')

#h('dc=-1')  # delete cell. This retains the cell, but makes all of its section diameter 0us, and it makes all projection weights to it 0.
#h('default_channels=0')
#h('plastic=1')
#h('get_dist=0')
#h('dlc=0')  # dilate centres.
#structural_plotting = 1

# I normaly use alphabet size of two symbols (0,1). There is no reason why this could not be increased. Using multilevel quantisation.

#prunenet = 0
#h('qd=0')  # quick and dirty for debugging
#h('lpy=1')  # declare lpy here, the value is not important, its only important that the declaration exists.

# h('if(qd==1){prunenet=900}')

#h('objref py')
#h('py = new PythonObject()')

###

#tstop = 1250  # 1500#550
#h('tstop1=1250')

#f = open('tstopfile', 'w')
#f.write(str(tstop))
#f.close()
###

from neuron import h
h('{load_file("stdgui.hoc")}')
#h('{load_file("nrngui.hoc")}')

'''

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
'''
# h('sprint(worm,"%s%s",workingdir,"/openworm/CNG version/")')

#h('objref py')
#h('py = new PythonObject()')

checkpoint_interval = 50000.
def gather(data):
    #h.load_file("stdlib.hoc")
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


#h('proc pprint(){ if (pc.id()==0){ print $s1 }}')


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



def prep_list(NFILE,allrows):
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


def spk_plt(RANK, NCELL, SIZE):
    idveci=[]
    tveci=[]
    idvece=[]
    tvece=[]
    idvecl= np.array(h.idvec.to_python(), dtype=np.int)
    tvecl= np.array(h.tvec.to_python(), dtype=np.float64)
    cnt=0
    h('objref target')

    for i in idvecl:
        #print i
        if h.pc.gid_exists(int(i)):
            h.target=h.pc.gid2cell(int(i))
            if h.target.polarity==1:
                idveci.append(int(i))
                tveci.append(tvecl[cnt])
                print tvecl[cnt]
            if h.target.polarity==0:
                idvece.append(int(i))
                tvece.append(tvecl[cnt])
        cnt+=1
 
    tveci = gather(np.array(tveci))
    idveci = gather(np.array(idveci))
    tvece = gather(np.array(tvece)) 
    idvece = gather(np.array(idvece))
    #COMM.Barrier()
    if RANK!=0: pass
    else:
        '''
        fig = plt.figure()
        fig.clf()
        plt.ylabel('Cell number')
        plt.xlabel('spike time (ms)')
        print tveci[:]
        print tvece[:]
        print idveci[:]
        print idvece[:]
        plt.plot(tveci, idveci, 'b.')
        plt.plot(tvece, idvece, 'g.')
        sfin = 'spikes_colour_coded' + str(NCELL) + str(SIZE) + '.png'
        fig.savefig(sfin)
    
        fig = plt.figure()
        fig.clf()
        plt.plot(tvec, idvec, 'b.')
        #idvec and tvec are not the same length for some reason.
        sfin = 'spikesbw' + str(NCELL) + str(SIZE) + '.png'
        fig.savefig(sfin)
        '''
    return 
#h('objref py')
#h('py = new PythonObject()')

def mbd(RANK, NCELL, SIZE, allrows2, gidvec, h,soff,boff,s1):
#make soff and boff optional parameters in python function.
    #pdb.set_trace()
    #h.load_file("stdlib.hoc")
    #global s1
    #global i, cnt, cnti
    #h('chdir(workingdir)')
    #os.chdir(os.getcwd() + '/main')
    h('py = new PythonObject()')
    h('cnt=0')
    ie0 = np.zeros((NCELL, 2))
    ie1 = np.zeros((NCELL, 2))

    i=0
    cnti = 0
    cnt=0
    s1=''
    storename=''
    h('cnt=pc.id')
       
    for i in range(RANK+soff, boff,SIZE):#NCELL-1, SIZE):  # 20 was int(len
        h.py.i=i
        s1 = allrows2[i]
        h.py.s1=s1

        storename = str(s1[3])  # //simply being inside a loop, may be the main problem
        if re.search('.swc', storename):
            h.cell = h.mkcell(storename)


            h('cell.geom_nseg()')
            h('cell.gid1=py.i')
            h('cell.gvpre.printf')
            h('cell.gvpost.printf')
            



            h.cell.nametype=str(s1[5])
            #h('cell.nametype=py.str(py.s1[5])')
            h.cell.num_type=int(s1[6])
            #h('cell.num_type=py.int(py.s1[6])')
            h('cell.population=py.str(py.s1[7])')
            #h('cell.reponame=py.str(py.storename)')
            h.cell.reponame=str(storename)
            h('cell.div.resize(py.int(py.NCELL))')
            h('cell.conv.resize(py.int(py.NCELL))')
            #h('cell.nametype=py.str(py.s1[5])')

            h('if(strcmp("pyramid",py.s1[5])==0){pyr_list.append(cell)}')
            h('if(strcmp(cell.population,"neocortex")==0){ if(strcmp("pyramid",py.s[5])==0){ cell.pyr() }}')


            h('if(strcmp(cell.population,"hippocampus")==0){ if(strcmp("pyramid",py.s[5])==0){ cell.pyr2() }}')

            #if 'hippocampus' in s1:
            #    h('cell.pyr2()')


            if 'interneuron' in s1:
                ie0[cnti] = i
                ie1[cnti] = 1
                cnti += 1
                h('cell.polarity=1')
            if 'pyramid' in s1:
                ie0[cnti] = i
                ie1[cnti] = 0
                cnti += 1
                h('cell.polarity=0')


            h('if(strcmp("interneuron",py.s1[5])==0){ inter_list.append(cell) }')
            h('if(strcmp("interneuron",py.s1[5])==0){ cell.basket()}')
            h('if(strcmp("interneuron",py.s1[5])==0){ print "inter neuron"}')

            h('if(strcmp("aspiny",py.s1[5])==0){ aspiny_list.append(cell) }')

            h('if(strcmp(cell.population,"hippocampus")==0){ hipp.append(cell)  }')

            h('if(strcmp(cell.population,"neocortex")==0){ neoc.append(cell) }')

            h('strdef cellposition')
            h('sprint(cellposition,"%s%d%s%d%s%d%s","cell.position(",py.float(py.s1[0]),",",py.float(py.s1[1]),",",py.float(py.s1[2]),")")')
            h('print cellposition')
            h('execute(cellposition)')
            
            #h.pc.set_gid2node(int(i),RANK)
            h('pc.set_gid2node(int(py.i), pc.id)')  # // associate gid i with this host

            h('cell.soma[0] nc =  new NetCon(&v(0.5), nil)')
            
            h('pc.cell(int(py.i), nc)')  # //
            #h.pc.cell(int(i),h.nc)

            h('cells.append(cell)')
            h('gidvec.append(py.i)')
            #h.gidvec.append(i)
            gidvec.append(i)
            print i
            h('print py.i, " py.i", cnt, " cnt"')
            h('cnt+=pc.nhost')

            cnt += 1
            #strangley i can update in this environment but not in wirecells.
    #destroy the allrows list to free up memory, then return the empty list.
    #I can't seem to destroy allrows here without wrecking something in rigp later.
    #allrows=[]
    return (h.cells, allrows2, ie0, ie1)

def mb(RANK, NCELL, SIZE, allrows2, gidvec, h,s1):
#make soff and boff optional parameters in python function.
    #pdb.set_trace()
    #h.load_file("stdlib.hoc")
    #global s1
    #global i, cnt, cnti
    #h('chdir(workingdir)')
    #os.chdir(os.getcwd() + '/main')
    h('py = new PythonObject()')
    h('cnt=0')
    ie0 = np.zeros((NCELL, 2))
    ie1 = np.zeros((NCELL, 2))

    i=0
    cnti = 0
    cnt=0
    s1=''
    storename=''
    h('cnt=pc.id')
       
    for i in range(RANK, NCELL,SIZE):#NCELL-1, SIZE):  # 20 was int(len
        h.py.i=i
        s1 = allrows2[i]
        h.py.s1=s1

        storename = str(s1[3])  # //simply being inside a loop, may be the main problem
        if re.search('.swc', storename):
            h.cell = h.mkcell(storename)


            h('cell.geom_nseg()')
            #h.cell.geom_nseg()
            h('cell.gid1=py.i')
            h('cell.gvpre.printf')
            h('cell.gvpost.printf')
            



            h.cell.nametype=str(s1[5])
            #h('cell.nametype=py.str(py.s1[5])')
            h.cell.num_type=int(s1[6])
            #h('cell.num_type=py.int(py.s1[6])')
            h('cell.population=py.str(py.s1[7])')
            #h('cell.reponame=py.str(py.storename)')
            h.cell.reponame=str(storename)
            h('cell.div.resize(py.int(py.NCELL))')
            h('cell.conv.resize(py.int(py.NCELL))')
            #h('cell.nametype=py.str(py.s1[5])')

            h('if(strcmp("pyramid",py.s1[5])==0){pyr_list.append(cell)}')

            h('if(strcmp(cell.population,"neocortex")==0){ if(strcmp("pyramid",py.s[5])==0){ cell.pyr() }}')



            h('if(strcmp(cell.population,"hippocampus")==0){ if(strcmp("pyramid",py.s[5])==0){ cell.pyr2() }}')

            #if 'hippocampus' in s1:
            #    h('cell.pyr2()')


            if 'interneuron' in s1:
                ie0[cnti] = i
                ie1[cnti] = 1
                cnti += 1
                h.cell.basket()
                h('cell.polarity=1')
            if 'pyramid' in s1:
                ie0[cnti] = i
                ie1[cnti] = 0
                cnti += 1
                h.cell.pyr()
                h('cell.polarity=0')
            #h.cell.basket()
            #h('cell.polarity=1')


            h('if(strcmp("interneuron",py.s1[5])==0){ inter_list.append(cell) }')
            h('if(strcmp("interneuron",py.s1[5])==0){ cell.basket()}')
            h('if(strcmp("interneuron",py.s1[5])==0){ print "inter neuron"}')

            h('if(strcmp("aspiny",py.s1[5])==0){ aspiny_list.append(cell) }')

            h('if(strcmp(cell.population,"hippocampus")==0){ hipp.append(cell)  }')

            h('if(strcmp(cell.population,"neocortex")==0){ neoc.append(cell) }')

            h('strdef cellposition')
            h('sprint(cellposition,"%s%d%s%d%s%d%s","cell.position(",py.float(py.s1[0]),",",py.float(py.s1[1]),",",py.float(py.s1[2]),")")')
            h('print cellposition')
            h('execute(cellposition)')
            
            #h.pc.set_gid2node(int(i),RANK)
            h('pc.set_gid2node(int(py.i), pc.id)')  # // associate gid i with this host

            h('cell.soma[0] nc =  new NetCon(&v(0.5), nil)')
            
            h('pc.cell(int(py.i), nc)')  # //
            #h.pc.cell(int(i),h.nc)

            h('cells.append(cell)')
            h('gidvec.append(py.i)')
            #h.gidvec.append(i)
            gidvec.append(i)
            print i
            h('print py.i, " py.i", cnt, " cnt"')
            h('cnt+=pc.nhost')

            cnt += 1
            #strangley i can update in this environment but not in wirecells.
    #destroy the allrows list to free up memory, then return the empty list.
    #I can't seem to destroy allrows here without wrecking something in rigp later.
    #allrows=[]
    #pdb.set_trace()
    return (h.cells, allrows2, ie0, ie1)


import math
r=0.
def wirecells(RANK,NCELL,SIZE,h,icm,ecm):
    pdb.set_trace()

    #h('objref py')
    #h('py = new PythonObject()')
    #global s,j,i,test,test2,r,gidn
    gidn=0
    s=0
    j=0
    h.py.j=0
    i=0
    test=0
    test2=0
    gidcompare = ''
    r=0
    secnames = ''# sec.name()
    cellind =0 #int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.

    polarity = 0
    #z=0
    y=0

    for s in xrange(0, SIZE):#Was full SIZE, not SIZE-1
        for j in xrange(0,NCELL):
            print s, j
            if RANK == s:
                coordlist = [[0 for x in xrange(0, 2)] for x in xrange(0,2)]

                test = 0

                print h.py.j, 'is this variable updating?'
                h.py.j=j
                h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
                #test=int(h.pc.gid_exists(int(j)))
                if test != 0:
                    
                    h('source=pc.gid2cell(py.int(py.j))')
                    #h.source=h.pc.gid2cell(int(j))
                    for sec in h.source.spk_trig_ls:
                        for seg in sec:
                            get_cox = str('coords.x[0]=x_xtra('+ str(seg.x) + ')')
                            h(get_cox)


                            get_coy = str('coords.x[1]=y_xtra('+ str(seg.x) + ')')
                            h(get_coy)
                            get_coz = str('coords.x[2]=z_xtra('+ str(seg.x) + ')')
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
                    h.py.i=i
                    gidn = int(data[0][1][3])
                    if i != int(gidn):  # if the gids are not the same
                        #test2 = 0
                        #test2=int(h.pc.gid_exists(int(i)))
                        h('py.test2=py.int(pc.gid_exists(py.int(py.i)))')
                        if test2 != 0:
                            h('target=pc.gid2cell(py.int(py.i))')
                            #h.target=h.pc.gid2cell(int(i))
                            for sec in h.target.spk_rx_ls:
                                for seg in sec:
                                    h(str('coords2.x[2]=') + str('z_xtra(')+ str(seg.x) + ')')
                                    h(str('coords2.x[1]=') + str('y_xtra(')+ str(seg.x) + ')')
                                    h(str('coords2.x[0]=') + str('x_xtra(')+ str(seg.x) + ')')

                                    h('coordsx=0.0')
                                    h.coordsx = data[0][1][0]
                                    h('coordsy=0.0')
                                    h.coordsy = data[0][1][1]
                                    h('coordsz=0.0')
                                    h.coordsz = data[0][1][2]


                                    
                                    coordsx = float(data[0][1][0])
                                    coordsy = float(data[0][1][1])
                                    coordsz = float(data[0][1][2])
                                    type(coordsx)
                                    r = 0.
                                    r=math.sqrt((h.coords2.x[0] - coordsx)**2+(h.coords2.x[1] - coordsy)**2+(h.coords2.x[2] - coordsz)**2)
                                    #h('r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)')
                                    #testr=0
                                    #if 
                                    #r = float(h.r)
                                    #print r
                                    if r < 15:  
                                        #print 'got here'
                                        gidcompare = ''

                                        secnames = sec.name()
                                        cellind = int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.

                                        polarity = 0
                                   

                                        #cellind is a cell index, that is relative to the host. So the identifier repeats on different hosts.
                                        #gidn is a global identifier. These numbers are not repeated on different hosts.
                                        polarity=int(h.Cell[int(cellind)].polarity)
                                        h('objref syn_')
                                        print polarity       # (
                                        if int(polarity) == int(1):
                                            post_syn = secnames + ' ' + 'syn_ = new GABAa(' + str(seg.x) + ')'
                                            icm[i][gidn] = icm[i][gidn] + 1
                                            h('gidn=0')
                                            h.gidn = int(data[0][1][3])
                                            h(post_syn)
                                            h('print syn_')
                                            h.syn_.cid=int(i)
                                            h('print syn_.cid')

                                        #h('syn_.cid=py.int(py.i)')  # This makes the gid a property of the synapse.
                                            h.Cell[cellind].gabalist.append(h.syn_)

                                        else:

                                            post_syn = secnames + ' ' + 'syn_ = new AMPA(' + str(seg.x) + ')'
                                            ecm[i][gidn] = ecm[i][gidn] + 1
                                            h('gidn=0')
                                            h.gidn = int(data[0][1][3])
                                            h(post_syn)
                                            h('print syn_')
                                            h.syn_.cid=int(i)
                                            h('print syn_.cid')

                                        #h('syn_.cid=py.int(py.i)')  # This makes the gid a property of the synapse.
                                            h.Cell[cellind].ampalist.append(h.syn_)


                                        h.synlist.append(h.syn_) 
                                        #h('synlist.append(syn_)')
                                        h('gidn=0')
                                        h.gidn = int(data[0][1][3])
                                        h.Cell[cellind].div.append(h.gidn)
                                        h.Cell[cellind].gvpre.append(h.gidn)


                                        h('nc = pc.gid_connect(int(gidn), syn_)')
                                        #h.nc=h.pc.gid_connect(int(gidn),h.syn_)
             
                                        print h.nc,' ',h.nc.srcgid(),' ',h.gidn
                                        h('print nc," ", nc.srcgid()," ", gidn')
                                        print h.nc,' ', h.nc.srcgid(), ' ',gidn
                                        print 'the gids of the connected cells are: ', i, ' ', gidn, '\n'
                                        print int(h.gidn), gidn
                                        #y+=1
                                        #h('print py.y, "is this variable updating?"')
                                        #h('print py.gidn, "is this variable updating?"')
                                        #print h.py.j
                                        h.nc.delay=1+float(r/3000)
                                        #h('nc.delay = 1+((py.r)/delay)')
                                        #h('nc.weight = py.r/iw')
                                        #print h.nc.weight
                                        #print r/h.iw
                                        h.nc.weight[0]=float(r/h.iw)#zeroth element refers to the synaptic weight.

                                        h('nclist.append(nc)')



    COMM.Barrier()

    return h.nclist, ecm, icm


#i=0
#j=0
#s=0
#gidn=0
def wirecells2(RANK,NCELL,SIZE,h,icm,ecm):
#def wirecells(RANK,NCELL,SIZE,h,icm,ecm):
    pc=h.ParallelContext()
    global s,j,i,test,test2,r
    s=0
    j=0
    i=0
    gidcompare = ''
    '''
    h('py.j=0')
    h('py.i=0')
    h('py.s=0')
    h('py.test=0')
    h('py.test2=0')
    h('py.r=0')
    h('py.gidcompare=""')
    '''

    secnames = ''# sec.name()
    cellind =0 #int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.

    polarity = 0
    
    h.py.j=0

    for s in xrange(0, SIZE):#Was full SIZE, not SIZE-1
        for j in xrange(0,NCELL):
            h.py.j=int(j)
            #print s, j
            if RANK == s:
                coordlist = [[0 for x in xrange(0, 2)] for x in xrange(0,
                         2)]

                test = 0
                #h('print py.j, py.test')
                test=int(pc.gid_exists(j))
                print 'test ', test
                #h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
                if test != 0:
                    h.source=pc.gid2cell(j)
                    print h.source
                    #h('source=pc.gid2cell(py.int(py.j))')
                    #h('print source')
                    h.source=h.pc.gid2cell(int(j))
                    h('print source')

                    for sec in h.source.spk_trig_ls:
                        for seg in sec:
                            #print sec, seg
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
                #pdb.set_trace()
 
                #h.py.i=0
                for i in xrange(0, NCELL-1):
                    h.py.i=i
                    gidn = int(data[0][1][3])
                    if i != int(gidn):  # if the gids are not the same
                        test2 = 0
                        test2=pc.gid_exists(i)
                        #h('py.test2=py.int(pc.gid_exists(py.int(py.i)))')
                        print test2
                        if test2 != 0:
                            h.target=pc.gid2cell(i)
                            print h.target
                            h('target=pc.gid2cell(py.int(py.i))')
                            h('print target')
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
                      
                                    coordsx = float(data[0][1][0])
                                    coordsy = float(data[0][1][1])
                                    coordsz = float(data[0][1][2])
                                    r = 0.
                                    r=math.sqrt((h.coords2.x[0] - coordsx)**2+(h.coords2.x[1] - coordsy)**2+(h.coords2.x[2] - coordsz)**2)

                                    #h('py.r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)')
                                    r = float(r)
                                    if r < 10:  
                                        print r, 'this is not hopefuly wiring everything to everything'
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
                                        print post_syn
                                        h('print syn_')
                                        #h.syn_.cid=i
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


                                        h.py.r=r
                                        h('nc.delay = 1+((py.r)/delay)')
                                        h('nc.weight = py.r/iw')
                                        h('nclist.append(nc)')


    COMM.Barrier()

    return (h.nclist, ecm, icm)
#h.nclist=wirecells(RANK,NCELL,SIZE,h,icm,ecm)


def wirecells3(RANK,NCELL,SIZE,h,icm,ecm):
#def wirecells(RANK,NCELL,SIZE,h,icm,ecm):
    pc=h.ParallelContext()
    global s,j,i,test,test2,r
    s=0
    j=0
    i=0
    gidcompare = ''
    '''
    h('py.j=0')
    h('py.i=0')
    h('py.s=0')
    h('py.test=0')
    h('py.test2=0')
    h('py.r=0')
    h('py.gidcompare=""')
    '''

    secnames = ''# sec.name()
    cellind =0 #int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.

    polarity = 0
    
    h.py.j=0


    "crash here please"
    for s in xrange(0, SIZE):#Was full SIZE, not SIZE-1
        for j in xrange(0,NCELL):
            h.py.j=int(j)
            #print s, j
            if RANK == s:
                coordlist = [[0 for x in xrange(0, 2)] for x in xrange(0,
                         2)]

                test = 0
                #h('print py.j, py.test')
                test=int(pc.gid_exists(j))
                print 'test ', test
                #h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
                if test != 0:
                    h.source=pc.gid2cell(j)
                    print h.source
                    #h('source=pc.gid2cell(py.int(py.j))')
                    #h('print source')
                    h.source=h.pc.gid2cell(int(j))
                    h('print source')

                    for sec in h.source.spk_trig_ls:
                        for seg in sec:
                            #print sec, seg
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
                #pdb.set_trace()
 
                #h.py.i=0
                for i in xrange(0, NCELL-1):
                    h.py.i=i
                    gidn = int(data[0][1][3])
                    if i != int(gidn):  # if the gids are not the same
                        test2 = 0
                        test2=pc.gid_exists(i)
                        #h('py.test2=py.int(pc.gid_exists(py.int(py.i)))')
                        print test2
                        if test2 != 0:
                            h.target=pc.gid2cell(i)
                            print h.target
                            h('target=pc.gid2cell(py.int(py.i))')
                            h('print target')
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
                      
                                    coordsx = float(data[0][1][0])
                                    coordsy = float(data[0][1][1])
                                    coordsz = float(data[0][1][2])
                                    r = 0.
                                    r=math.sqrt((h.coords2.x[0] - coordsx)**2+(h.coords2.x[1] - coordsy)**2+(h.coords2.x[2] - coordsz)**2)

                                    #h('py.r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)')
                                    r = float(r)
                                    if r < 10:  
                                        print r, 'this is not hopefuly wiring everything to everything'
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
                                        print post_syn
                                        h('print syn_')
                                        #h.syn_.cid=i
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


                                        h.py.r=r
                                        h('nc.delay = 1+((py.r)/delay)')
                                        h('nc.weight = py.r/iw')
                                        h('nclist.append(nc)')


    COMM.Barrier()

    return (h.nclist, ecm, icm)



h('objref storval')
h('storval = new Vector()')
h('objref vb[2]')
h('vb[0]=new Vector()')
h('vb[1]=new Vector()')
h('objref vec')
h('vec=new Vector()')
#h('objref vec1, vec2 ')
#h('vec1=new Vector()')
#h('vec2=new Vector()')

##
# I think the whole problem is that Python resets tstop to 5ms by default.
# Thus breaking my approach to histogramming.
##

fn = open('tstopfile')
tstopl = [line.strip() for line in open('tstopfile', 'r')]
tstop = int(tstopl[0])

h('objref work1')
h('objref work2')
h('work1 = new Vector()')
h('work2 = new Vector()')
h('objref source')
COMM.Barrier()

histogram=[]
ainfplt=[]



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




def t_e(SIZE, NCELL, RANK, trentm,np):
    #The HOC variable does not need to passed to the function. Its all always visible. Don't ask me why.
    #if these variables are not declared as global, they simply function to work.
    #h('objref py')
    #h('py = new PythonObject()')

    global test,testn,test3
    global gidn
    global i, ent
    global k, j, i
    gidn=0
    ent=0
    k=0
    j=0
    i=0
    histogram=[]
    lists=[]
    gidn=0
    test=0
    testn=0
    ent=0
    COMM.Barrier()
    h('objref listing')
    h('listing = new List()')
    h('objref storval')
    h('storval = new Vector()')
    h('objref vb[2]')
    h('vb[0]=new Vector()')
    h('vb[1]=new Vector()')
    h('objref vec')
    h('vec=new Vector()')
    h('objref source')
    h('j=0')

#
# I think the unusual parts of the transfer entropy matrix
# are caused by SIZE and NCELL being incremented one too many times.
# SIZE should only be incremented to be s...SIZE-1
#
    pdb.set_trace()
    for k in xrange(0, SIZE):
        for j in xrange(0, int(NCELL-1)):
            if RANK == k:
                histogram = []
                #print 'j, s', j, ' ', s
                test = 0
                h.py.j=j
                h.j=j
                #h('py.test=py.int(pc.gid_exists(py.j))')
                #h.py.test=h.pc.gid_exists(j)
                test=h.pc.gid_exists(j)
                if test != 0:
                    h('source=pc.gid2cell(j)')
                    print j
                    h('print source, "source"')
                    h('print j')
                    h('vec=source.recvec1')

                    h('vec.x[0]=py.j')

                    histogram = h.source.recvec1.to_python()
                if test == 0:
                    histogram = None
                    h('vec.fill(0)')
            else:

                histogram = None
                h('vec.fill(0)')
            histogram = COMM.bcast(histogram, root=k)  # ie root = rank
            #assert isinstance(histogram, list)
            h.py.k=k
            h('pc.broadcast(vec,py.int(py.k))')
            #h.py.gidn=h.vec.x[0]
            h('py.gidn=int(vec.x[0])')
            h('vec.x[0]=0')
            gidn = int(gidn)
            h('py.testn=int(vec.sum())')
            testn = h.vec.sum()
            #print 'testn = ', testn, ' histogram= ', histogram


            if int(testn) != 0:
                print 'cell number ', gidn
                #h('vec.printf')
                for i in xrange(0, int(NCELL - 1)):
                    h.py.i=i
                    h.i=i
                    test3 = 0
                    #print i, ' i '
                    h('py.test3=py.int(pc.gid_exists(py.int(py.i)))')
                    #test3=h.pc.gid_exists(i)
                    if test3 != 0:
                        h('target=pc.gid2cell(py.int(py.i))')
                        h('print target, " target"')
                        h('print i')

                      # # Here is the problem. Histograming creates a vector with only 1, 0 element. The problem is in the Vector data structure, and its inability to get updated in this context.
                      # # I will need to work around it somehow. Possibly by using different data structures
                      # Its only the binning that does not work here, so why not do it somewhere else, when then the spike bins are first attributed to the cells.
                        #print 'got here'
                        h('storval=normte(target.recvec1,vec,20)')
                        #h('print storval')
                        #h('storval.printf')
                        ent = 0.0
                        h('py.ent=py.float(storval.x[2])')
                        #h('listing.append(py.float(storval.x[2]))')
                        trentm[i][gidn] = ent
                        
                        trentm[i][gidn] =float(h.storval.x[2])
                        #lists.append(float(h.storval.x[2]))
                        #print i, j, k
                        #h('target.recvec1.printf')
                        #h('vec.printf')
                        
    #COMM.Barrier()
    '''
    h('print "listingcnt", listing.count')
    print trentm[:]
    print np.sum(trentm), 'np.sum'
    print lists[:]
    print k
    h('print "listingcnt", listing.count')
    print trentm[:]
    print np.sum(trentm), 'np.sum'
    print lists[:]
    print k
    '''
    #print i, j, k
    #h('target.recvec1.printf')
    #h('vec.printf')

    return trentm           
'''
trentm=t_e(SIZE, NCELL, RANK, trentm,np)

my_trentm = np.zeros_like(trentm)

COMM.Reduce([trentm, MPI.DOUBLE], [my_trentm, MPI.DOUBLE], op=MPI.SUM,
            root=0)

COMM.Barrier()
#pdb.set_trace()
############

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

############


COMM.Barrier()
'''

def t_e2(SIZE, NCELL, RANK, histogram, trentm, visited, h):
 
    
    #I think the unusual parts of the transfer entropy matrix
    #are caused by SIZE and NCELL being incremented one too many times.
    #SIZE should only be incremented to be s...SIZE-1
    
    
    
    global s,j,i,test,testn,r
    s=0
    j=0
    i=0

    global gidn
    gidn=0
    COMM.Barrier()



    for s in xrange(0, SIZE-1):
        for j in xrange(0, int(NCELL-1)):
            h.py.j=j
            if RANK == s:
                histogram = []
                #print 'j, s', j, ' ', s
                test = 0
                
                h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
            
                if test != 0:
                    h('source=pc.gid2cell(py.int(py.j))')
                    h('vec=source.recvec1')
                    h('vec.x[0]=py.j')
                    histogram = h.source.spk_trainb.to_python()
                    #h('source=pc.gid2cell(py.int(py.j))')
                    #h('vec=source.spk_trainb')
                    #h('vec.x[0]=py.j')
                    #histogram = h.source.spk_trainb.to_python()
                if test == 0:
                    histogram = None
                    h('vec.fill(0)')
            else:

                histogram = None
                h('vec.fill(0)')
            histogram = COMM.bcast(histogram, root=s)  # ie root = rank
            h('pc.broadcast(vec,py.int(py.s))')

            h('py.gidn=int(vec.x[0])')
            h('vec.x[0]=0')
            gidn = int(gidn)
            h('py.testn=int(vec.sum())')
            testn = h.vec.sum()
            print 'testn = ', testn, ' histogram= ', histogram


            if int(testn) != 0:
                print 'cell number ', gidn
                h('vec.printf')
                for i in xrange(0, int(NCELL - 1)):
                    h.py.i=i
                    test3 = 0
                    print i, ' i '
                    h('py.test3=py.int(pc.gid_exists(py.int(py.i)))')
                    if test3 != 0:
                        visited[i][gidn] = visited[i][gidn] + 1

                        h('target=pc.gid2cell(py.int(py.i))')
                        h('target.spk_train.printf')

                        h('storval=normte(vec,target.spk_trainb,20)')
                        #ent = 0.0
                        #h('py.ent=py.float(storval.x[2])')
                        #trentm2[i][gidn] = ent
                        trentm[i][gidn] =float(h.storval.x[2])


                        #h('target=pc.gid2cell(py.int(py.i))')

                      # # Here is the problem. Histograming creates a vector with only 1, 0 element. The problem is in the Vector data structure, and its inability to get updated in this context.
                      # # I will need to work around it somehow. Possibly by using different data structures
                      # Its only the binning that does not work here, so why not do it somewhere else, when then the spike bins are first attributed to the cells.

                        #h('storval=normte(target.spk_trainb,vec,20)')
                        #ent = 0.0
                        #h('py.ent=py.float(storval.x[2])')
                        #trentm[i][gidn] =float(h.storval.x[2])
                        #trentm[i][gidn] = ent
    COMM.Barrier()
    print trentm[:]
    print np.sum(trentm), np.sum
    return (trentm, visited)    


'''
trentm2,visited=t_e2(SIZE, NCELL, RANK, histogram, trentm2, visited, h)
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

    im = plt.imshow(my_trentm2, interpolation='nearest')

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
