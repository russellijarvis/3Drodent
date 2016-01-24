#!/usr/bin/python
# -*- coding: utf-8 -*-

#from bsmart import granger  # Load the Granger calculation tool
import numpy as np
from matplotlib.colors import LogNorm
import numpy as np
from mpi4py import MPI
import pdb


from neuron import h

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
#h('time_start=pc.time()')

#execfile('analysisp3.py')
#trentm2 = np.zeros((NCELL, NCELL))
#visited = np.zeros((NCELL, NCELL))
#msgcv = np.zeros((NCELL, NCELL))

h.xopen('analysis2.hoc')

# h('xopen("isi.hoc")')

#surface = np.empty((NCELL, NCELL, NCELL))
#surfaceg = np.empty((NCELL, NCELL, NCELL))

#from pyhoc import downsample
h('storval = new Vector()')
h('objref vb[2]')
h('vb[0]=new Vector()')
h('vb[1]=new Vector()')
h('objref vec')
h('vec=new Vector()')

##
# I think the whole problem is that Python resets tstop to 5ms by default.
# Thus breaking my approach to histogramming.
##

h('objref work1')
h('objref work2')
h('work1 = new Vector()')
h('work2 = new Vector()')
h('objref source')
COMM.Barrier()



import numpy as np
import unittest
#global NCELL
#NCELL= int(np.shape(allmemb)[0])

class transfer():
    def init(self):       
        self.test
        self.testn
        self.test3
        self.gidn
        self.ent
        self.k
        self.j
        self.i
        
    def t_e(SIZE, NCELL, RANK, trentm,np):
    ''' 
        
        
    '''    
 
 
 
        self.test=0
        self.testn=0
        self.test3=0
        self.gidn=0
        self.ent=0
        self.k=0
        self.j=0
        self.i=0
        histogram=[]
        lists=[]
        h('objref listing')
        h('listing = new List()')
        h('storval = new Vector()')
        h('objref vb[2]')
        h('vb[0]=new Vector()')
        h('vb[1]=new Vector()')
        h('objref vec')
        h('vec=new Vector()')
        COMM.Barrier()
    
        #pdb.set_trace()#PDB can be used inside functions no worries.
      
        for k in xrange(0, SIZE):
            for j in xrange(0, int(NCELL-1)):
                if RANK == k:
                    histogram = []
                    test = 0
                    #h('print py.j, " py.j"')
                    test=int(h.pc.gid_exists(int(j)))
                    if test != 0:
                        h.source=h.pc.gid2cell(int(j))
                        h.vec=h.source.recvec1
                        h('print source, " source"')
                        h.vec.x[0]=j
                        
                        histogram = h.source.recvec1.to_python()
                    if test == 0:
                        histogram = None
                        h('vec.fill(0)')
                else:
    
                    histogram = None
                    h('vec.fill(0)')
                histogram = COMM.bcast(histogram, root=k)  # ie root = rank
                #assert isinstance(histogram, list)
                h.pc.broadcast(h.vec,int(k))
                gidn=int(h.vec.x[0])
                h('vec.x[0]=0')
                h.vec.x[0]=0
                gidn = int(gidn)
                testn = h.vec.sum()
    
    
                if int(testn) != 0:
                    print 'cell number ', gidn
                   
                    for i in xrange(0, int(NCELL - 1)):
                       
                        test3 = 0
                        test3=h.pc.gid_exists(i)
                        if test3 != 0:
    
                            h('target=pc.gid2cell(py.int(py.i))')
                            h('print target, " target"')
                            
                            h('print py.i, " py.i"')
                            h('vb[0].copy(target.recvec1,vec,20)')
                            h('storval=normte(vb[0],vec,20)')
                            ent = 0.0
                            ent=float(h.storval.x[2])
                            trentm[i][gidn] = ent
    
        my_trentm = np.zeros_like(trentm)
    
        COMM.Reduce([trentm, MPI.DOUBLE], [my_trentm, MPI.DOUBLE], op=MPI.SUM,
                root=0)
    
        COMM.Barrier()
        return my_trentm           
    
    
    
    import numpy as np
    #global NCELL
    #NCELL= int(np.shape(allmemb)[0])
    
    
    def t_e2(SIZE, NCELL, RANK, trentm,np):
    

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
        h('storval = new Vector()')
        h('objref vb[2]')
        h('vb[0]=new Vector()')
        h('vb[1]=new Vector()')
        h('objref vec')
        h('vec=new Vector()')
    
    
        for k in xrange(0, SIZE):
            for j in xrange(0, int(NCELL-1)):
                if RANK == k:
                    histogram = []
                    test = 0
                    
                    test=int(h.pc.gid_exists(int(j)))
                    if test != 0:
                        h.source=h.pc.gid2cell(int(j))
                        h.vec=h.source.spk_trainb
                        h.vec.x[0]=j
                
                        histogram = h.source.spk_trainb.to_python()
                    if test == 0:
                        histogram = None
                        h('vec.fill(0)')
                else:
    
                    histogram = None
                    h('vec.fill(0)')
                histogram = COMM.bcast(histogram, root=k)  # ie root = rank
                assert isinstance(histogram, list)
                h.pc.broadcast(h.vec,int(k))
                gidn=int(h.vec.x[0])
                h('vec.x[0]=0')
                gidn = int(gidn)
                testn = h.vec.sum()
    
    
                if int(testn) != 0:
                    print 'cell number ', gidn
                   
                    for i in xrange(0, int(NCELL - 1)):
                       
                        test3 = 0
                        test3=h.pc.gid_exists(i)
                        if test3 != 0:
                            h('target=pc.gid2cell(py.int(py.i))')
                    
                            h('storval=normte(target.spk_trainb,vec,20)')
                            ent = 0.0
                            ent=float(h.storval.x[2])
                            trentm[i][gidn] = ent
    
        my_trentm = np.zeros_like(trentm)
    
        COMM.Reduce([trentm, MPI.DOUBLE], [my_trentm, MPI.DOUBLE], op=MPI.SUM,
                root=0)
    
        COMM.Barrier()
        return my_trentm


t=transfer()
t.te()
t.te_2()


#test=0
#histogram=[]
#ent=0.0
#trentm = np.zeros((NCELL, NCELL))
#trentm=t_e(SIZE, NCELL, RANK, trentm,np)
'''
          #(SIZE, NCELL, RANK, histogram, trentm,h)
#trentm=t_e(SIZE, NCELL, RANK, histogram, trentm,h)
#(SIZE, NCELL, RANK, trentm,h)
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

h('objref vec1, vec2 ')
h('vec1=new Vector()')
h('vec2=new Vector()')

COMM.Barrier()
'''
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
            if RANK == s:
                histogram = []
                #print 'j, s', j, ' ', s
                test = 0
                
                h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
            
                if test != 0:
                    h('source=pc.gid2cell(py.int(py.j))')
                    h('vec=source.recvec1')
                    h('vec.x[0]=py.j')
                    histogram = h.source.recvec1.to_python()
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
                    test3 = 0
                    print i, ' i '
                    h('py.test3=py.int(pc.gid_exists(py.int(py.i)))')
                    if test3 != 0:
                        visited[i][gidn] = visited[i][gidn] + 1

                        h('target=pc.gid2cell(py.int(py.i))')
                        h('target.spk_train.printf')

                        h('storval=normte(vec,target.recvec1,20)')
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
#trentm2,visited=t_e2(SIZE, NCELL, RANK, histogram, trentm2, visited, h)
'''
'''
for s in xrange(0, SIZE-1):
    for j in xrange(0, int(NCELL-1)):
        if RANK == s:
            histogram = []

            test = 0
            h('py.test=py.int(pc.gid_exists(py.int(py.j)))')

            if test != 0:
                h('source=pc.gid2cell(py.int(py.j))')

                h('vec=source.recvec1')

                h('vec.x[0]=py.j')

                histogram = h.source.recvec1.to_python()
            else:
                histogram = None
                h('vec.fill(0)')
        else:

            histogram = None
            h('vec.fill(0)')
        histogram = COMM.bcast(histogram, root=s)  # ie root = rank
        h('pc.broadcast(vec,py.int(py.s))')

        testn = 0
        h('py.gidn=int(vec.x[0])')
        h('vec.x[0]=0')
        gidn = int(gidn)
        h('py.testn=int(vec.sum())')
        testn = h.vec.sum()
        if int(testn) != 0:
            for i in xrange(0, int(NCELL - 1)):
                test3 = 0
                h('py.test3=py.int(pc.gid_exists(py.int(py.i)))')
                if test3 != 0:
                    visited[i][gidn] = visited[i][gidn] + 1

                    h('target=pc.gid2cell(py.int(py.i))')
                    h('target.spk_train.printf')

                    h('storval=normte(vec,target.recvec1,20)')
                    ent = 0.0
                    h('py.ent=py.float(storval.x[2])')
                    trentm2[i][gidn] = ent

                    # if(i!=int(gidn)):#if the gids are not the same

                    # vec1 = downsample(histogram, oldrate=40000,newrate=200)
                    # vec2 = downsample(h.target.recvec.to_python(), oldrate=40000,newrate=200)
                    # I think its related to the maximum number of machine operations that permit synchronicity between nodes, being exceeded.
                    # h('vec.resample(vec,1/20)')
                    # h('target.recvec.resample(target.recvec,1/20)')
                    # h('storval=normte(target.recvec.add(-target.recvec.min(),target.recvec),vec.add(-vec.min(),vec),10)')
                    # TO DO UNCOMMENT

                    # ##

COMM.Barrier()
'''
'''
# trentml2=COMM.gather(np.matrix(trentm2, dtype=np.float64),root=0)
# trentm2out=gather(trentm2)
# visitedout=gather(visited)

trentm2,visited=t_e2(SIZE, NCELL, RANK, histogram, trentm2, visited, h)
my_trentm2 = np.zeros_like(trentm2)
COMM.Reduce([trentm2, MPI.DOUBLE], [my_trentm2, MPI.DOUBLE],
            op=MPI.SUM, root=0)

my_visited = np.zeros_like(visited)
COMM.Reduce([visited, MPI.DOUBLE], [my_visited, MPI.DOUBLE],
            op=MPI.SUM, root=0)

# trentm2=COMM.Reduce([np.matrix(trentm2), MPI.DOUBLE], [np.matrix(trentm), MPI.DOUBLE],
#            op=MPI.SUM, root=0)
# msgcv=COMM.Reduce([np.matrix(msgcv), MPI.DOUBLE], [np.matrix(trentm), MPI.DOUBLE],
#            op=MPI.SUM, root=0)
# visited=COMM.Reduce([np.matrix(visited), MPI.DOUBLE], [np.matrix(trentm), MPI.DOUBLE],
#            op=MPI.SUM, root=0)

# visited=COMM.gather(np.matrix(visited, dtype=np.float64),root=0)

COMM.Barrier()

if RANK == 0:

   # summ2 = np.zeros((NCELL,NCELL))
   # summ2 = np.matrix(summ, dtype=np.float64)

   # print np.shape(trentm), ' ', np.shape(summ)
   # visited2 =np.zeros((NCELL,NCELL))
   # for k in xrange(0,SIZE):
   #     np.add(trentml2[k],summ2,summ2)
   #     np.add(visited[k],visited2,visited2)

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

# memory leak or memory intensive?
#
# Same thing, but substitute variables in nTE for SGC and vice versa
#
#
#


'''
