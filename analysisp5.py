#!/usr/bin/python
# -*- coding: utf-8 -*-

from bsmart import granger  # Load the Granger calculation tool
import numpy as np
from matplotlib.colors import LogNorm
import numpy as np
from mpi4py import MPI



#import neuron as h

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

h.xopen('analysis2.hoc')

h('objref storval')
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


ainfplt = gather(np.array(ainf))
'''
Negative Active information means that each cell is a bad predictor of its own activity.

'''  
if RANK == 0:

    fig = plt.figure()
    fig.clf()
    cell_ax = np.linspace(0.0, len(ainfplt), len(ainfplt))#sampling frequency
    plt.plot(cell_ax,ainfplt,'.')
    plt.xlabel('nTE')
    plt.ylabel('cells')
    plt.title('Active Information')
    sfin =  'Active_Informationst'+str(SIZE) + str(NCELL) +'.png'
    fig.savefig(sfin)
    

'''
The reason that this does not work is not because of the h object.
I have seen the h object work in the def prun without being entered as a parameter or a global variable. I think the reason it did not work was more likely to be a problem with updating hoc objects inside python definitions. But I bet I could find an example of this inside pynn.
'''


trentm = np.zeros((NCELL, NCELL))


def broadcastiter(SIZE, NCELL, RANK, trentm, function_x):
  self.SIZE=SIZE
  self.NCELL=NCELL
  self.RANK=RANK
  self.trentm

    storval = h.Vector()
    vec = h.Vector()
    pc = h.ParallelContext()
    rank = int(pc.id())
    nhost = int(pc.nhost())

    for k in xrange(0, self.SIZE):
      for j in xrange(0, int(self.NCELL)):
          test=int(pc.gid_exists(j))
  
                if test != 0:
                    h('source=pc.gid2cell(py.int(py.j))')
                
                    h('vec=source.recvec1')

                    h('vec.x[0]=py.j')

                    histogram = h.source.recvec1.to_python()
                if test == 0:
                    histogram = None
            else:
                histogram = None
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
                    test3 = 0
                    #print i, ' i '
                    h('py.test3=py.int(pc.gid_exists(py.int(py.i)))')
                    #test3=h.pc.gid_exists(i)
                    if test3 != 0:
                        h('target=pc.gid2cell(py.int(py.i))')
                   
                        function_x(target,histogram)


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
    h('storval = new Vector()')
    h('objref vb[2]')
    h('vb[0]=new Vector()')
    h('vb[1]=new Vector()')
    h('objref vec')
    h('vec=new Vector()')


#
# I think the unusual parts of the transfer entropy matrix
# are caused by SIZE and NCELL being incremented one too many times.
# SIZE should only be incremented to be s...SIZE-1
#
    for k in xrange(0, SIZE):
        for j in xrange(0, int(NCELL-1)):
            if RANK == k:
                histogram = []
                #print 'j, s', j, ' ', s
                test = 0
                h.py.j=j
                h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
            
                if test != 0:
                    h('source=pc.gid2cell(py.int(py.j))')
                
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
                    test3 = 0
                    #print i, ' i '
                    h('py.test3=py.int(pc.gid_exists(py.int(py.i)))')
                    #test3=h.pc.gid_exists(i)
                    if test3 != 0:
                        h('target=pc.gid2cell(py.int(py.i))')
                   
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


trentm2 = np.zeros((NCELL, NCELL))
visited = np.zeros((NCELL, NCELL))

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

# memory leak or memory intensive?
#
# Same thing, but substitute variables in nTE for SGC and vice versa
#
#
#


'''
'''
