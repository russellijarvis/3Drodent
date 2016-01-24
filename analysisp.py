#!/usr/bin/env python

from bsmart import granger  # Load the Granger calculation tool
from matplotlib.colors import LogNorm

def gather(data):
    #https://github.com/NeuralEnsemble/PyNN/blob/master/examples/distrib_example.py
    #Andrew Davison

    assert isinstance(data, np.ndarray)
    # first we pass the data size
    size = data.size
    sizes = COMM.gather(size, root=0) or []
    # now we pass the data
    displacements = [sum(sizes[:i]) for i in range(len(sizes))]
    print(COMM.rank, "sizes=", sizes, "displacements=", displacements)
    gdata = np.empty(sum(sizes))
    COMM.Gatherv([data, size, MPI.DOUBLE], [gdata, (sizes,displacements), MPI.DOUBLE], root=0)
    return gdata

trentm =np.zeros((NCELL,NCELL))
trentm2 =np.zeros((NCELL,NCELL))
visited=np.zeros((NCELL,NCELL))
msgcv = np.zeros((NCELL,NCELL))


h.xopen("analysis2.hoc")
#h('xopen("isi.hoc")')

surface=np.empty((NCELL,NCELL,NCELL))
surfaceg=np.empty((NCELL,NCELL,NCELL))

from pyhoc import downsample
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

for s in xrange(0,SIZE):
    for j in xrange(0,int(NCELL)):
        if (RANK==s):     
            histogram=[]  

            test=0
            h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
            
            if (test!=0):
                h('source=pc.gid2cell(py.int(py.j))')
                
                h('vec=source.spk_trainb')
                
                h('vec.x[0]=py.j')
                
                #histogram=h.source.recvec.to_python()
                histogram=h.source.spk_trainb.to_python()
            if (test==0):
                histogram=None
                h('vec.fill(0)')
    
        else: 
            histogram=None
            h('vec.fill(0)')
        histogram = COMM.bcast(histogram, root=s)#ie root = rank
        h('pc.broadcast(vec,py.int(py.s))')

#       testn=0
        h('py.gidn=int(vec.x[0])')
        h('vec.x[0]=0')
        gidn=int(gidn)
        h('py.testn=int(vec.sum())')
        testn=h.vec.sum()
        print 'testn = ', testn, ' histogram= ', histogram
        #print histogram==testn
        #if histogram!=None:
        if int(testn)!=0:    
            print 'cell number ', gidn
            h('vec.printf')
            for i in xrange(0,int(NCELL-1)):
                test3=0
                print i, ' i '
                h('py.test3=py.int(pc.gid_exists(py.int(py.i)))')
                if (test3!=0):
                    h('target=pc.gid2cell(py.int(py.i))')
        
                      ## Here is the problem. Histograming creates a vector with only 1, 0 element. The problem is in the Vector data structure, and its inability to get updated in this context.
                      ## I will need to work around it somehow. Possibly by using different data structures
                      #Its only the binning that does not work here, so why not do it somewhere else, when then the spike bins are first attributed to the cells.

                      

                      
                      
                    h('storval=normte(target.spk_trainb,vec,20)')
                    ent=0.0
                    h('py.ent=py.float(storval.x[2])')
 
                    trentm[i][gidn]=ent

 


COMM.Barrier()

trentmout=gather(trentm)

#trentml=COMM.gather(np.matrix(trentm, dtype=np.float64),root=0)
COMM.Barrier()
############


print 'np.sum(trentm)= ', np.sum(trentm)

if RANK==0:
   print np.shape(surface)

   sumsurface=np.empty((NCELL,NCELL,NCELL))
   
   
   summ = np.zeros((NCELL,NCELL))
   summ2 = np.zeros((NCELL,NCELL))

   ident = np.matrix(np.identity(NCELL))
   actinf = summ*ident

   print np.shape(trentm), ' ', np.shape(summ)

   for k in xrange(0,SIZE):
        np.add(trentml[k],summ,summ)
        
   fig = plt.figure()
   fig.clf()
   im = plt.imshow(summ, interpolation='nearest',norm=LogNorm())
   plt.autoscale(True)
   plt.colorbar(im)
   plt.xlabel('columns = targets')
   plt.ylabel('rows = sources')
   plt.title('Transfer Entropy Matrix')
   plt.grid(True)
   sfin = str(SIZE)+str(NCELL)+'Transfer_Entropy_matrix.png'
   fig.savefig(sfin)


############

h('objref vec1, vec2 ')
h('vec1=new Vector()')
h('vec2=new Vector()')

COMM.Barrier()

for s in xrange(0,SIZE):
    for j in xrange(0,int(NCELL)):
        if (RANK==s):     
            histogram=[]  
        
            test=0
            h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
            
            if (test!=0):
                h('source=pc.gid2cell(py.int(py.j))')
                
                h('vec=source.recvec1')
  
                h('vec.x[0]=py.j')
                
                histogram=h.source.recvec.to_python()
            else :
                histogram=None
                h('vec.fill(0)')
    
        else: 
            histogram=None
            h('vec.fill(0)')
        histogram = COMM.bcast(histogram, root=s)#ie root = rank
        h('pc.broadcast(vec,py.int(py.s))')

        testn=0
        h('py.gidn=int(vec.x[0])')
        h('vec.x[0]=0')
        gidn=int(gidn)
        h('py.testn=int(vec.sum())')
        testn=h.vec.sum()
        if int(testn)!=0:
            for i in xrange(0,int(NCELL-1)):
                test3=0
                h('py.test3=py.int(pc.gid_exists(py.int(py.i)))')
                if (test3!=0):
                    visited[i][gidn]=visited[i][gidn]+1
 
                    h('target=pc.gid2cell(py.int(py.i))')
                    h('target.spk_train.printf')

                    vec1 = downsample(histogram, oldrate=40000,newrate=200)
                    vec2 = downsample(h.target.recvec.to_python(), oldrate=40000,newrate=200)
                    h('storval=normte(vec,target.recvec1,20)')
                    ent=0.0
                    h('py.ent=py.float(storval.x[2])')
                    trentm2[i][gidn]=ent
                    if(i!=int(gidn)):#if the gids are not the same

 
                       
                      # I think its related to the maximum number of machine operations that permit synchronicity between nodes, being exceeded.
                      #h('vec.resample(vec,1/20)')
                      #h('target.recvec.resample(target.recvec,1/20)')
                      #h('storval=normte(target.recvec.add(-target.recvec.min(),target.recvec),vec.add(-vec.min(),vec),10)')
                      #TO DO UNCOMMENT
                      
                      ###
                          if (np.std(vec1)!=0) and (np.std(vec2)!=0) :
                              order = 10
                              rate = 200
                              maxfreq = 0


                              npts = len(vec1)#len(h.vec1.to_python())
                              fs = 200
                              n = len(vec2)#len(h.vec1.to_python())
                              freq = 100
                              p = 15
                              ntrls = 1
                              (  # lag = 2 bins of 5ms each.
                                  F,
                                  pp,
                                  cohe,
                                  Fx2y,
                                  Fy2x,
                                  Fxy,
                              ) = granger(vec1, vec2, 1)

                              msgcv[i][gidn] = np.mean(Fx2y)
                      




COMM.Barrier()
############



#trentml2=COMM.gather(np.matrix(trentm2, dtype=np.float64),root=0)
trentm2out=gather(trentm2)
visitedout=gather(visited)
#visited=COMM.gather(np.matrix(visited, dtype=np.float64),root=0)

COMM.Barrier()

if RANK==0:

   summ2 = np.zeros((NCELL,NCELL))
   summ2 = np.matrix(summ, dtype=np.float64)

   #print np.shape(trentm), ' ', np.shape(summ)
   visited2 =np.zeros((NCELL,NCELL))
   for k in xrange(0,SIZE):
        np.add(trentml2[k],summ2,summ2)
        np.add(visited[k],visited2,visited2)


   fig = plt.figure()
   fig.clf()


   im = plt.imshow(visited2, interpolation='nearest')

   plt.autoscale(True)
   plt.colorbar(im)
   plt.xlabel('columns = targets')
   plt.ylabel('rows = sources')
   plt.title('visited')
   plt.grid(True)
        
   sfin = str(SIZE)+str(NCELL)+'Visited.png'
   fig.savefig(sfin)
   fig = plt.figure()
   fig.clf()


   im = plt.imshow(summ2, interpolation='nearest',norm=LogNorm())

   plt.autoscale(True)
   plt.colorbar(im)
   plt.xlabel('columns = targets')
   plt.ylabel('rows = sources')
   plt.title('Transfer Entropy Matrix Two')
   plt.grid(True)
        
   sfin = str(SIZE)+' '+str(NCELL)+'Transfer_Entropy_matrix2.png'
   fig.savefig(sfin)


COMM.Barrier()

msgcv=gather(msgcv)

#msgcv2=COMM.gather(np.matrix(msgcv, dtype=np.float64),root=0)
if RANK==0:
   summg = np.zeros((NCELL,NCELL))  
   for k in xrange(0,SIZE):
       summg=summg+msgcv2[k]
  
   fig = plt.figure()
   fig.clf()
   im = plt.imshow(np.matrix(summg), interpolation='nearest',norm=LogNorm())
   plt.autoscale(True)
   plt.colorbar(im)
   plt.xlabel('columns = targets')
   plt.ylabel('rows = sources')
   plt.title('SGC Matrix')
   plt.grid(True)
        
   sfin='sgctest.png'
   fig.savefig(sfin)
#   h('system("eog sgctest.png")')
   #h('system("eog SGC_matrix.png")')
        #timec2=timec2+timec[k]#Concatonate the lists so that they are all the size.
surfaceg=COMM.gather(surfaceg,root=0)

'''
work1=[]
work2=[]
COMM.Barrier()

for j in xrange(0,int(NCELL)):
    h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
    if (test!=0):
        h('source=pc.gid2cell(py.int(py.j))')
        #h('vec=source.spk_trainb')
        #h('vec.x[0]=py.j')                
        histogram=h.source.recvec.to_python()
        work1.append(histogram)

COMM.Barrier()

#work2=COMM.Alltoall(work1)
'''
#
# Same thing, but substitute variables in nTE for SGC and vice versa
#
#
#

'''

trentm2= [[0 for x in xrange(0,int(NCELL))] for x in xrange(0,int(NCELL))]
#h.xopen("/home/zaza3/trunk/examples/expericomp20140421/analysis2.hoc")
#surface=np.empty((NCELL,NCELL,NCELL))
#surfaceg=np.empty((NCELL,NCELL,NCELL))

msgcv2 = [[0 for x in xrange(0, int(NCELL))] for x in xrange(0, int(NCELL))]

for s in xrange(0,SIZE):
    for j in xrange(0,int(NCELL)):
        if (RANK==s):     
            histogram=[]  
            #histogram2=[]  

            test=0
            #test2=0
            h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
            #h('py.test2=gidvec.contains(py.j)')
            
            if (test!=0):
                h('source=pc.gid2cell(py.int(py.j))')
                
                #print test==test2
                #h('vb[0] = new Vector() //vb stands for vector binned.')
                #h('source.spk_train.printf')

                #h('vb[0].hist(source.spk_train,0,(tstop1+binsz-1)/binsz,binsz)')
                #h('vb[0].printf')
                h('vec=source.spk_train')
  
                h('vec.x[0]=py.j')
                #h('py.histogram=source.spk_train.to_python()')#.hist(source.spk_train,0,(tstop1+binsz-1)/binsz,binsz).to_python()')
                #histogram2=h.source.spk_train.to_python()
                #
                #
                #THIS LINE WORKS
                #histogram=h.recvectors[int(j)].to_python()
                
                histogram=h.source.recvec.to_python()
                #print histogram, 'not indexed'
            else :
                histogram=None
                h('vec.fill(0)')
    
        else: 
            histogram=None
            h('vec.fill(0)')
        data = COMM.bcast(histogram, root=s)#ie root = rank
        h('pc.broadcast(vec,py.int(py.s))')

        testn=0
        h('py.gidn=int(vec.x[0])')
        h('vec.x[0]=0')
        gidn=int(gidn)
        h('py.testn=int(vec.sum())')
        testn=h.vec.sum()
        if histogram!=None:
        #if (testn!=0):
            print 'cell number ', gidn
            h('vec.printf')
            #break
            #i=0
            for i in xrange(0,int(NCELL-1)):
                
                if(i!=int(gidn)):#if the gids are not the same
                    test3=0
                    print i, ' i '
                    h('py.test3=py.int(pc.gid_exists(py.int(py.i)))')
                    if (test3!=0):
                      h('target=pc.gid2cell(py.int(py.i))')
                      h('target.spk_train.printf')

                      vec1 = downsample(histogram, oldrate=40000,newrate=200)
                      vec2 = downsample(h.target.recvec.to_python(), oldrate=40000,newrate=200)
                      vec1[:]=vec1[:]+np.min(vec1)
                      vec2[:]=vec2[:]+np.min(vec2)
  
                      order = 10

                      rate = 200
                      maxfreq = 0
                      #vec1=h.vec1.to_python()
                      #vec2=h.vec2.to_python()


   

                      #trentm[i][gidn]=trentm[i][gidn]+ent
      

                      #h('target.spk_train.x[0]=0')
                      h('vb[1] = new Vector() //vb stands for vector binned.')
                      h('vb[1].hist(target.spk_train,0,(tstop1+binsz-1)/binsz,binsz)')

                      h('vb[1].x[0]=0')   
                      h('vb[0].hist(vec,0,(tstop1+binsz-1)/binsz,binsz)')
                      
                      h('vb[0].x[0]=0')#Restore the Vector to having no spikes at time zero. In fact the gid vector was being stored in this element for convenience.
                
                      h('storval=normte(py.np.abs(py.vec1),py.np.abs(vec2),20)')
                      #The only way these vectors can have the same length is through histograming?
                      #h('storval=normte(vb[0],vb[1],20)')
                      h('print vb[0]')
                      h('vb[0].printf')
                      h('print vb[1]')
                      h('vb[1].printf')

                      npts = len(h.vb[0].to_python())#len(h.vec1.to_python())
                      fs = 200
                      n = len(h.vb[0].to_python())#len(h.vec1.to_python())
                      freq = 100
                      p = 15
                      ntrls = 1


                      (  # lag = 2 bins of 5ms each.
                          F,
                          pp,
                          cohe,
                          Fx2y,
                          Fy2x,
                          Fxy,
                          ) = granger(h.vb[0].to_python(), h.vb[1].to_python(), 1)

                      msgcv2[i][gidn] = np.mean(Fx2y)
                      #for k in xrange(1,5):
                      #    (  # lag = 2 bins of 5ms each.
                      #        F,
                      #        pp,
                      #        cohe,
                      #        Fx2y,
                      #        Fy2x,
                      #        Fxy,
                      #        ) = granger(h.vb[0].to_python(), h.vb[1].to_python(), int(k))
                          #surfaceg[i][gidn][k]=surfaceg[i][gidn][k]+np.mean(Fx2y)



                     
                      #h('storval=GetTENQ(vb[0],vb[1],5,py.int(py.i),py.int(py.gidn))')
                      
                      #h('// from(0) to(1) TE(2) NTE(3) HX2|X2P(4) prefdir(5) TEshufavg(6) TEshufstd(7) sig(8)')	  
                      ent=0.0
                      h('py.ent=py.float(storval.x[2])')
                      #h('if (storval.x[2]!=0){ print storval.x[2]}')
                      trentm2[i][gidn]=ent
                      #trentm[i][gidn]=trentm[i][gidn]+ent
                      #print trentm[i][gidn], ' trenm element ', COMM.rank,' ',n_cell
                      #for k in xrange(0,NCELL):
                      #    h('storval1=normte(vb[0],vb[1],py.int(py.k)*2)')
                      #    h('storval=GetTENQ(vb[0],vb[1],py.int(py.k)*2,py.int(py.i),py.int(py.gidn))')
                     
                       #   h('py.ent2=py.float(storval1.x[2])')
                       #   surface[i][gidn][k]=surface[i][gidn][k]+ent2

COMM.Barrier()

'''
