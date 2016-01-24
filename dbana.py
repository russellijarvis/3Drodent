COMM.Barrier()
h('spikeout()')
h('vout()')
COMM.Barrier()
h('xopen("post_analysisp2.hoc")')
h('objref vb[2]')
h('binsz=10//ms')

COMM.Barrier()

n_cell=0
h('py.n_cell=idvec.size')
print 'completed'

print NCELL, ' NCELL', n_cell, ' n_cell'
trentm= [[0 for x in xrange(0,int(NCELL))] for x in xrange(0,int(NCELL))]
h.xopen("/home/zaza3/trunk/examples/expericomp20140421/analysis2.hoc")

h('for i=0, cells.count-1{ cells.o(i).spk_train.printf}')

h('objref vb[2]')
h('vb[0]=new Vector()')
h('vb[1]=new Vector()')
h('objref vec')
h('vec=new Vector()')
for s in xrange(0,SIZE):
    for j in xrange(0,int(NCELL)):
        if (RANK==s):     
            histogram=[]  
            test=0
            test2=0
            h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
            h('py.test2=gidvec.contains(py.j)')
            
            if (test!=0):
                h('source=pc.gid2cell(py.int(py.j))')

                print test==test2
                h('vb[0] = new Vector() //vb stands for vector binned.')
                h('source.spk_train.printf')

                h('vb[0].hist(source.spk_train,0,(tstop+binsz-1)/binsz,binsz)')
                h('vb[0].printf')
                h('vec=source.spk_train')
  
                h('vec.x[0]=py.j')
                h('py.histogram=source.spk_train.to_python()')#.hist(source.spk_train,0,(tstop+binsz-1)/binsz,binsz).to_python()')

                print histogram, 'not indexed'
            else :
                histogram=None
                h('vec.fill(0)')
    
        else: 
            histogram=None
            h('vec.fill(0)')
        data = COMM.bcast(histogram, root=s)#ie root = rank
        h('pc.broadcast(vec,py.int(py.s))')

        testn=0
        h('py.testn=vec.sum()')
        if (testn!=0):
        #if (data!=None):#and(RANK!=s):
            for i in xrange(0,int(NCELL)):
                h('print "data recieved"')
                h('vec.printf')
                h('py.gidn=vec.x[0]')
                h('vec.x[0]=0')
                
                if(i!=int(gidn)):#if the gids are not the same
                    test3=0
                    h('py.test4=gidvec.contains(py.i)')
                    if (test3!=0):
                      h('target=pc.gid2cell(py.int(py.i))')
                      #print test4==test3
                      #print 'target.spk_train.printf, was it empty, or just one element?'
                      
                      h('target.spk_train.printf')
                      h('target.spk_train.x[0]=0')
                      #print 'target.spk_train.printf, was it empty, or just one element?'
                      
                      h('vb[1] = new Vector() //vb stands for vector binned.')
                      h('vb[1].hist(target.spk_train,0,(tstop+binsz-1)/binsz,binsz)')
                      h('vb[1].x[0]=0')   
                      h('vb[0].hist(vec,0,(tstop+binsz-1)/binsz,binsz)')
                      
                      h('vb[0].x[0]=0')#Restore the Vector to having no spikes at time zero. In fact the gid vector was being stored in this element for convenience.
                      #This may confuse subsequent developers.
                      #h('storval=normte(vec,target.spk_train,10)')
                      #h('storval=normte(vb[0],vb[1],5)')
                      h('vec.printf')
                      h('target.spk_train.printf')
                      h('print py.int(py.i)," ",py.int(py.gidn)')
                      h('storval=GetTENQ(vec,target.spk_train,5,py.int(py.i),py.int(py.gidn))')
                     
                      h('storval=GetTENQ(vb[0],vb[1],5,py.int(py.i),py.int(py.gidn))')
                      h('print storval')
                      
                      #h('// from(0) to(1) TE(2) NTE(3) HX2|X2P(4) prefdir(5) TEshufavg(6) TEshufstd(7) sig(8)')	  
                      ent=0.0
                      h('py.ent=py.float(storval.x[2])')
                      trentm[i][gidn]=trentm[i][gidn]+dirin
                      print trentm[i][gidn], ' trenm element ', COMM.rank,' ',n_cell
                      for k in xrange(0,NCELL):
                          h('storval1=normte(vb[0],vb[1],py.int(py.k)*2)')
                          h('storval=GetTENQ(vb[0],vb[1],py.int(py.k)*2,py.int(py.i),py.int(py.gidn))')
                     
                          h('py.ent2=py.float(storval1.x[2])')
                          surface[i][gidn][k]=surface[i][gidn][k]+ent2
COMM.Barrier()
print NCELL, ' NCELL ', n_cell, ' n_cell '   

trentm=COMM.gather(np.matrix(trentm, dtype=np.float32),root=0)


if RANK==0:
   print np.shape(surface)

   summ = [[0 for x in xrange(0,int(NCELL))] for x in xrange(0,int(NCELL))]
   sumsurface=np.empty((NCELL,NCELL,NCELL))


   summ = np.matrix(summ, dtype=np.float32)
   print np.shape(trentm), ' ', np.shape(summ)

   for k in xrange(0,SIZE):
        np.add(trentm[k],summ,summ)
        #print np.shape(surface[k])

        #np.add(surface[k],sumsurface,sumsurface)
   fin =0
    #fin = fin + 1
   fig = plt.figure()
   fig.clf()


   im = plt.imshow(summ, interpolation='nearest')

   plt.autoscale(True)
   plt.colorbar(im)
   plt.xlabel('columns = targets')
   plt.ylabel('rows = sources')
   plt.title('Transfer Entropy Matrix')
   plt.grid(True)
        
   sfin = 'Transfer_Entropy_matrix.png'
   fig.savefig(sfin)
   h('system("eog Transfer_Entropy_matrix.png")')
        #timec2=timec2+timec[k]#Concatonate the lists so that they are all the size.
surface=COMM.gather(surface,root=0)

h('vout()')

if RANK==0:
   print np.shape(surface)

#   summ = [[0 for x in xrange(0,int(NCELL))] for x in xrange(0,int(NCELL))]
   sumsurface=np.empty((NCELL,NCELL,NCELL))


   summ = np.matrix(summ, dtype=np.float32)
   print np.shape(trentm), ' ', np.shape(summ)

   for k in xrange(0,SIZE):

        np.add(surface[k],sumsurface,sumsurface)

   from mpl_toolkits.mplot3d import Axes3D
   from matplotlib import cm
   from matplotlib.ticker import LinearLocator, FormatStrFormatter
   fig = plt.figure()
   ax = fig.gca(projection='3d')
   indexs=int(np.linspace(1,NCELL,NCELL))#cast to int.
   #surf = ax.plot_surface(sumsurface[indexs][:][:],sumsurface[:][indexs][:],sumsurface[:][:][indexs] rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
   ax.set_zlim(-1.01, 1.01)
   ax.zaxis.set_major_locator(LinearLocator(10))
   ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
   fig.colorbar(surf, shrink=0.5, aspect=5)
   fig.savefig(sfin)
   sfin = 'Transfer_Entropy_surface.png'

   h('system("eog Transfer_Entropy_surface.png")')

