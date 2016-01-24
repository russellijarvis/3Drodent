
#def spkplt():

for k in xrange(0,SIZE):
    np.add(trentml[k],summ,summ)
    print k
    print 'this value should be increasing: ', np.sum(summ)
    print np.shape(trentml)
    #print summ[0][:]
    #print summ[1][:]
    print 'this value should not be dependent on the previous value: ', np.sum(trentml[k])

if RANK==0:
    print 'idvec[:] ', idvec1[:]
    print 'tvec[:] ', tvec1[:]
    timec2=[]
    tvec=[]
    tvece=[]
    tveci=[]
    idvec=[]
    isivec=[]
    idveci=[]
    isiveci=[]
    idvece=[]
    isivece=[]
    for k in xrange(0,SIZE):
            #timec2=timec2+timec[k]#Concatonate the lists so that they are all the size.

        tvec=tvec+tvec1[k]
        idvec=idvec+idvec1[k]
            #isivec=isivec+isivec1[k]
            
    i=0
    for i in xrange(0,len(idvec)):
        i=int(i)
        for n in xrange(0,len(ie2)-1):
            if int(ie2[n][0])==int(idvec[i]):#if we have selected the idvec from the list

                print ie2[n][0], ie2[n][1]
                #print idvec[i]
                
                    #of cell polarities.
                    #if the polarity is + add the corresponding 
                    #identity and time to the appropriate lists.
                if int(ie2[n][1])==int(0):
                    
                    #print ie2[n][1], 'excitatory'
                    idvece.append(idvec[i])
                    tvece.append(tvec[i])
                if int(ie2[n][1])==int(1):
                    #print ie2[n][1], 'inhibitory'
                    
                    idveci.append(idvec[i])
                    tveci.append(tvec[i])
                        
                    #All the polarities are 0 Excitatory

                    #ie[i][0]=gidn
                    #ie[i][1]=polarity
                    
    h('idvec=py.idvec')
    h('tvec=py.tvec')
    fig=plt.figure()
    fig.clf()
    plt.ylabel('Cell number')
    plt.xlabel('spike time (ms)')
    plt.plot(tveci,idveci,'bo')
    plt.plot(tvece,idvece,'go')
    sfin='spikes'+str(NCELL)+str(SIZE)+'.png'
    fig.savefig(sfin)
    fig=plt.figure()
    fig.clf()
    plt.plot(tvec,idvec,'bo')
        
    sfin='spikesbw'+str(NCELL)+str(SIZE)+'.png'
    fig.savefig(sfin)
    print np.shape(tvece), np.shape(idvece)
    print np.shape(tveci), np.shape(idveci)  
#    h('isivec=py.isivec')
#COMM.Barrier()
