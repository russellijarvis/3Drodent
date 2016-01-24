
def spkplt2(tvec,idvec, ie):
    #global tvec, idvec, ie
#intended only to be executed on RANK 0.
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
    fig = plt.figure()
    fig.clf()
    plt.ylabel('Cell number')
    plt.xlabel('spike time (ms)')
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
    #print np.shape(tvece), np.shape(idvece)
    #print np.shape(tveci), np.shape(idveci)
    #print vecs[:]
    #print vecs2[:]
    return vecs




#GABA neuro transmitters are made by this code.
def mkcells(RANK,NCELL,SIZE,allrows,ie0,ie1,gidvec,cnt):
    '''
    A function for making the cells in the network
    
    '''

    for i in range(RANK, NCELL, SIZE):  # 20 was int(len(allrows))
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
            return h.cells

#print np.shape(allrows2)[0], NCELL, h('cells.count')


