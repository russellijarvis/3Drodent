
ent_isicp=[]
actinf=[]
actinf2=[]
h('py.ent_isicp=ent_isicp.to_python()')
h('py.actinf=actinf.to_python()')
h('py.actinf2=actinf2.to_python()')
input_isi_e=np.empty(numcell,float)
input_isi_e[:]=veisi


  
fig = p1.figure()
fig.clf()
p1.hold(True)
p1.scatter(across,ent_isicp)
p1.plot(across,input_isi_e)
p1.hold(False)
p1.xlabel('cell number')
p1.ylabel('bits')
sfin = 'ISI_based_ent_' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin)
 
""" 
fig = p1.figure()
fig.clf()
p1.scatter(across,actinf)
p1.title('ISI active information')
sfin = 'ISI_based_activ_inf_' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin)
"""
fig = p1.figure()
fig.clf()
p1.scatter(across,(actinf2))
p1.title('time binned active information')
sfin = 'Active Information' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin)



#print 'is the array symetric?'
#(np.array(entfa).transpose(1, 0) == arr).all()


#ax = p1.subplot(111)
fig = p1.figure()
fig.clf()
p1.hold(True)
#p1.bar(np.arange(numcell)-0.2,actinf2,width=0.2,color='b',align='center')
p1.bar(np.arange(int(numcell/4))-0.2, (actinf2[:int(numcell/4)]),label='active inf',width=0.2,color='g',align='center')
p1.bar(np.arange(int(numcell/4)), impinging[:int(numcell/4)],label='recieving inf',width=0.2,color='b',align='center')
p1.bar(np.arange(int(numcell/4))+0.2, predictive[:int(numcell/4)],label='transmitting inf',width=0.2,color='r',align='center')
p1.hold(False)
p1.xlabel('cell number')
p1.ylabel('active information/impinging information/predictive information')
p1.title('Information transfer versus cell')
p1.legend(bbox_to_anchor=(0.0, 1.02, 1.0, .102), loc=9, ncol=2, \
          mode='expand', borderaxespad=0.0)
sfin = 'Active_impinging_predictive1' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin)

fig = p1.figure()
fig.clf()
p1.hold(True)
#p1.bar(np.arange(numcell)-0.2,actinf2,width=0.2,color='b',align='center')
p1.bar(across[int(numcell/4):int(2*numcell/4)]-0.2, (actinf2[int(numcell/4):int(2*numcell/4)]),label='active inf',width=0.2,color='g',align='center')
p1.bar(across[int(numcell/4):int(2*numcell/4)], impinging[int(numcell/4):int(2*numcell/4)],label='recieving inf',width=0.2,color='b',align='center')
p1.bar(across[int(numcell/4):int(2*numcell/4)]+0.2, predictive[int(numcell/4):int(2*numcell/4)],label='transmitting inf',width=0.2,color='r',align='center')
p1.hold(False)
p1.xlabel('cell number')
p1.ylabel('active information/impinging information/predictive information')
p1.title('Information transfer versus cell')
p1.legend(bbox_to_anchor=(0.0, 1.02, 1.0, .102), loc=9, ncol=2, \
          mode='expand', borderaxespad=0.0)
sfin = 'Active_impinging_predictive2' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin)


fig = p1.figure()
fig.clf()
p1.hold(True)
#p1.bar(np.arange(numcell)-0.2,actinf2,width=0.2,color='b',align='center')
p1.bar(across[int(2*numcell/4):int(3*numcell/4)]-0.2, (actinf2[int(2*numcell/4):int(3*numcell/4)]),label='active inf',width=0.2,color='g',align='center')
p1.bar(across[int(2*numcell/4):int(3*numcell/4)], impinging[int(2*numcell/4):int(3*numcell/4)],label='recieving inf',width=0.2,color='b',align='center')
p1.bar(across[int(2*numcell/4):int(3*numcell/4)]+0.2, predictive[int(2*numcell/4):int(3*numcell/4)],label='transmitting inf',width=0.2,color='r',align='center')
p1.hold(False)
p1.xlabel('cell number')
p1.ylabel('active information/impinging information/predictive information')
p1.title('Information transfer versus cell')
p1.legend(bbox_to_anchor=(0.0, 1.02, 1.0, .102), loc=9, ncol=2, \
          mode='expand', borderaxespad=0.0)
sfin = 'Active_impinging_predictive3' + str(int(h.prunenet)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin)



fig = p1.figure()
fig.clf()
p1.hold(True)
#p1.bar(np.arange(numcell)-0.2,actinf2,width=0.2,color='b',align='center')
p1.bar(across[int(floor(3*numcell/4)):int(numcell)]-0.2, (actinf2[int(floor(3*numcell/4)):int(numcell)]),label='active inf',width=0.2,color='g',align='center')
p1.bar(across[int(floor(3*numcell/4)):int(numcell)], impinging[int(floor(3*numcell/4)):int(numcell)],label='recieving inf',width=0.2,color='b',align='center')
p1.bar(across[int(floor(3*numcell/4)):int(numcell)]+0.2, predictive[int(floor(3*numcell/4)):int(numcell)],label='transmitting inf',width=0.2,color='r',align='center')
p1.legend(bbox_to_anchor=(0.0, 1.02, 1.0, .102), loc=9, ncol=2, \
          mode='expand', borderaxespad=0.0)
p1.hold(False)
p1.xlabel('cell number')
p1.ylabel('active information/impinging information/predictive information')
p1.title('Information transfer versus cell')

sfin = 'Active_impinging_predictive4' + str(int(h.prunenet))+str(int(h.iw)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin)



fig = p1.figure()
fig.clf()
p1.scatter(across,predictive[:len(across)])
sfin = 'Predictive_information' + str(int(h.prunenet))+str(int(h.iw)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin)


fig = p1.figure()
fig.clf()
p1.scatter(across,impinging[:len(across)])
sfin = 'Impinging_information2' + str(int(h.prunenet))+str(int(h.iw)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin)





#h.ent_table_sq_isi()

fin = fin + 1 
fig = p1.figure()
fig.clf()
im = p1.imshow(dirin2, interpolation='nearest')
p1.colorbar(im)


p1.xlabel('columns = targets')
p1.ylabel('rows = sources')
p1.title('ISI based NTE')
p1.autoscale(True)
p1.grid(True)
sfin = 'ISI_based_NTE_' + str(int(h.prunenet))+str(int(h.iw)) \
    + str(tstop) + str(h.plastic) + str(int(h.ff)) + str(ncell) \
    + str(int(fin)) + str(int(h.minr)) + '.png'
fig.savefig(sfin)

