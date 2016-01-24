



fig01 = drawmatrix_channels(np.matrix(msgcv), np.array(roi_names), size= [10., 10.], color_anchor=0)
sfin='granger_matrixfinal'+str(h.prunenet)+str(fin)+'.png' 
fig01.savefig(sfin)
fin=fin+1




#low_values_indices = np.mean(msgcv)+np.std(msgcv) < 0#///  # Where values are low
#msgcv[low_values_indices] = 0#//  # All low values set to 0
 
#fig02 = drawgraph_channels(np.matrix(msgcv), np.array(roi_names), size= [10., 10.], color_anchor=0)
#fig02.savefig('chanel_type_graph.png')

#fig03 = draw_graph(nx.to_networkx_graph(np.matrix(msgcv),create_using=nx.DiGraph()))
#sfin='grange_graph'+str(h.prunenet)+str(fin)+'.png' 
#fig02.savefig(sfin)
#sfin='grange_graphfinal'+str(h.prunenet)+str(fin)+'.png' 

#fig03.savefig(sfin)
fig=p1.figure(fin)
fig.clf()
p1.imshow(msgcv, interpolation='nearest')
sfin='granger_matrix_imshow'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000)

  
#msgcv=np.array(msgcv) 
#low_values_indices = np.mean(msgcv)+np.std(msgcv) < 0#///  # Where values are low
#msgcv[low_values_indices] = 0#//  # All low values set to 0





dirfinals =nx.DiGraph() 
dirfinals2 =nx.DiGraph() 
#dirfinals.add_nodes_from(msgcv) 
for i in xrange(0,len(msgcv)):
 for j in xrange(0,len(msgcv)):
  if msgcv[i][j]>(np.mean(msgcv)+2.5*np.std(msgcv)):
    dirfinals.add_edge(i,j,weight=msgcv[i][j]) 
  if msgcv2[i][j]>(np.mean(msgcv)+2.5*np.std(msgcv)):
    dirfinals2.add_edge(i,j,weight=msgcv2[i][j]) 

fin+=1
fig = p1.figure(fin)
nx.draw(dirfinals)#, roi_names)
sfin='granger_graph_final_withoutsubtraction'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000)
#fig05 = draw_graph(dirfinals)
#fig05.savefig('experiment.png')


fin+=1
fig = p1.figure(fin)
nx.draw(dirfinals)#, roi_names)
sfin='granger_graph_final_withsubtraction'+str(h.prunenet)+str(fin)+'.png' 

fig=p1.figure(fin)
fig.clf()
#interactive(True)
#p1.plot((xf,yf1))
#p1.title("Raster Plot")
p1.hold(True)
colors=array([[0.42,0.67,0.84],[0.50,0.80,1.00],[0.90,0.32,0.00],[0.34,0.67,0.67],[0.42,0.82,0.83],[0.90,0.59,0.00],[0.33,0.67,0.47],[0.42,0.83,0.59],[0.90,0.76,0.00],[1.00,0.85,0.00],[0.71,0.82,0.41],[0.57,0.67,0.33],[1.00,0.38,0.60]]) # Colors for each cell population
j=len(colors)-1#12
p1.plot(tvec,intervec,'bo',label='inhibitory interneuron')#,c=colors[j], markeredgecolor = 'none')
j-=1
p1.plot(tvec,pyra,'g^')#o',c=colors[j], markeredgecolor = 'none')

p1.plot(tvec,zerovec,'g^', label='pyramidal cell')#,c=colors[j], markeredgecolor = 'none')
j-=1
p1.plot(record_inputs,vecin,'ro',linewidth=15, label='synapse input stim')#, markeredgecolor = 'none')
p1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=9,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(1000)
#manrnpxtime=h.max(h.tvec)
p1.xlim(0,maxtime) # To convert to seconds
p1.ylim(-2,int(numcell)) # Just larger than the number of cells in the model
p1.ylabel("Cell number")
p1.xlabel("spike time (ms)")
#Save fig hangs the program for some reason.
#p1.savefig("ff.png") # generates 'libpng error: zlib error' under nrniv
#p1.savefig("ff.eps")
#p1.axis.set_xticks(r_[0:maxtime+1]) # Don't show half seconds

p1.hold(False)
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000)   
fin+=1
spc=0 # Only one column, so pick it

	





sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(fin)
fig.clf()
fin+=1
##
p1.hold(True)
tc=np.array(time_courses[int(1)])
N=len(tc)
t = np.linspace(0.0, 0.025*N, N)
t=np.array(t)

tc=np.array(time_courses[int(10)])
str3='cell number= '+str(10)
p1.plot(t[0:7555],tc[0:7555],linewidth=1.5,label=str3)
#p1.plot(tc[0:25],out1[0:25])
p1.title('pyramidal neuron membrane potential')
p1.xlabel("ms")
p1.ylabel("mV")
p1.legend()

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(fin)
fig.clf()
fin+=1


for i in xrange(0,int(len(time_courses))):
 string='voltages'
 bc=' cells i= '+str(int(i)) +'ff= '+str(int(h.ff))+'prune net= '+str(int(h.prunenet))
 string=string+str(i)+bc 
 tc=np.array(time_courses[int(i)])
 p1.plot(t[500:1500],tc[500:1500],linewidth=3)
 #p1.plot(tc[0:25],out1[0:25])
 p1.title('cells index membrane potential')
 p1.xlabel("ms")
 p1.ylabel("mV")

 #out1=downsample(time_courses[int(i)],oldrate=40000,newrate=200)
 #downsampled.append(out1)





numcell=int(h.ncell)
"""
p1.hold(False)
vsum1 = np.array(vsum1)
vsum2 = np.array(vsum2)
#vtotal=vtotal/2 #so its an average again.
in_trace = np.array(in_trace)
out_trace = np.array(in_trace)
in_trace=in_trace/int(numcell)

"""

###
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
#belongs to above
fig=p1.figure(fin)
fig.clf()
fin+=1
###
p1.title("Entropy Plot")
p1.hold(True)
j=len(colors)-1

across=np.arange(0,int(numcell),1)
#p1.plot(tvec,intervec,'bs')#p1.plot(tvec,pyra,'g^')
p1.plot(across,input_marke,'r-',linewidth=2,label='H(x) of synapse input stim')#,c=colors[j], markeredgecolor = 'none')
j-=1
p1.plot(across,pspke,'go', label='H(x) of cell num spike train')
p1.xlabel('cell number')
p1.ylabel('bits/sec')
p1.hold(False)
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)



j=len(colors)-1

across=np.arange(0,int(numcell),1)
#p1.plot(tvec,intervec,'bs')#p1.plot(tvec,pyra,'g^')
divin=np.zeros(int(numcell))
#need to use numpy to remove inf
np.divide(input_marke,ratein,divin)
winfs=isinf(divin) # replace zeros with infinity
divin[winfs]=40
j-=1
divout=np.zeros(int(numcell))
np.divide(pspke,rates,divout)

p1.hold(True)

##
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(fin)
fig.clf()
fin+=1
##
p1.plot(across,divin,'r-',linewidth=2,label='H(x) of synapse input stim')#,'b|',c=colors[j], markeredgecolor = 'none')

p1.plot(input_marke,ratein,'go',label='H(x) of cell num spike train')

#p1.plot(across,divout,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)
p1.title("Entropy divided by rate")
p1.xlabel('cell number')
p1.ylabel('bits/sec')

###
p1.hold(False)
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(fin)
fin+=1
p1.hold(True)
###

p1.title("Lempel Ziv")
j=len(colors)-1

p1.plot(across,input_markl,'r-',linewidth=2,label='H(x) of synapse input stim')#,'|',c=colors[j],linewidth=5)
j-=1
p1.plot(across,pspkl,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(h.tstop)

p1.xlabel("Neuron number")
p1.ylabel("Bits/sec")
p1.hold(False)
p1.title("Lempel Ziv Divided by rate")
p1.hold(True)
j=len(colors)-1

divin=np.zeros(int(numcell))
#need to use numpy to remove inf
np.divide(input_markl,ratein,divin) #dividing a large num by small number creates
#approaching infinitely large number.
winfs=isinf(divin) # replace zeros with infinity
divin[winfs]=40
j-=1
divout=np.zeros(int(numcell))
np.divide(pspkl,rates,divout)


p1.plot(across,divin,'r-',linewidth=2,label='H(x) of synapse input stim')#,'|',c=colors[j],linewidth=5)
j-=1
p1.plot(divout,rates,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(h.tstop)

p1.xlabel("Neuron number")
p1.ylabel("Bits/sec")

###
p1.hold(False)
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fin+=1
fig=p1.figure(fin)
#p1.hold(True)
# Two subplots, the axes array is 1-d
f, axarr = p1.subplots(2, sharex=True)
axarr[0].plot(across,crlts,'go', markeredgecolor = 'none')
axarr[0].set_title('correlations between input and cell number')
axarr[1].plot(across,fhv,'go', markeredgecolor = 'none')
axarr[1].set_title('nTE between input and cell number')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)
p1.xlabel('cell number')


p1.hold(False)
maxtime=int(h.tstop)
#p1.xlabel("Neuron number")
#p1.ylabel("nTE")
#p1.show()
#axarr[1].xlabel('neuron number')
""" 
"""
p1.title("nTE and Correlations")
###
j=len(colors)-1
p1.plot(across,crlts,'o',c=colors[j], markeredgecolor = 'none')
j-=1
p1.plot(across,fhv,'o',c=colors[j], markeredgecolor = 'none')

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000)
fin=fin+1
fig=p1.figure(sfin) 


"""


fig01 = drawmatrix_channels(coh, roi_names, size=[10., 10.], color_anchor=0)
"""

fig01.savefig('sgc_matrix..png')

mbin2=np.array(mbin) 
mbin2=nx.to_networkx_graph(mbin2,create_using=nx.DiGraph())   #directed graph. 




"""		 
im=p1.imshow(WeightsAb,interpolation='nearest') 

p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Covariance of weight before training') 
# p1.autoscale(True) 
p1.grid(True)  """

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(sfin) 
fin=fin+1

"""
im=p1.imshow(nTEin,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('The difference between weights before and after training') 
# p1.autoscale(True) 
p1.grid(True) 
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(sfin) 
fin=fin+1
"""
"""
im=p1.imshow(nTE,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('nTE after running') 
p1.grid(True) 
#  p1.autoscale(True) 
#}
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(sfin) 
fin=fin+1




im=p1.imshow(mbin,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Adjacency matrix both transmitters') 
p1.autoscale(True) 
p1.grid(True) 
# p1.figure(5) 
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(sfin) 
fin=fin+1


im=p1.imshow(Matrix, norm=LogNorm(),interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix Excitatory and Inhibitory Connections') 
p1.autoscale(True) 
p1.grid(True) 
# p1 .grid(b='on',which='minor')   ##, which='major', axis='both', **kwargs)	
# p1.ax.grid(color='r',linewidth=1) 

# im=p1.imshow(Matrix)
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(sfin) 
fin=fin+1
# im=p1.imshow(MatrixG) 
im=p1.imshow(MatrixG, norm=LogNorm(),interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix GABA') 
p1.grid(True) 
p1.autoscale(True) 

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(sfin) 
fin=fin+1
# im=p1.imshow(MatrixA) 
im=p1.imshow(MatrixA, norm=LogNorm(),interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix AMPA') 
p1.grid(True) 
p1.autoscale(True) 


sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(sfin) 
fin=fin+1
#import networkx as nx 
#print ma
#Ga=nx.from_numpy_matrix(ma) 
nx.draw_networkx(Ga,label='Structural connectivity, AMPA, degree n') 
p1.title('Structural connectivity, AMPA, degree n') 
p1.draw()

#nx.draw_networkx_labels(Ga,'Structural connectivity, AMPA, degree n') 


sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(sfin) 
fin=fin+1
#   import networkx as nx 
Gg=nx.from_numpy_matrix(mg) 
nx.draw_networkx(Gg,label='Structural connectivity, GABA, degree n') 

p1.title('Structural connectivity, GABA, degree n')
p1.draw()
#p1.show()
# nx.draw_networkx(Gg) 
# p1.title('Structural connectivity, GABA, degree n') 
#sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(sfin) 
#fin=fin+1
#   import networkx as nx 
# Gg=nx.from_numpy_matrix(mg) 

#   draw_networkx_labels(Ga,'Structural connectivity, AMPA, degree n') 
#p1.draw()

#nx.draw_networkx_labels(mbin2,label='Structural connectivity,degree 1') 
# # # #/
#Show for every graph.
# #/




# p1.show() 


ps='pruned by= '+str(h.prunenet)+' (ums) '+ 'feed forward (true/false)= '+str(h.ff)
# pos=nx.spring_layout(mbin) # positions for all nodes 
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(sfin) 
fin=fin+1

#mbin2=nx.from_numpy_matrix(mbin) 
nx.draw_networkx(mbin2,label='Structural connectivity, degree 1') 
p1.title(ps+'Structural connectivity, degree 1') 
p1.draw() 
p1.savefig('Structural_connectivity_degree 1..png')



sfin=str(h.prunenet)+str(fin)+'.png' 
#sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin,dpi=1000) 
fig=p1.figure(sfin) 
fin=fin+1
"""

