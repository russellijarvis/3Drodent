
#Now that I have limited the number of connections these colours for the degree matrices are realistic.

#import pyhoc    
#import statsmodels as statsmodels
from matplotlib.colors import LogNorm 
#import matplotlib
#matplotlib.use('Agg') 




#nrnpython("medr=np.arange(0,int(numcell),0)")
#nrnpython("np.divide(input_marke,rates,medr)")
#nrnpython("meanedr=np.mean(medr)")
#nrnpython("print meanedr, 'mean entropy divided by firing ratem, this should be low for FF, high for FB')


fin=0
p1.figure(fin)
p1.clf()


nx.draw(structuralg,node_size=5,node_color='r',edge_color='b',alpha=.2)

sfin='Structural Graph of the Network'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(numcell)+str(int(fin))+'.png'
p1.savefig(sfin)


out_degrees = stucturalg.out_degree()
# dictionary node:degree
out_values = sorted(set(out_degrees.values()))
out_hist = [out_degrees.values().count(x) for x in out_values]


in_degrees = stucturalg.in_degree()
# dictionary node:degree
in_values = sorted(set(in_degrees.values()))
in_hist = [in_degrees.values().count(x) for x in in_values]
fin+=1
p1.figure(fin)
p1.clf()
p1.hold(True)
p1.plot(in_values,in_hist,'ro-')
# in-degree
p1.plot(out_values,out_hist,'bv-')
p1.xscale('log')
p1.yscale('log')
p1.hold(False)
# out-degree
p1.legend(['In-degree','Out-degree'])
p1.xlabel('Degree Scale free network')
p1.ylabel('Number of Neurons that have this degree value')
p1.title('Degree Scale free network')
sfin='dirg_degree_distribution.png'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(numcell)+str(int(fin))+'.png'
p1.savefig(sfin)



fin=0
h('print nclist.count, "nclist.count", minr, "minr"')

fin+=1
fig=p1.figure(fin)
fig.clf()
if int(h.get_dist)==1:
  #h('py.the_length=py.int(h.adjnqs.v[8].size)'); across=np.arange(0,the_length,0)
  h('py.dist=adjnqs.v[7].to_python()')
  h('py.transfer=adjnqs.v[8].to_python()')
  #p1.plot(across,transfer)
  #p1.scatter(transfer,dist)
  xdata=transfer
  ydata=dist

  p1.scatter(xdata,ydata)
  p1.hold(True)
  slope, yint = p1.polyfit(xdata,ydata,1)
  xline = p1.xticks()[0]
  yline = map(lambda x: slope*x+yint,xline)
  p1.plot(xline,yline,ls='--',color='b')
  # Set new x- and y-axis limits
  p1.xlim((0.0,max(xdata)+(.15*max(xdata))))
  p1.ylim((0.0,max(ydata)+(.15*max(ydata))))
  
  p1.title('arc length versus transfer impedence at synapse positions 100Hz')
  p1.xlabel("um")
  p1.ylabel("Mohm")
  sfin='arc length '+str(int(h.prunenet))+str(tstop)+str (h.plastic)+str(int(h.ff))+str(numcell)+str(int(fin))+'.png'
  p1.savefig(sfin, format='png') 




fin+=1
fig=p1.figure(fin)
fig.clf()
if int(h.get_dist)==1:

  h('py.targets=adjnqs.v[0].to_python()')
  targets=np.array(targets)
  transfer=np.array(transfer)
  dist=np.array(dist)
  np.where(targets==indegree)
  #h('py.the_length=py.int(h.adjnqs.v[8].size)'); across=np.arange(0,the_length,0)

  #p1.plot(across,transfer)
  xdata=transfer[np.where(targets==indegree)]
  ydata=dist[np.where(targets==indegree)]
  #p1.scatter(transfer[np.where(targets==indegree)],dist[np.where(targets==indegree)])
  p1.scatter(xdata,ydata)
  p1.hold(True)
  slope, yint = p1.polyfit(xdata,ydata,1)
  xline = p1.xticks()[0]
  yline = map(lambda x: slope*x+yint,xline)
  p1.plot(xline,yline,ls='--',color='b')
  # Set new x- and y-axis limits
  p1.xlim((0.0,max(xdata)+(.15*max(xdata))))
  p1.ylim((0.0,max(ydata)+(.15*max(ydata))))
  p1.title('arc length versus transfer impedence at synapse positions 100Hz')
  p1.xlabel("um")
  p1.ylabel("Mohm")
  sfin='arc length versus distance but only on indegree cell'+str(int(h.prunenet))+str(tstop)+str (h.plastic)+str(int(h.ff))+str(numcell)+str(int(fin))+'.png'
  p1.savefig(sfin, format='png') 
  
  #h.


fin=fin+1
fig=p1.figure(fin) 
fig.clf()
im=p1.imshow(mbin,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Adjacency matrix both transmitters') 
p1.autoscale(True) 
p1.grid(True) 
# p1.figure(5) 
fig=p1.figure(fin) 
fin=fin+1
sfin='Both_transmitters_adjacency_matrix'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(numcell)+str(int(fin))+'.png' 
fig.savefig(sfin)

fin=fin+1
fig=p1.figure(fin) 
fig.clf()
im=p1.imshow(Matrix, interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix Excitatory and Inhibitory Connections') 
p1.autoscale(True) 
p1.grid(True) 
sfin='Both_transmitters_degree_matrix'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(numcell)+str(int(fin))+'.png' 
fig.savefig(sfin)
# p1 .grid(b='on',which='minor')   ##, which='major', axis='both', **kwargs)	
# p1.ax.grid(color='r',linewidth=1) 

# im=p1.imshow(Matrix)
fin=fin+1
fig=p1.figure(fin) 
fig.clf()
# im=p1.imshow(MatrixG) 
im=p1.imshow(MatrixG, interpolation='nearest') 
p1.autoscale(True) 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix GABA') 
p1.grid(True) 
sfin='GABA_degree_matrix'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(numcell)+str(int(fin))+'.png' 
fig.savefig(sfin)

fin=fin+1
fig=p1.figure(fin) 
fig.clf()
# im=p1.imshow(MatrixA) 
im=p1.imshow(MatrixA, interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix AMPA') 
p1.grid(True) 
p1.autoscale(True) 
sfin='AMPA_degree_matrix'+str(int(h.prunenet))+str(tstop)+str(h.plastic)+str(int(h.ff))+str(numcell)+str(int(fin))+'.png' 
fig.savefig(sfin)




#j=len(colors)-1
numcell=numcell
numcell=numcell
roi_names=[0 for x in xrange(0,int(numcell))]
for i in xrange(0,int(numcell-1)):
 roi_names[i]=str(i)+' '+str(h.cells.o(i).nametype)
 print roi_names[i]
 
#In the future it would be nice to label every xtick with 
#the elements of roi_names

inimp0=[]
inimp20=[]
inimp100=[]
h('for i=0,cells.count-1{ py.inimp0.append(cells.o(i).inimp0) }')
h('for i=0,cells.count-1{ py.inimp20.append(cells.o(i).inimp20) }')
h('for i=0,cells.count-1{ py.inimp100.append(cells.o(i).inimp100) }')

across=np.arange(0,int(numcell),1)
#p1.plot(tvec,intervec,'bs')#p1.plot(tvec,pyra,'g^')
fin+=1

p1.figure(fin)
p1.hold(True)

p1.plot(across,np.array(inimp0),'ro',linewidth=1,label='input impedence at DC')
p1.plot(across,np.array(inimp20),'bo',linewidth=1,label='input impedence at 20Hz')#,c=colors[j], markeredgecolor = 'none')
j-=1
p1.plot(across,np.array(inimp100),'go',linewidth=1, label='input impedence at 100Hz')
p1.title("Input Impedence Versus cell")
p1.xlabel('cell number')
p1.ylabel('MOhm')
p1.hold(False)
#p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8, ncol=2, mode="expand", borderaxespad=0.), label='input impedence at 20Hz')
sfin='Input Impedence Versus cell'+str(int(h.prunenet))+str(int(h.ff))+str(int(fin))
p1.savefig(sfin)

