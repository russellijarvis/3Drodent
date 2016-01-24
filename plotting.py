
mbin2=np.array(mbin) 
mbin2=nx.to_networkx_graph(mbin2,create_using=nx.DiGraph())   #directed graph. 



fin=0 
from matplotlib.colors import LogNorm 
import matplotlib 
matplotlib.use('Agg') 
import pylab as p1
import pylab as plt
from scipy import signal 
#and user make sure that param6.hoc has run.
# fig=p1.figure(fin) 
#fin=fin+1

"""		 im=p1.imshow(WeightsAb,interpolation='nearest') 

p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Covariance of weight before training') 
# p1.autoscale(True) 
p1.grid(True)  """

fig=p1.figure(fin) 
fin=fin+1

"""
im=p1.imshow(nTEin,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('The difference between weights before and after training') 
# p1.autoscale(True) 
p1.grid(True) 
fig=p1.figure(fin) 
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
fig=p1.figure(fin) 
fin=fin+1




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
fig=p1.figure(fin) 
fin=fin+1
# im=p1.imshow(MatrixG) 
im=p1.imshow(MatrixG, norm=LogNorm(),interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix GABA') 
p1.grid(True) 
p1.autoscale(True) 

fig=p1.figure(fin) 
fin=fin+1
# im=p1.imshow(MatrixA) 
im=p1.imshow(MatrixA, norm=LogNorm(),interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Degree matrix AMPA') 
p1.grid(True) 
p1.autoscale(True) 


fig=p1.figure(fin) 
fin=fin+1
#import networkx as nx 
#print ma
#Ga=nx.from_numpy_matrix(ma) 
nx.draw_networkx(Ga,label='Structural connectivity, AMPA, degree n') 
plt.title('Structural connectivity, AMPA, degree n') 
plt.draw()

#nx.draw_networkx_labels(Ga,'Structural connectivity, AMPA, degree n') 


fig=p1.figure(fin) 
fin=fin+1
#   import networkx as nx 
Gg=nx.from_numpy_matrix(mg) 
nx.draw_networkx(Gg,label='Structural connectivity, GABA, degree n') 

p1.title('Structural connectivity, GABA, degree n')
p1.draw()
#p1.show()
# nx.draw_networkx(Gg) 
# plt.title('Structural connectivity, GABA, degree n') 
#fig=p1.figure(fin) 
#fin=fin+1
#   import networkx as nx 
# Gg=nx.from_numpy_matrix(mg) 

#   draw_networkx_labels(Ga,'Structural connectivity, AMPA, degree n') 
#plt.draw()

#nx.draw_networkx_labels(mbin2,label='Structural connectivity,degree 1') 
# # # #/
#Show for every graph.
# #/




# p1.show() 
"""

ps='pruned by= '+str(h.prunenet)+' (ums) '+ 'feed forward (true/false)= '+str(h.ff)
# pos=nx.spring_layout(mbin) # positions for all nodes 
fig=fig=p1.figure(fin) 
fin=fin+1

#mbin2=nx.from_numpy_matrix(mbin) 
nx.draw_networkx(mbin2,label='Structural connectivity, degree 1') 
plt.title(ps+'Structural connectivity, degree 1') 
plt.draw() 

fig.savefig('Structural_connectivity_degree_1.png')



#proc disp_graphs(){
# fig=p1.figure(fin) 
fig=p1.figure(fin) 
fin=fin+1

pos=nx.spring_layout(dir2) # positions for all nodes 
nx.draw_networkx_nodes(dir2,pos,node_size=700) 
nx.draw_networkx_edges(dir2,pos,width=2.5) 
nx.draw_networkx_labels(dir2,pos,font_size=20,font_family='sans-serif')
nx.draw_networkx(dir2,label=ps+'Effective connectivity via pref direction, degree 1') 
p1.title(ps+'Effective connectivity, via pref direction, degree 1') 
p1.draw() 
fig=fig=p1.figure(fin)
fig.savefig('Effective_connectivity_degree.png')

# plt.show() 

"""

fig=p1.figure(fin) 
fin=fin+1
im=p1.imshow(entfa,interpolation='nearest') 
im=p1.grid(b='True',color='white',which='both', axis='both') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
ff=h.ff
prunenet=h.prunenet
title_str='Between cell nTE across run'+'ff='+str(ff)+'prunenet='+str(prunenet) 
p1.grid(True) 
p1.title(title_str) 
#

fig=p1.figure(fin) 
fin=fin+1
im=p1.imshow(entfa, norm=LogNorm(),interpolation='nearest') 
p1.autoscale(True) 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 

prunenet=prunenet
title_str='scaled Between cell nTE across run'+'ff='+str(ff)+'prunenet='+str(prunenet) 
p1.grid(True) 
p1.title(title_str) 

#print idvec.max*idvec.max, cells.count*cells.count, idvec.max*idvec.max
fig=p1.figure(fin) 
fin=fin+1
im=p1.imshow(corr,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
ff=ff
prunenet=prunenet
title_str='Between cell correlations for run'+'ff='+str(ff)+'prunenet='+str(prunenet) 
p1.title(title_str) 

fig=p1.figure(fin) 
fin=fin+1
im=p1.imshow(dir, norm=LogNorm(),interpolation='nearest') 
# p1.grid(True) 

p1.grid(b='True',color='black',which='both',axis='both',linestyle='--', linewidth=0.4) 
p1.autoscale(True) 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
ff=ff
prunenet=prunenet
title_str='Between cell prefered direction for run'+'ff='+str(ff)+'prunenet='+str(prunenet) 
p1.title(title_str) 

p1.show() 
"""
entfa=np.array(entfa) 

""" efc = nx.DiGraph() 
efc.add_nodes_from(range(0,numcell+2)) 
efc.add_edges_from(entfa)  """
entfa2=nx.to_networkx_graph(entfa,create_using=nx.DiGraph())   #directed graph. 
# nx.add_weighted_edges_from(dir) 
#nx.draw_networkx(entfa2) 
nx.draw_networkx(entfa2,label='Effective Connectivity') 
#plt.title('Effective Connectivity') 
plt.draw()

# this is a linear filter, I probably want gaussian smoothing, and butterworth bandpass.

"""
theta1 = signal.lfilter(0, 30, float(ecpv1))
gamma1 = signal.lfilter(0.80, 0.100, ecpv)


theta2 = signal.lfilter(0, 30, ecpv2)
gamma2 = signal.lfilter(80, 100, ecpv2)


fin+=1
fig=p1.figure(fin)
p1.hold(True)
p1.plot(theta1)
p1.plot(gamma1)
p1.plot(theta2)
p1.plot(gamma2)
p1.hold(False)
p1.show()  

if(h.ff!=1):
	outout=downsample(time_courses[int(outdegree)],oldrate=40000,newrate=200)
	outin=downsample(time_courses[int(indegree)],oldrate=40000,newrate=200)
	synapse=downsample(recinsyn[int(indegree)],oldrate=40000,newrate=200)
	F,pp,cohe,Fx2y,Fy2x,Fxy=granger(outout,outin,20)#,200,100)# nfs2[0], nfs2[1] ,
	fin+=1
	fig=plt.figure(fin) 
	# Plot Granger spectra
	plt.hold(True) #Matlab style hold
	labels=list()
	labels=['Fx2y','Fy2x','Fxy','cohe','pp']
	alldata=list()
	colors=[[0,0,1],[1,0,0],[0,0.5,0.5],[0.5,1,0],[1,0.5,1]]
	lb1='causality between cells i= '+str(int(outdegree)) +' and cell j= '+str(int(indexj))
	lb2='causality between cells i= '+str(int(indegree)) +' and cell j= '+str(int(indexi))
	plt.plot(F,Fxy,label=lb1,linewidth=1.5,c=colors[0])
	plt.plot(F,Fy2x,label=lb2,linewidth=1.5,c=colors[1])
	#plt.xlim(0,2)
	#plt.ylim(-0.4,0.4)
	plt.xlabel("Hz")
	plt.xlim(0,100)
	plt.ylim(-2,2)
	plt.title('GSC components'+bc)
	plt.legend()#bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
	#      ncol=2, mode="expand", borderaxespad=0.)
	plt.hold(False)

	F,pp,cohe,Fx2y,Fy2x,Fxy=granger(outout,synapse,20)#,200,100)# nfs2[0], nfs2[1] ,
	fig=plt.figure(fin) 
	# Plot Granger spectra
	plt.hold(True) #Matlab style hold
	labels=list()
	labels=['Fx2y','Fy2x','Fxy','cohe','pp']
	alldata=list()
	colors=[[0,0,1],[1,0,0],[0,0.5,0.5],[0.5,1,0],[1,0.5,1]]
	lb1='causality between cells i= '+str(int(outdegree)) +' and cell j= '+str(int(indexj))
	lb2='causality between cells i= '+str(int(indegree)) +' and cell j= '+str(int(indexi))
	plt.plot(F,Fxy,label=lb1,linewidth=1.5,c=colors[0])
	plt.plot(F,Fy2x,label=lb2,linewidth=1.5,c=colors[1])
	#plt.xlim(0,2)
	#plt.ylim(-0.4,0.4)
	plt.xlabel("Hz")
	plt.xlim(0,100)
	plt.ylim(-2,2)
	plt.title('GSC components'+bc)
	plt.legend()#bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
	#      ncol=2, mode="expand", borderaxespad=0.)
	plt.hold(False)

#F,pp,cohe,Fx2y,Fy2x,Fxy=granger(outout,synapse,20)#,200,100)# nfs2[0], nfs2[1] ,
plt.show() 
"""	
