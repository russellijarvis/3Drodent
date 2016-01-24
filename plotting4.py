import nitime
import nitime.analysis as nta
import nitime.timeseries as ts
import nitime.utils as tsu
from nitime.viz import drawmatrix_channels
#from nitime.viz import drawmatrix_channels2
from nitime.viz import drawgraph_channels
from nitime.viz import draw_graph

import statsmodels as statsmodels
from matplotlib.colors import LogNorm 
import matplotlib 
#matplotlib.use('Agg') 
import pylab as p2
import pyhoc    
from scipy import signal 
from bsmart import granger # Load the Granger calculation tool
#and user make sure that param6.hoc has run.
# sfin=str(h.prunenet)+str(fin)+'.png' 
#fig.savefig(sfin) 

evr=[]
for i in xrange(0, len(rates)-1):
 evr.append(rates[i]*tstop*np.log2(e/(rates[i]*0.01)))

fin=0
fig=p1.figure(fin)
fig.clf()
#fig = drawmatrix_channelsnn(np.matrix(entfa), np.array(roi_names), size= [10., 10.],color_anchor=0)
im=p1.imshow(entfa,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('transfer entropy') 
# p1.autoscale(True) 
p1.grid(True) 
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1

#p1.show()


    
fin+=1
fig = p1.figure(fin)
#p1.clf()
#nx.draw_networkx(dirfinal,nodelist=dicten,edgelist=dictee)#, roi_names)
#nx.draw(dirfinal)
#nx.draw(poot)
#fig = drawgraph_channels(dirfinal, roi_names)
#sfin='Transfer Entropy graph centrality 
#sfin='now 2'+str(h.prunenet)+str(fin)+'.png' 
#fig.savefig(sfin)

execfile('pyhoc.py')

""" 
dictc=nx.degree_centrality(dirfinals2)
#nx.draw(dirfinals2[dictc[0,10]])#, roi_names)
    
fin+=1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(dirfinals,nodelist=dictc)#, roi_names)
sfin='subset'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)

dicten=nx.betweenness_centrality(dirfinal)
dictee=nx.edge_betweenness_centrality(dirfinal)
#nx.draw(dirfinals2[dictc[0,10]])#, roi_names)
listout=sorted(dicte.items(), key=itemgetter(1))
l2=listout[len(listout)-15:]
"""
from bsmart import timefreq, pwcausalr
from scipy import array, size
import pyhoc    
n_samples = int(tstop/h.dt)+2 #data_rec.shape[0]
nseq=int(numcell)

roi_names=[0 for x in xrange(0,int(numcell))]
for i in xrange(0,int(numcell-1)):
 roi_names[i]=str(i)+' '+str(h.cells.o(i).nametype)
 print roi_names[i]
#Make an empty container for the data
msgcv= [[0 for x in xrange(0,int(int(int(numcell))))] for x in xrange(0,int(int(int(numcell))))]
msgcv2= [[0 for x in xrange(0,int(int(int(numcell))))] for x in xrange(0,int(int(int(numcell))))]

#cant even really use a basic example, because if its full of nans it wont evaluate subsequentely.


#  for i=0,cells.count-1{//count up
#    for(j=cells.count-1;j>0;j-=1){//count down 
for indexi in xrange(0,int(numcell-1)):
#  for indexj in xrange(0,int(numcell-1)):
  indexj=0
  if(indexi!=indexj):    
    maxt=h.times[indexj].max
    #printf("maxt=%g\n",maxt)
    if(h.times[indexi].max>maxt): maxt=h.times[indexi].max
    #printf("maxt=%g\n",maxt)
    binsz=10
    #if(len(h.times[i].to_python())>0):
    #if(len(h.times[j].to_python())>0):
    #if(maxt>0):
from bsmart import timefreq, pwcausalr
from scipy import array, size
import pyhoc  
def get_sgc(indexi,indexj):
    vec1=downsample(time_courses[int(indexi)],oldrate=40000,newrate=200)
    vec2=downsample(time_courses[int(indexj)],oldrate=40000,newrate=200)
        
    """
    order=10
    rate=200
    maxfreq=0
    
    if maxfreq==0: F=timefreq(outj,rate) # Define the frequency points
    else: F=array(range(0,maxfreq+1)) # Or just pick them
    npts=size(F,0)

   
    h("for k=0,1 vb[k] = new Vector() //vb stands for vector binned.")
    h("vb[0].hist(times[py.indexi],0,(py.maxt+py.binsz-1)/py.binsz,py.binsz)")
    h("vb[1].hist(times[py.indexj],0,(py.maxt+py.binsz-1)/py.binsz,py.binsz)")

    h("vo1=normte(vb[0],vb[1],20)")
    #return F,pp[0,:],cohe[0,:],Fx2y[0,:],Fy2x[0,:],Fxy[0,:]
    
    """

    npts=len(vec1)
    fs=200
    n=len(vec1)
    freq=100
    p=15
    ntrls=1
                          #(x,ntrls,npts,p,fs,freq)
    x=array([vec1,vec2])
    F,pp,cohe,Fx2y,Fy2x,Fxy=pwcausalr(x,1,npts,order,rate,maxfreq)
    print "GSC components", np.mean(Fx2y), np.mean(Fy2x), 
    print "GSC difference", np.mean(Fx2y)-np.mean(Fy2x)
    print indexi
    #print h.ntefa.v[2].x[indexi]
    #F,pp,cohe,Fx2y,Fy2x,Fxy=pwcausalr(np.array(x),ntrls,npts,p,fs,freq)                      
#F,pp,cohe,Fx2y,Fy2x,Fxy=granger(x,ntrls,npts,p,fs,freq)#,200,100)# nfs2[0], nfs2[1] ,
    #F,pp,cohe,Fx2y,Fy2x,Fxy=granger(time_courses[int(indexj)],time_courses[int(indexi)],20)#,200,100)# nfs2[0], nfs2[1] ,
    #the two refers to a projection of x onto y, not an order of 2. Fill a matrix with
        #These values. 

    if(~np.isnan(np.mean(Fx2y))):
      msgcv2[indexi][indexj]=(np.mean(Fx2y)-np.mean(Fy2x))
      msgcv[indexi][indexj]=np.mean(Fx2y)
#Since granger values are x->y and are not constrained, it's more informative to look at the values for x -> y minus y -> x (to reduce the confound of simultaneous influence) We'll create that difference matrix and plot it as well.
    else:
      msgcv[indexi][indexj]=0

      
  #print 'pacificier', indexi 
  
for indexi in xrange(0,int(numcell-1)):
  for indexj in xrange(0,int(numcell-1)):
  

     if(~np.isnan(msgcv[indexi][indexj])):
       msgcv[indexi][indexj]= msgcv[indexi][indexj]
     else:
       msgcv[indexi][indexj]=0
     if(~np.isnan(msgcv[indexi][indexj])):
       msgcv2[indexi][indexj]= msgcv[indexi][indexj]
     else:
       msgcv2[indexi][indexj]=0
   
trains[0:]   
prefdire=[]
prefdirs=[]

mx=msgcv2
low_values_indices = np.mean(msgcv)> 0#///  # Where values are low
mx[low_values_indices] = 0#





fin+1
fig=p1.figure(fin)
fig.clf()
im=p1.imshow(entfa,interpolation='nearest') 
p1.colorbar(im) 
#p1.autoscale(True) 
#p1.grid(True) 
p1.grid(True) 
p1.imshow(dir, interpolation='nearest')
sfin='x4_pref_dir_imshow'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)   

   
""" 
   
   
   
   
   
def cdc(graph):
 g=graph()
 dc=nx.degree_centrality(g)
 nx.set_node_attributes(g,'degree_cent',dc)
 dcsorted=sorted(dc.items(),key=itemgetter(1),reverse=True)
 for key,value in degcent_sorted[0:10]:
  print "highest degree", key, value
 return graph, dcsorted

listb=[]
for line in nx.generate_edgelist(msgc3,data=['weight']):
#if (data>0):
 print line
listb.append(line)
#data=['weight'])

def highest_centrality(cent_dict):
# Create ordered tuple of centrality data
cent_items=[(b,a) for (a,b) in cent_dict.iteritems()]
# Sort in descending order
cent_items.sort()
#cent_items.reverse()
return tuple((cent_items[:10]))

nodelist=highest_centrality()

# p1 .grid(b='on',which='minor')   ##, which='major', axis='both', **kwargs)	
# p1.ax.grid(color='r',linewidth=1)    

#First pref dir NTE
fin+=1
p1.clf()
fig.clf()
fig=p1.figure(sfin)
dir_=nx.to_networkx_graph(np.matrix(dir),create_using=nx.DiGraph())
#dir_=nx.to_networkx_graph(np.matrix(dir),create_using=nx.DiGraph())
bet_c = nx.degree_centrality(dir_)
nodelist=bet_d = highest_centrality(bet_c)

pos,nx.spring_layout(dir_)
nx.draw_networkx(dir_,pos=pos,node_list=nodelist)#,label=ps+'Effective connectivity via pref direction, degree 1') 
p1.title('Effective connectivity, via pref direction, degree 1') 
p1.draw() 
sfin='1_'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)#+sfin) 




#First pref dir sgc
fin+=1
msgc3=nx.to_networkx_graph(np.matrix(msgcv2),create_using=nx.DiGraph())
p1.clf()
fig.clf()
fig=p1.figure(sfin)
#msg_=nx.to_networkx_graph(np.matrix(dir),create_using=nx.DiGraph())
bet_msg = nx.degree_centrality(msgc3)
bet_ms = highest_centrality(bet_msg)


for i in xrange(0,59):
 if(sum(pos[i])<1): print i
#[if(sum(pos[v])<1): for v in (msgc3)]

sub=msgc3.subgraph(bet_ms)


H=nx.Graph()
for v in msgc3:
 H.add_node(v)
 for (u,v,d) in msgc3.edges(data=True):
  if d['weight'] > 0.1:
   H.add_edge(u,v)

# with nodes colored by degree sized by population
node_color=[msgc3.degree(v) for v in msgc3]
node_size=[msgc3.degree(v)*20 for v in msgc3]
pos=nx.spring_layout(msgc3)

nx.draw_networkx(msgc3,pos=pos,node_size=node_size,node_color=node_color)

# scale the axes equally
plt.savefig("knuth_miles.png")


#nx.draw_networkx(msgc3)#,label=ps+'Effective connectivity via pref direction, degree 1') 
pos=nx.spring_layout(sub)
nx.draw_networkx_(msgc3,nodelist=bet_ms)
nx.draw_networkx_edges(msgc3,pos,edgelist=bet_ms)
# specifiy edge labels explicitly
edge_labels=dict([((u,v,),d['weight'])
             for u,v,d in G.edges(data=True)])
nx.draw_networkx_edge_labels(msgc3,pos,edge_labels=edge_labels)




ncenter, _ = min(pos.items(), key = lambda (node, (x, y)): (x-2)**2+(y-2)**2)

# color by path length from node near center
p = {node:length
     for node, length in nx.length(G, ncenter).items()
      if length =1
    }

plt.figure(figsize = (8, 8))
node_color = p.values()
H = G.subgraph(p.keys())    
nx.draw_networkx_edges(H, pos, alpha = 0.4)
nx.draw_networkx_nodes(H, pos, node_size = 80, node_color = node_color)
nx.draw_networkx_edge_labels(H,pos,edge_labels=edge_labels)




plt.xlim(-0.05, 1.05)
plt.ylim(-0.05, 1.05)
plt.axis('off')
plt.savefig('random_geometric_graph.png')
#plt.show()


pos = nx.get_node_attributes(msgc3, 'pos')
pos=nx.spring_layout(msgc3)
# find node near center (0.5,0.5)
ncenter, _ = min(pos.items(), key = lambda (node, (x, y)): (x-0.5)**2+(y-0.5)**2)

# color by path length from node near center
p = {length
     for node, length in nx.length(G, ncenter).items()
     if length<2:
    }

plt.clf()
plt.figure(figsize = (8, 8))
node_color = p.values()
H = msgc3.subgraph(p.keys()) 
#H=nx.subgraph(nx.strongly_connected_components(msgc3))   
nx.draw_networkx_edges(H, pos, alpha = 0.4)
nx.draw_networkx_nodes(H, pos, node_size = 80, node_color='#0F1C95',
                   cmap=plt.cm.Reds_r)
nx.draw_networkx_edge_labels(H,pos,edge_labels=edge_labels)
plt.savefig('random_geometric_graph.png')


"""




#p1.title(ps+'Effective connectivity, via pref direction, degree 1') 
p1.draw() 
sfin='2_'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig('working')#+sfin) 
#Third imshow

#fig.savefig('Effective_connectivity_degree_1.png')




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
fig.savefig(sfin)

  
#msgcv=np.array(msgcv) 
#low_values_indices = np.mean(msgcv)+np.std(msgcv) < 0#///  # Where values are low
#msgcv[low_values_indices] = 0#//  # All low values set to 0





dirfinals =nx.DiGraph() 
dirfinals2 =nx.DiGraph() 
#dirfinals.add_nodes_from(msgcv) 
for i in xrange(0,len(msgcv)):
 for j in xrange(0,len(msgcv)):
   if (int(h.prunenet>0)):
     if msgcv[i][j]>(np.mean(msgcv)+2.5*np.std(msgcv)):
       dirfinals.add_edge(i,j,weight=msgcv[i][j]) 
     if msgcv2[i][j]>(np.mean(msgcv)+2.5*np.std(msgcv)):
       dirfinals2.add_edge(i,j,weight=msgcv2[i][j]) 
   else:
     if msgcv[i][j]>(np.mean(msgcv)+1.5*np.std(msgcv)):
       dirfinals.add_edge(i,j,weight=msgcv[i][j]) 
     if msgcv2[i][j]>(np.mean(msgcv)+1.5*np.std(msgcv)):
       dirfinals2.add_edge(i,j,weight=msgcv2[i][j]) 
    
    
dictc=nx.degree_centrality(dirfinals2)
#nx.dr aw(dirfinals2[dictc[0,10]])#, roi_names)
    
fin+=1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(dirfinals)#,nodelist=dictc)#, roi_names)
sfin='nowsgc'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)
#fig05 = draw_graph(dirfinals)
#fig05.savefig('experiment.png')




fin+=1
fig = p1.figure(fin)
nx.draw(dirfinals2)#, roi_names)
p1.draw() 




fin+=1
fig = p1.figure(fin)
nx.draw(dir2)#, roi_names)
p1.draw() 

#fig.savefig(sfin) 
fig=p1.figure(sfin)
fig.savefig('core Effective_connectivity_degree_1'+sfin)
#s = g.subgraph(nx.shortest_path(g.to_undirected(),55))
#dir2.subgraph(
#core=nx.strongly_connected_component_subgraphs(dir2)

#for h in nx.connected_component_subgraphs(dir2.to_undirected()):
#   nx.draw(h,pos,node_color='red')
 
#if target_node in h:
# else:
#  nx.draw(h,pos,node_color='white')

#plt.show()
#pos=nx.spring_layout(core) # positions for all nodes 
#nx.draw_networkx_nodes(core,pos,node_size=700) 
#nx.draw_networkx_edges(core,pos,width=2.5) 
#nx.draw_networkx_labels(core,pos,font_size=20,font_family='sans-serif')

#nx.draw_networkx(dir2,nodes=nx.strongly_connected_component_subgraphs(dir2),label=sfin+'Effective connectivity via pref direction, degree 1') 
#sfin=str(h.prunenet)+str(fin)+'.png' 
#p1.title(sfin+'Effective connectivity, via pref direction, degree 1') 





sfin='granger_graph_final_withsubtraction'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)


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
fig.savefig(sfin)   
fin+=1
spc=0 # Only one column, so pick it

	





sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
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
fig.savefig(sfin) 
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
fig.savefig(sfin) 
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
fig.savefig(sfin) 
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
fig.savefig(sfin) 
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
p1.plot(across,divout,'go',label='H(x) of cell num spike train')
p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
       ncol=2, mode="expand", borderaxespad=0.)

maxtime=int(h.tstop)

p1.xlabel("Neuron number")
p1.ylabel("Bits/sec")

###
p1.hold(False)
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
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
fig.savefig(sfin)
fin=fin+1
fig=p1.figure(sfin) 



"""


fig01 = drawmatrix_channels(coh, roi_names, size=[10., 10.], color_anchor=0)
"""

#fig01.savefig('sgc_matrix..png')

#mbin2=np.array(mbin) 
#mbin2=nx.to_networkx_graph(mbin2,create_using=nx.DiGraph())   #directed graph. 




"""		 
im=p1.imshow(WeightsAb,interpolation='nearest') 

p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('Covariance of weight before training') 
# p1.autoscale(True) 
p1.grid(True)  """

#sfin=str(h.prunenet)+str(fin)+'.png' 
#fig.savefig(sfin) 
#fig=p1.figure(sfin) 
#fin=fin+1

"""
im=p1.imshow(nTEin,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
p1.title('The difference between weights before and after training') 
# p1.autoscale(True) 
p1.grid(True) 
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
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
fig.savefig(sfin) 
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
fig.savefig(sfin) 
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
fig.savefig(sfin) 
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
fig.savefig(sfin) 
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
fig.savefig(sfin) 
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
fig.savefig(sfin) 
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
fig.savefig(sfin) 
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
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1

#mbin2=nx.from_numpy_matrix(mbin) 
nx.draw_networkx(mbin2,label='Structural connectivity, degree 1') 
p1.title(ps+'Structural connectivity, degree 1') 
p1.draw() 
p1.savefig('Structural_connectivity_degree 1..png')



sfin=str(h.prunenet)+str(fin)+'.png' 
#sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1


#execfile('gct.py')

"""

# p1.show() 

"""

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
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

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
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
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
fin=fin+1
im=p1.imshow(corr,interpolation='nearest') 
p1.colorbar(im) 
p1.xlabel('columns = targets') 
p1.ylabel('rows = sources') 
ff=ff
prunenet=prunenet
title_str='Between cell correlations for run'+'ff='+str(ff)+'prunenet='+str(prunenet) 
p1.title(title_str) 

sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin) 
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

 
""" efc = nx.DiGraph() 
efc.add_nodes_from(range(0,int(numcell)+2)) 
efc.add_edges_from(entfa)  """
#entfa2=nx.to_networkx_graph(entfa,create_using=nx.DiGraph())   #directed graph. 

#nx.draw_networkx(entfa2,label='significance based Effective Connectivity') 
#p1.title('Effective Connectivity') 
#p1.draw()

# this is a linear filter, I probably want gaussian smoothing, and butterworth bandpass.

"""
theta1 = signal.lfilter(0, 30, float(ecpv1))
gamma1 = signal.lfilter(0.80, 0.100, ecpv)


theta2 = signal.lfilter(0, 30, ecpv2)
gamma2 = signal.lfilter(80, 100, ecpv2)


fin+=1
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(sfin)
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
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin) 
# Plot Granger spectra
p1.hold(True) #Matlab style hold
labels=list()
labels=['Fx2y','Fy2x','Fxy','cohe','pp']
alldata=list()
colors=[[0,0,1],[1,0,0],[0,0.5,0.5],[0.5,1,0],[1,0.5,1]]
lb1='causality between cells i= '+str(int(outdegree)) +' and cell j= '+str(int(indexj))
lb2='causality between cells i= '+str(int(indegree)) +' and cell j= '+str(int(indexi))
p1.plot(F,Fxy,label=lb1,linewidth=1.5,c=colors[0])
p1.plot(F,Fy2x,label=lb2,linewidth=1.5,c=colors[1])
#p1.xlim(0,2)
#p1.ylim(-0.4,0.4)
p1.xlabel("Hz")
p1.xlim(0,100)
p1.ylim(-2,2)
p1.title('GSC components'+bc)
p1.legend()#bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
#      ncol=2, mode="expand", borderaxespad=0.)
p1.hold(False)

F,pp,cohe,Fx2y,Fy2x,Fxy=granger(outout,synapse,20)#,200,100)# nfs2[0], nfs2[1] ,
sfin=str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin) 
fig=p1.figure(fin) 
# Plot Granger spectra
p1.hold(True) #Matlab style hold
labels=list()
labels=['Fx2y','Fy2x','Fxy','cohe','pp']
alldata=list()
colors=[[0,0,1],[1,0,0],[0,0.5,0.5],[0.5,1,0],[1,0.5,1]]
lb1='causality between cells i= '+str(int(outdegree)) +' and cell j= '+str(int(indexj))
lb2='causality between cells i= '+str(int(indegree)) +' and cell j= '+str(int(indexi))
p1.plot(F,Fxy,label=lb1,linewidth=1.5,c=colors[0])
p1.plot(F,Fy2x,label=lb2,linewidth=1.5,c=colors[1])
#p1.xlim(0,2)
#p1.ylim(-0.4,0.4)
p1.xlabel("Hz")
p1.xlim(0,100)
p1.ylim(-2,2)
p1.title('GSC components'+bc)
p1.legend()#bbox_to_anchor=(0., -0.106, 1., .102), loc=8,
#      ncol=2, mode="expand", borderaxespad=0.)
p1.hold(False)

#F,pp,cohe,Fx2y,Fy2x,Fxy=granger(outout,synapse,20)#,200,100)# nfs2[0], nfs2[1] ,
p1.show() 
"""	



"""
def getdirm():	  
 dirsgc=nx.DiGraph()#directed graph. 
 dirsgc=np.array(sgcdir)")
 indexi= 0
 indexj=0
 for i in xrange(0,int(h.ntefa.v[0].size-1)): 

   sources=ntefa.v[0].x[iter]
   targets=ntefa.v[1].x[iter]

  ent=ntefa.v[2].x[iter] //get the nte convert to pythonic form
  //print py.ent
  
  if(ntefa.v[2].x[iter]>thr):
      print py.ent
      cnt+=1
      py.indexi= py.sources
      py.indexj= py.targets
      py.fingure_cnt=py.int(cnt)
  

   entfa[int(sources)][int(targets)] = ent

  if(ent>thr/10):      
      dir2.add_edge(int(sources),int(targets),weight=ent)
  


get_dirm()
"""
def draw_graph(G,
               labels=None,
               node_colors=None,
               node_shapes=None,
               node_scale=1.0,
               edge_style='solid',
               edge_cmap=None,
               colorbar=False,
               vrange=None,
               layout=None,
               title=None,
               font_family='sans-serif',
               font_size=9,
               stretch_factor=1.0,
               edge_alpha=True,
               fig_size=None):
    """
  
    Draw a weighted graph with options to visualize link weights.

    The resulting diagram uses the rank of each node as its size, and the
    weight of each link (after discarding thresholded values, see below) as the
    link opacity.

    It maps edge weight to color as well as line opacity and thickness,
    allowing the color part to be hardcoded over a value range (to permit valid
    cross-figure comparisons for different graphs, so the same color
    corresponds to the same link weight even if each graph has a different
    range of weights).  The nodes sizes are proportional to their degree,
    computed as the sum of the weights of all their links.  The layout defaults
    to circular, but any nx layout function can be passed in, as well as a
    statically precomputed layout.

    Parameters
    ----------
    G : weighted graph
      The values must be of the form (v1,v2), with all v2 in [0,1].  v1 are
      used for colors, v2 for thickness/opacity.

    labels : list or dict, optional.
      An indexable object that maps nodes to strings.  If not given, the
      string form of each node is used as a label.  If False, no labels are
      drawn.

    node_colors : list or dict, optional.
      An indexable object that maps nodes to valid matplotlib color specs.  See
      matplotlib's plot() function for details.

    node_shapes : list or dict, optional.
      An indexable object that maps nodes to valid matplotlib shape specs.  See
      matplotlib's scatter() function for details.  If not given, circles are
      used.

    node_scale : float, optional
      A scale factor to globally stretch or shrink all nodes symbols by.

    edge_style : string, optional
      Line style for the edges, defaults to 'solid'.

    edge_cmap : matplotlib colormap, optional.
      A callable that returns valid color specs, like matplotlib colormaps.
      If not given, edges are colored black.

    colorbar : bool
      If true, automatically add a colorbar showing the mapping of graph weight
      values to colors.

    vrange : pair of floats
      If given, this indicates the total range of values that the weights can
      in principle occupy, and is used to set the lower/upper range of the
      colormap.  This allows you to set the range of multiple different figures
      to the same values, even if each individual graph has range variations,
      so that visual color comparisons across figures are valid.

    layout : function or layout dict, optional
      A NetworkX-like layout function or the result of a precomputed layout for
      the given graph.  NetworkX produces layouts as dicts keyed by nodes and
      with (x,y) pairs of coordinates as values, any function that produces
      this kind of output is acceptable.  Defaults to nx.circular_layout.

    title : string, optional.
      If given, title to put on the main plot.

    font_family : string, optional.
      Font family used for the node labels and title.

    font_size : int, optional.
      Font size used for the node labels and title.

    stretch_factor : float, optional
      A global scaling factor to make the graph larger (or smaller if <1).
      This can be used to separate the nodes if they start overlapping.

    edge_alpha: bool, optional
      Whether to weight the transparency of each edge by a factor equivalent to
      its relative weight

    fig_size: list of height by width, the size of the figure (in
    inches). Defaults to [6,6]

    Returns
    -------
    fig
      The matplotlib figure object with the plot.
    """
    if fig_size is None:
        figsize = [6, 6]

    scaler = figsize[0] / 6.
    # For the size of the node symbols
    node_size_base = 1000 * scaler
    node_min_size = 200 * scaler
    default_node_shape = 'o'
    # Default colors if none given
    default_node_color = 'r'
    default_edge_color = 'k'
    # Max edge width
    max_width = 13 * scaler
    min_width = 2 * scaler
    font_family = 'sans-serif'

    # We'll use the nodes a lot, let's make a numpy array of them
    nodes = np.array(sorted(G.nodes()))
    nnod = len(nodes)

    # Build a 'weighted degree' array obtained by adding the (absolute value)
    # of the weights for all edges pointing to each node:
    amat = nx.adj_matrix(G).A  # get a normal array out of it
    degarr = abs(amat).sum(0)  # weights are sums across rows

    # Map the degree to the 0-1 range so we can use it for sizing the nodes.
    try:
        odegree = rescale_arr(degarr, 0, 1)
        # Make an array of node sizes based on node degree
        node_sizes = odegree * node_size_base + node_min_size
    except ZeroDivisionError:
        # All nodes same size
        node_sizes = np.empty(nnod, float)
        node_sizes.fill(0.5 * node_size_base + node_min_size)

    # Adjust node size list.  We square the scale factor because in mpl, node
    # sizes represent area, not linear size, but it's more intuitive for the
    # user to think of linear factors (the overall figure scale factor is also
    # linear).
    node_sizes *= node_scale ** 2

    # Set default node properties
    if node_colors is None:
        node_colors = [default_node_color] * nnod

    if node_shapes is None:
        node_shapes = [default_node_shape] * nnod

    # Set default edge colormap
    if edge_cmap is None:
        # Make an object with the colormap API, that maps all input values to
        # the default color (with proper alhpa)
        edge_cmap = (lambda val, alpha:
                     colors.colorConverter.to_rgba(default_edge_color, alpha))

    # if vrange is None, we set the color range from the values, else the user
    # can specify it

    # e[2] is edge value: edges_iter returns (i,j,data)
    gvals = np.array([e[2]['weight'] for e in G.edges(data=True)])
    gvmin, gvmax = gvals.min(), gvals.max()

    gvrange = gvmax - gvmin
    if vrange is None:
        vrange = gvmin, gvmax
    # Now, construct the normalization for the colormap
    cnorm = mpl.colors.Normalize(vmin=vrange[0], vmax=vrange[1])

    # Create the actual plot where the graph will be displayed
    figsize = np.array(figsize, float)
    figsize *= stretch_factor

    fig = plt.figure(figsize=figsize)
    ax_graph = fig.add_subplot(1, 1, 1)
    fig.sca(ax_graph)

    if layout is None:
        layout = nx.circular_layout
    # Compute positions for all nodes - nx has several algorithms
    if callable(layout):
        pos = layout(G)
    else:
        # The user can also provide a precomputed layout
        pos = layout

    # Draw nodes
    for nod in nodes:
        nx.draw_networkx_nodes(G,
                               pos,
                               nodelist=[nod],
                               node_color=node_colors[nod],
                               node_shape=node_shapes[nod],
                               node_size=node_sizes[nod])
    # Draw edges
    if not isinstance(G, nx.DiGraph):
        # Undirected graph, simple lines for edges
        # We need the size of the value range to properly scale colors
        vsize = vrange[1] - vrange[0]
        gvals_normalized = G.metadata['vals_norm']
        for (u, v, y) in G.edges(data=True):
            # The graph value is the weight, and the normalized values are in
            # [0,1], used for thickness/transparency
            alpha = gvals_normalized[u, v]
            # Scale the color choice to the specified vrange, so that
            ecol = (y['weight'] - vrange[0]) / vsize
            #print 'u,v:',u,v,'y:',y,'ecol:',ecol  # dbg

            if edge_alpha:
                fade = alpha
            else:
                fade = 1.0

            edge_color = [tuple(edge_cmap(ecol, fade))]
            #dbg:
            #print u,v,y
            draw_networkx_edges(G,
                                pos,
                                edgelist=[(u, v)],
                                width=min_width + alpha * max_width,
                                edge_color=edge_color,
                                style=edge_style)
    else:
        # Directed graph, use arrows.
        # XXX - this is currently broken.
        raise NotImplementedError("arrow drawing currently broken")

        ## for (u,v,x) in G.edges(data=True):
        ##     y,w = x
        ##     draw_arrows(G,pos,edgelist=[(u,v)],
        ##                 edge_color=[w],
        ##                 alpha=w,
        ##                 edge_cmap=edge_cmap,
        ##                 width=w*max_width)

    # Draw labels.  If not given, we use the string form of the nodes.  If
    # labels is False, no labels are drawn.
    if labels is None:
        labels = map(str, nodes)

    if labels:
        lab_idx = range(len(labels))
        labels_dict = dict(zip(lab_idx, labels))
        nx.draw_networkx_labels(G,
                                pos,
                                labels_dict,
                                font_size=font_size,
                                font_family=font_family)

    if title:
        plt.title(title, fontsize=font_size)

    # Turn off x and y axes labels in pylab
    plt.xticks([])
    plt.yticks([])

    # Add a colorbar if requested
    if colorbar:
        divider = make_axes_locatable(ax_graph)
        ax_cb = divider.new_vertical(size="20%", pad=0.2, pack_start=True)
        fig.add_axes(ax_cb)
        cb = mpl.colorbar.ColorbarBase(ax_cb,
                                    cmap=edge_cmap,
                                    norm=cnorm,
                                    #boundaries = np.linspace(min((gvmin,0)),
                                    #                         max((gvmax,0)),
                                    #                         256),
                                    orientation='horizontal',
                                    format='%.2f')

    # Always return the MPL figure object so the user can further manipulate it
    return fig


def lab2node(labels, labels_dict):
    return [labels_dict[ll] for ll in labels]


### Patched version for networx draw_networkx_edges, sent to Aric.
def draw_networkx_edges(G, pos,
                        edgelist=None,
                        width=1.0,
                        edge_color='k',
                        style='solid',
                        alpha=None,
                        edge_cmap=None,
                        edge_vmin=None,
                        edge_vmax=None,
                        ax=None,
                        arrows=True,
                        **kwds):
    """Draw the edges of the graph G

    This draws only the edges of the graph G.

    pos is a dictionary keyed by vertex with a two-tuple
    of x-y positions as the value.
    See networkx.layout for functions that compute node positions.

    edgelist is an optional list of the edges in G to be drawn.
    If provided, only the edges in edgelist will be drawn.

    edgecolor can be a list of matplotlib color letters such as 'k' or
    'b' that lists the color of each edge; the list must be ordered in
    the same way as the edge list. Alternatively, this list can contain
    numbers and those number are mapped to a color scale using the color
    map edge_cmap.  Finally, it can also be a list of (r,g,b) or (r,g,b,a)
    tuples, in which case these will be used directly to color the edges.  If
    the latter mode is used, you should not provide a value for alpha, as it
    would be applied globally to all lines.

    For directed graphs, 'arrows' (actually just thicker stubs) are drawn
    at the head end.  Arrows can be turned off with keyword arrows=False.

    See draw_networkx for the list of other optional parameters.

    """
    try:
        import matplotlib.pylab as pylab
        import matplotlib.cbook as cb
        from matplotlib.colors import colorConverter, Colormap
        from matplotlib.collections import LineCollection
    except ImportError:
        raise ImportError("Matplotlib required for draw()")
    except RuntimeError:
        pass  # unable to open display

    if ax is None:
        ax = pylab.gca()

    if edgelist is None:
        edgelist = G.edges()

    if not edgelist or len(edgelist) == 0:  # no edges!
        return None

    # set edge positions
    edge_pos = np.asarray([(pos[e[0]], pos[e[1]]) for e in edgelist])

    if not cb.iterable(width):
        lw = (width,)
    else:
        lw = width

    if not cb.is_string_like(edge_color) \
           and cb.iterable(edge_color) \
           and len(edge_color) == len(edge_pos):
        if np.alltrue([cb.is_string_like(c)
                         for c in edge_color]):
            # (should check ALL elements)
            # list of color letters such as ['k','r','k',...]
            edge_colors = tuple([colorConverter.to_rgba(c, alpha)
                                 for c in edge_color])
        elif np.alltrue([not cb.is_string_like(c)
                           for c in edge_color]):
            # If color specs are given as (rgb) or (rgba) tuples, we're OK
            if np.alltrue([cb.iterable(c) and len(c) in (3, 4)
                             for c in edge_color]):
                edge_colors = tuple(edge_color)
                alpha = None
            else:
                # numbers (which are going to be mapped with a colormap)
                edge_colors = None
        else:
            e_s = 'edge_color must consist of either color names or numbers'
            raise ValueError(e_s)
    else:
        if len(edge_color) == 1:
            edge_colors = (colorConverter.to_rgba(edge_color, alpha),)
        else:
            e_s = 'edge_color must be a single color or list of exactly'
            e_s += 'm colors where m is the number or edges'
            raise ValueError(e_s)
    edge_collection = LineCollection(edge_pos,
                                     colors=edge_colors,
                                     linewidths=lw,
                                     antialiaseds=(1,),
                                     linestyle=style,
                                     transOffset=ax.transData,
                                     )

    # Note: there was a bug in mpl regarding the handling of alpha values for
    # each line in a LineCollection.  It was fixed in matplotlib in r7184 and
    # r7189 (June 6 2009).  We should then not set the alpha value globally,
    # since the user can instead provide per-edge alphas now.  Only set it
    # globally if provided as a scalar.
    if cb.is_numlike(alpha):
        edge_collection.set_alpha(alpha)

    # need 0.87.7 or greater for edge colormaps
    if edge_colors is None:
        if edge_cmap is not None:
            assert(isinstance(edge_cmap, Colormap))
        edge_collection.set_array(np.asarray(edge_color))
        edge_collection.set_cmap(edge_cmap)
        if edge_vmin is not None or edge_vmax is not None:
            edge_collection.set_clim(edge_vmin, edge_vmax)
        else:
            edge_collection.autoscale()
        pylab.sci(edge_collection)

#    else:
#        sys.stderr.write(\
#            """matplotlib version >= 0.87.7 required for colormapped edges.
#        (version %s detected)."""%matplotlib.__version__)
#        raise UserWarning(\
#            """matplotlib version >= 0.87.7 required for colormapped edges.
#        (version %s detected)."""%matplotlib.__version__)

    arrow_collection = None

    if G.is_directed() and arrows:

        # a directed graph hack
        # draw thick line segments at head end of edge
        # waiting for someone else to implement arrows that will work
        arrow_colors = (colorConverter.to_rgba('k', alpha),)
        a_pos = []
        p = 1.0 - 0.25  # make head segment 25 percent of edge length
        for src, dst in edge_pos:
            x1, y1 = src
            x2, y2 = dst
            dx = x2 - x1  # x offset
            dy = y2 - y1  # y offset
            d = np.sqrt(float(dx ** 2 + dy ** 2))  # length of edge
            if d == 0:  # source and target at same position
                continue
            if dx == 0:  # vertical edge
                xa = x2
                ya = dy * p + y1
            if dy == 0:  # horizontal edge
                ya = y2
                xa = dx * p + x1
            else:
                theta = np.arctan2(dy, dx)
                xa = p * d * np.cos(theta) + x1
                ya = p * d * np.sin(theta) + y1

            a_pos.append(((xa, ya), (x2, y2)))

        arrow_collection = LineCollection(a_pos,
                                colors=arrow_colors,
                                linewidths=[4 * ww for ww in lw],
                                antialiaseds=(1,),
                                transOffset=ax.transData,
                                )

    # update view
    minx = np.amin(np.ravel(edge_pos[:, :, 0]))
    maxx = np.amax(np.ravel(edge_pos[:, :, 0]))
    miny = np.amin(np.ravel(edge_pos[:, :, 1]))
    maxy = np.amax(np.ravel(edge_pos[:, :, 1]))

    w = maxx - minx
    h = maxy - miny
    padx, pady = 0.05 * w, 0.05 * h
    corners = (minx - padx, miny - pady), (maxx + padx, maxy + pady)
    ax.update_datalim(corners)
    ax.autoscale_view()

    edge_collection.set_zorder(1)  # edges go behind nodes
    ax.add_collection(edge_collection)
    if arrow_collection:
        arrow_collection.set_zorder(1)  # edges go behind nodes
        ax.add_collection(arrow_collection)

    return ax


def mkgraph(cmat, threshold=0.0, threshold2=None):
    """Make a weighted graph object out of an adjacency matrix.

    The values in the original matrix cmat can be thresholded out.  If only one
    threshold is given, all values below that are omitted when creating edges.
    If two thresholds are given, then values in the th2-th1 range are
    ommitted.  This allows for the easy creation of weighted graphs with
    positive and negative values where a range of weights around 0 is omitted.

    Parameters
    ----------
    cmat : 2-d square array
      Adjacency matrix.
    threshold : float
      First threshold.
    threshold2 : float
      Second threshold.

    Returns
    -------
    G : a NetworkX weighted graph object, to which a dictionary called
    G.metadata is appended.  This dict contains the original adjacency matrix
    cmat, the two thresholds, and the weights
    """

    # Input sanity check
    nrow, ncol = cmat.shape
    if nrow != ncol:
        raise ValueError("Adjacency matrix must be square")

    row_idx, col_idx, vals = threshold_arr(cmat, threshold, threshold2)
    # Also make the full thresholded array available in the metadata
    cmat_th = np.empty_like(cmat)
    if threshold2 is None:
        cmat_th.fill(threshold)
    else:
        cmat_th.fill(-np.inf)
    cmat_th[row_idx, col_idx] = vals

    # Next, make a normalized copy of the values.  For the 2-threshold case, we
    # use 'folding' normalization
    if threshold2 is None:
        vals_norm = minmax_norm(vals)
    else:
        vals_norm = minmax_norm(vals, 'folding', [threshold, threshold2])

    # Now make the actual graph
    G = nx.Graph(weighted=True)
    G.add_nodes_from(range(nrow))
    # To keep the weights of the graph to simple values, we store the
    # normalize ones in a separate dict that we'll stuff into the graph
    # metadata.

    normed_values = {}
    for i, j, val, nval in zip(row_idx, col_idx, vals, vals_norm):
        if i == j:
            # no self-loops
            continue
        G.add_edge(i, j, weight=val)
        normed_values[i, j] = nval

    # Write a metadata dict into the graph and save the threshold info there
    G.metadata = dict(threshold1=threshold,
                      threshold2=threshold2,
                      cmat_raw=cmat,
                      cmat_th=cmat_th,
                      vals_norm=normed_values,
                      )
    return G


"""
