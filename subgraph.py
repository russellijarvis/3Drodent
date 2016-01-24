
   
   
   
   
   
def cdc(graph):
 g=graph()
 dc=nx.degree_centrality(g)
 nx.set_node_attributes(g,'degree_cent',dc)
 dcsorted=sorted(dc.items(),key=itemgetter(1),reverse=True)
 for key,value in dcsorted[0:10]:
 #for key,value in degcent_sorted[0:10]:
  print "highest degree", key, value
 return graph, dcsorted[0:10]

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
 cent_items.reverse()
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

