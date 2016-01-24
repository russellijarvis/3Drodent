
from operator import itemgetter

def plotx(args,args2):
 #pos=nx.get_node_attributes(args,'pos')
 pos=nx.spring_layout(args)
 # find node near center (0.5,0.5)
 dmin=2
 ncenter=0
 for n in pos:
     x,y=pos[n]
     d=(x-0.5)**2+(y-0.5)**2
     if d<dmin:
         ncenter=n
         dmin=d

 # color by path length from node near center
 
 
 """
 node_and_degree=msgc3.degree()
 (largest_hub,degree)=sorted(node_and_degree.items(),key=itemgetter(1))[-1]
 # Create ego graph of main hub
 hub_ego=nx.ego_graph(msgc3,largest_hub)
 pos=nx.spring_layout(hub_ego)

 nx.draw(hub_ego,pos,node_color='b',node_size=50,with_labels=False)
 # Draw ego as large and red
 nx.draw_networkx_nodes(hub_ego,pos,node_color='r')
 plt.savefig('ego_graph.png')
 p=nx.single_source_shortest_path_length(args,ncenter)
 """ 
 

 plt.figure(1)
 plt.clf()
 nx.draw_networkx(args)
 #nx.draw_networkx_edges(args,pos)#,nodelist=[ncenter],alpha=0.54)
 #nx.draw_networkx_nodes(args,pos)#,nodelist=p.keys(),node_color=p.values(),cmap=plt.cm.Reds_r)
 #nx.draw_networkx_labels(args,pos)
 

 s=args2+'randomgraph.png'
 plt.savefig(s)
 return

plotx(dirfinal,'dirfinal')
