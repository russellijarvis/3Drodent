#!/usr/bin/python
# -*- coding: utf-8 -*-

##
# A adjacency list is deliberately chosen over a degree list, such that the degree of convergent connections does not
# bias the result, and only cells with the most number of different inputs/outputs are considered to be high degree.
# This file merely finds the vertrices with the highest number of degree 1 neighbours, for both indegree and outdegree.
# A more comprehensive test, would test for number degree n>=2 neighbours also.
##
# import py_compile
# py_compile.compile("mymodule.py")

import numpy as np
import networkx as nx
import pylab as p1
from matplotlib import patches
from matplotlib.colors import LogNorm

import community as louvain# works on trifid.
from collections import defaultdict # works on trifid.
#import brainx
# import re
'''
# It would be smarter to write this as a function that deals with a generic matrix, not a specific matrix.
# These methods are broken, I am not sure why.
# ma=np.array(conn_matrix)

mg =np.array(conm) #np.array(MatrixG)
ma =np.array(conm) #np.array(MatrixA)
structural = np.array(conm)

type(mg)
type(ma)
type(structural)
# ri2=0
# ci2=0
'''
def get_in(ma):
    old_row = 0
    row_index = 0
    old_j = 0
    for j in xrange(0, int(ma.shape[0])):
        if sum(ma[j, :]) > old_row:  # old :
            old_row = sum(ma[j, :])
            row_index = j
            #print row_index, 'j= ', j, old_row, old_j
            old_j = j
    return row_index


def get_out(ma):
    old_column = 0
    column_index = 0
    old_i = 0
    for i in xrange(0, int(ma.shape[1])):
        if sum(ma[:, i]) > old_column:  # old :
            old_column = sum(ma[:, i])
            column_index = i
            #print column_index, 'i= ', i, old_column, old_i
            old_i = i
            #print ma
    return column_index


# scndo_deg=old_column
# scndi_deg=old_row

outdegree=0
indegree=0
outdegreeg=0
indegreeg=0

# outdegree=column_index
def obtain_values(q,r):
   global outdegree, indegree, indegreeg, outdegreeg
   ma=q
   mg=r
   #structural=s
   """
   if int(h.name_declared("lpy")) != 5:
       outdegree = get_out(ma)
       indegree = get_in(ma)
       outdegreeg = get_out(mg)
       indegreeg = get_in(mg)
       return
   """
   if int(h.ff) == 1:
       outdegree = get_out(structural)
       indegree = get_out(structural)
       outdegreeg = get_out(structural)
       indegreeg = get_in(structural)
       mg = structural
       ma = structural
   if int(h.ff) != 1:
       outdegree = get_out(ma)
       indegree = get_in(ma)
       outdegreeg = get_out(mg)
       indegreeg = get_in(mg)

   print outdegreeg
   print indegreeg

   print outdegree
   print indegree
   return

#obtain_values(ma,mg)
structural_plotting=1

def draw_adjacency_matrix(G, node_order=None, partitions=[], colors=[]):
    """
    - G is a netorkx graph
    - node_order (optional) is a list of nodes, where each node in G
          appears exactly once
    - partitions is a list of node lists, where each node in G appears
          in exactly one node list
    - colors is a list of strings indicating what color each
          partition should be
    If partitions is specified, the same number of colors needs to be
    specified.
    """
    adjacency_matrix = nx.to_numpy_matrix(G)#, dtype=np.bool, nodelist=node_order)

    #Plot adjacency matrix in toned-down black and white
    fig = p1.figure(figsize=(5, 5)) # in inches
    
    p1.imshow(adjacency_matrix, interpolation='nearest')
    h('py.titlestr=nclist.count/cells.count')

    #titlestr='Louavain Partitions of Adj Matrix'+str(int(h.nclist.count))
    
    #titlestr=titlestr+str(int(h.cells.count))
    p1.title(titlestr)
    # The rest is just if you have sorted nodes by a partition and want to
    # highlight the module boundaries
    assert len(partitions) == len(colors)
    ax = p1.gca()
    for partition, color in zip(partitions, colors):
        current_idx = 0
        for module in partition:
            ax.add_patch(patches.Rectangle((current_idx, current_idx),
                                          len(module), # Width
                                          len(module), # Height
                                          facecolor="none",
                                          edgecolor=color,
                                          linewidth="1"))
            current_idx += len(module)
    fin=0
    sfin='community_partitions'+ str(int(h.prunenet)) + str(tstop) \
    + str(h.plastic) + str(int(h.ff)) + str(numcell) + str(int(fin)) \
    + str(int(h.minr)) + '.png'
    p1.savefig(sfin)   

#if structural_plotting==1:
Ga = nx.to_networkx_graph(ma, create_using=nx.DiGraph())  # directed graph.
Gg = nx.to_networkx_graph(mg, create_using=nx.DiGraph())  # from_numpy_matrix(mg)
#
structuralg = nx.to_networkx_graph(structural,
							  create_using=nx.DiGraph())
#





# Run louvain community finding algorithm
louvain_community_dict = louvain.best_partition(structuralg.to_undirected())
# Convert community assignmet dict into list of communities
louvain_comms = defaultdict(list)
for node_index, comm_id in louvain_community_dict.iteritems():
    louvain_comms[comm_id].append(node_index)
louvain_comms = louvain_comms.values()
nodes_louvain_ordered = [node for comm in louvain_comms for node in comm]

"""
Method is drawn here.
"""
#draw_adjacency_matrix(structuralg, nodes_louvain_ordered, [louvain_comms], ["red"])   


# import matplotlib
# matplotlib.use('Agg')

# nrnpython("medr=np.arange(0,int(numcell),0)")
# nrnpython("np.divide(input_marke,rates,medr)")
# nrnpython("meanedr=np.mean(medr)")
# nrnpython("print meanedr, 'mean entropy divided by firing ratem, this should be low for FF, high for FB')

"""Plot of Degree Histograms to check for scale free structural properties"""
'''
fin = 0
fig=p1.figure(fin)
p1.clf()

nx.draw(structuralg, node_size=5, node_color='r', edge_color='b',
   alpha=.2)

sfin = 'Structural Graph of the Network' + str(int(h.prunenet)) \
+ str(tstop) + str(h.plastic) + str(int(h.ff)) + str(numcell) \
+ str(int(fin)) + '.png'
p1.savefig(sfin)

out_degrees = structuralg.out_degree()

# dictionary node:degree

out_values = sorted(set(out_degrees.values()))
out_hist = [out_degrees.values().count(x) for x in out_values]

in_degrees = structuralg.in_degree()

# dictionary node:degree

in_values = sorted(set(in_degrees.values()))
in_hist = [in_degrees.values().count(x) for x in in_values]
fin += 1
p1.figure()
p1.clf()
p1.hold(True)
p1.plot(in_values, in_hist, 'ro-')

# in-degree

p1.plot(out_values, out_hist, 'bv-')
p1.xscale('log')
p1.yscale('log')
p1.hold(False)

# out-degree

p1.legend(['In-degree', 'Out-degree'])
p1.xlabel('Degree Scale free network')
p1.ylabel('Number of Neurons that have this degree value')
p1.title('Degree Scale free network')
sfin = 'dirg_degree_distribution.png' + str(int(h.prunenet)) \
+ str(tstop) + str(h.plastic) + str(int(h.ff)) + str(numcell) \
+ str(int(fin)) + '.png'
p1.savefig(sfin)

fin = 0
#h('print nclist.count, "nclist.count", minr, "minr"')

fin += 1
fig = p1.figure()
fig.clf()

"""Plot of Transfer Impedences at real synapse locations"""


if int(h.get_dist) == 1:

		# h('py.the_length=py.int(h.adjnqs.v[8].size)'); across=np.arange(0,the_length,0)

		h('py.dist=adjnqs.v[7].to_python()')
		h('py.transfer=adjnqs.v[8].to_python()')

		# p1.plot(across,transfer)
		# p1.scatter(transfer,dist)

		xdata = transfer
		ydata = dist

		p1.scatter(xdata, ydata)
		p1.hold(True)
		(slope, yint) = p1.polyfit(xdata, ydata, 1)
		xline = p1.xticks()[0]
		yline = map(lambda x: slope * x + yint, xline)
		p1.plot(xline, yline, ls='--', color='b')

		# Set new x- and y-axis limits

		p1.xlim((0.0, max(xdata) + .15 * max(xdata)))
		p1.ylim((0.0, max(ydata) + .15 * max(ydata)))

		p1.title('arc length versus transfer impedence at synapse positions 100Hz'
				)
		p1.xlabel('um')
		p1.ylabel('Mohm')
		sfin = 'arc length ' + str(int(h.prunenet)) + str(tstop) \
		   + str(h.plastic) + str(int(h.ff)) + str(numcell) \
		   + str(int(fin)) + '.png'
		p1.savefig(sfin)

		fin += 1
		fig = p1.figure()
		fig.clf()


		"""
		The following is a scatter plot like above except only for a single cell.

		"""
		if int(h.get_dist) == 1:

				h('py.targets=adjnqs.v[0].to_python()')
				targets = np.array(targets)
				transfer = np.array(transfer)
				dist = np.array(dist)
				np.where(targets == indegree)

				# h('py.the_length=py.int(h.adjnqs.v[8].size)'); across=np.arange(0,the_length,0)

				# p1.plot(across,transfer)

				xdata = transfer[np.where(targets == indegree)]
				ydata = dist[np.where(targets == indegree)]

				# p1.scatter(transfer[np.where(targets==indegree)],dist[np.where(targets==indegree)])

				p1.scatter(xdata, ydata)
				p1.hold(True)
				(slope, yint) = p1.polyfit(xdata, ydata, 1)
				xline = p1.xticks()[0]
				yline = map(lambda x: slope * x + yint, xline)
				p1.plot(xline, yline, ls='--', color='b')

				# Set new x- and y-axis limits

				p1.xlim((0.0, max(xdata) + .15 * max(xdata)))
				p1.ylim((0.0, max(ydata) + .15 * max(ydata)))
				p1.title('arc length versus transfer impedence at synapse positions 100Hz'
						)
				p1.xlabel('um')
				p1.ylabel('Mohm')
				sfin = 'arc length versus distance but only on indegree cell' \
				   + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
				   + str(int(h.ff)) + str(numcell) + str(int(fin)) + '.png'
				p1.savefig(sfin)

				# h.

				fin = fin + 1
				fig = p1.figure()
				fig.clf()
				im = p1.imshow(mbin, interpolation='nearest')
				p1.colorbar(im)
				p1.xlabel('columns = targets')
				p1.ylabel('rows = sources')
				p1.title('Adjacency matrix both transmitters')
				p1.autoscale(True)
				p1.grid(True)

				# p1.figure(5)

				fig = p1.figure()
				fin = fin + 1
				sfin = 'Both_transmitters_adjacency_matrix' + str(int(h.prunenet)) \
				+ str(tstop) + str(h.plastic) + str(int(h.ff)) + str(numcell) \
				+ str(int(fin)) + '.png'
				fig.savefig(sfin)

				fin = fin + 1
				fig = p1.figure()
				fig.clf()
				im = p1.imshow(Matrix, interpolation='nearest')
				p1.colorbar(im)
				p1.xlabel('columns = targets')
				p1.ylabel('rows = sources')
				p1.title('Degree matrix Excitatory and Inhibitory Connections')
				p1.autoscale(True)
				p1.grid(True)
				sfin = 'Both_transmitters_degree_matrix' + str(int(h.prunenet)) \
				+ str(tstop) + str(h.plastic) + str(int(h.ff)) + str(numcell) \
				+ str(int(fin)) + '.png'
				fig.savefig(sfin)

				# p1 .grid(b='on',which='minor')   ##, which='major', axis='both', **kwargs)....
				# p1.ax.grid(color='r',linewidth=1)

				# im=p1.imshow(Matrix)

				fin = fin + 1
				fig = p1.figure()
				fig.clf()

				# im=p1.imshow(MatrixG)

				im = p1.imshow(MatrixG, interpolation='nearest')
				p1.autoscale(True)
				p1.colorbar(im)
				p1.xlabel('columns = targets')
				p1.ylabel('rows = sources')
				p1.title('Degree matrix GABA')
				p1.grid(True)
				sfin = 'GABA_degree_matrix' + str(int(h.prunenet)) + str(tstop) \
				+ str(h.plastic) + str(int(h.ff)) + str(numcell) + str(int(fin)) \
				+ '.png'
				fig.savefig(sfin)

				fin = fin + 1
				fig = p1.figure()
				fig.clf()

				# im=p1.imshow(MatrixA)

				im = p1.imshow(MatrixA, interpolation='nearest')
				p1.colorbar(im)
				p1.xlabel('columns = targets')
				p1.ylabel('rows = sources')
				p1.title('Degree matrix AMPA')
				p1.grid(True)
				p1.autoscale(True)
				sfin = 'AMPA_degree_matrix' + str(int(h.prunenet)) + str(tstop) \
				+ str(h.plastic) + str(int(h.ff)) + str(numcell) + str(int(fin)) \
				+ '.png'
				fig.savefig(sfin)

				# j=len(colors)-1

				numcell = numcell
				numcell = numcell
				roi_names = [0 for x in xrange(0, int(numcell))]
				for i in xrange(0, int(numcell - 1)):
						roi_names[i] = str(i) + ' ' + str(h.cells.o(i).reponame)
				#print roi_names[i]

				# In the future it would be nice to label every xtick with
				# the elements of roi_names
"""
inimp0 = []
inimp20 = []
inimp100 = []
h('for i=0,cells.count-1{ py.inimp0.append(cells.o(i).inimp0) }')
h('for i=0,cells.count-1{ py.inimp20.append(cells.o(i).inimp20) }')
h('for i=0,cells.count-1{ py.inimp100.append(cells.o(i).inimp100) }')

across = np.arange(0, int(numcell), 1)

# p1.plot(tvec,intervec,'bs')#p1.plot(tvec,pyra,'g^')

fin += 1

p1.figure()
p1.hold(True)

p1.plot(across, np.array(inimp0), 'ro', linewidth=1,
   label='input impedence at DC')
p1.plot(across, np.array(inimp20), 'bo', linewidth=1,
   label='input impedence at 20Hz')  # ,c=colors[j], markeredgecolor = 'none')
j -= 1
p1.plot(across, np.array(inimp100), 'go', linewidth=1,
   label='input impedence at 100Hz')
p1.title('Input Impedence Versus cell')
p1.xlabel('cell number')
p1.ylabel('MOhm')
p1.hold(False)

# p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8, ncol=2, mode="expand", borderaxespad=0.), label='input impedence at 20Hz')

sfin = 'Input Impedence Versus cell' + str(int(h.prunenet)) \
+ str(int(h.ff)) + str(int(fin))
p1.savefig(sfin)
"""
# nx.draw(Ga)

																																																																																																																																																																																																																																																																																																																																																																																																																											 # nx.cycle_graph(55,create_using=Ga)
# list(nx.simple_cycles(Ga)) These lines lead to memory errors.
# list(nx.simple_cycles(Gg))
# print b
# print a
# whos
# p1.show()
# print Ga.out_degree()
# dictionary node:degree
""" Duplicate code

out_values = sorted(set(out_degrees.values()))
out_hist = [out_degrees.values().count(x) for x in out_values]

in_degrees = structuralg.in_degree()

# dictionary node:degree

in_values = sorted(set(in_degrees.values()))
in_hist = [in_degrees.values().count(x) for x in in_values]
fin += 1
p1.figure()
p1.clf()
p1.hold(True)
p1.plot(in_values, in_hist, 'ro-')

# in-degree

p1.plot(out_values, out_hist, 'bv-')
p1.xscale('log')
p1.yscale('log')
p1.hold(False)

# out-degree

p1.legend(['In-degree', 'Out-degree'])
p1.xlabel('Degree Scale free network')
p1.ylabel('Number of Neurons that have this degree value')
p1.title('Degree Scale free network')
sfin = 'dirg_degree_distribution.png' + str(int(h.prunenet)) \
+ str(tstop) + str(h.plastic) + str(int(h.ff)) + str(numcell) \
+ str(int(fin)) + '.png'
p1.savefig(sfin,dpi=fig.dpi)

fin = 0
h('print nclist.count, "nclist.count", minr, "minr"')

fin += 1
fig = p1.figure()
fig.clf()
if int(h.get_dist) == 1:

# h('py.the_length=py.int(h.adjnqs.v[8].size)'); across=np.arange(0,the_length,0)

h('py.dist=adjnqs.v[7].to_python()')
h('py.transfer=adjnqs.v[8].to_python()')

# p1.plot(across,transfer)
# p1.scatter(transfer,dist)

xdata = transfer
ydata = dist

p1.scatter(xdata, ydata)
p1.hold(True)
(slope, yint) = p1.polyfit(xdata, ydata, 1)
xline = p1.xticks()[0]
yline = map(lambda x: slope * x + yint, xline)
p1.plot(xline, yline, ls='--', color='b')

# Set new x- and y-axis limits

p1.xlim((0.0, max(xdata) + .15 * max(xdata)))
p1.ylim((0.0, max(ydata) + .15 * max(ydata)))

p1.title('arc length versus transfer impedence at synapse positions 100Hz'
		)
p1.xlabel('um')
p1.ylabel('Mohm')
sfin = 'arc length ' + str(int(h.prunenet)) + str(tstop) \
   + str(h.plastic) + str(int(h.ff)) + str(numcell) \
   + str(int(fin)) + '.png'
p1.savefig(sfin,dpi=fig.dpi, format='png')

fin += 1
fig = p1.figure()
fig.clf()
if int(h.get_dist) == 1:

h('py.targets=adjnqs.v[0].to_python()')
targets = np.array(targets)
transfer = np.array(transfer)
dist = np.array(dist)
np.where(targets == indegree)

# h('py.the_length=py.int(h.adjnqs.v[8].size)'); across=np.arange(0,the_length,0)

# p1.plot(across,transfer)

xdata = transfer[np.where(targets == indegree)]
ydata = dist[np.where(targets == indegree)]

# p1.scatter(transfer[np.where(targets==indegree)],dist[np.where(targets==indegree)])

p1.scatter(xdata, ydata)
p1.hold(True)
(slope, yint) = p1.polyfit(xdata, ydata, 1)
xline = p1.xticks()[0]
yline = map(lambda x: slope * x + yint, xline)
p1.plot(xline, yline, ls='--', color='b')

# Set new x- and y-axis limits

p1.xlim((0.0, max(xdata) + .15 * max(xdata)))
p1.ylim((0.0, max(ydata) + .15 * max(ydata)))
p1.title('arc length versus transfer impedence at synapse positions 100Hz'
		)
p1.xlabel('um')
p1.ylabel('Mohm')
sfin = 'arc length versus distance but only on indegree cell' \
   + str(int(h.prunenet)) + str(tstop) + str(h.plastic) \
   + str(int(h.ff)) + str(numcell) + str(int(fin)) + '.png'
p1.savefig(sfin,dpi=fig.dpi, format='png')

# h.

fin = fin + 1
fig = p1.figure()
fig.clf()
im = p1.imshow(mbin, interpolation='nearest')
p1.colorbar(im)
p1.xlabel('columns = targets')
p1.ylabel('rows = sources')
p1.title('Adjacency matrix both transmitters')
p1.autoscale(True)
p1.grid(True)

# p1.figure(5)

fig = p1.figure()
fin = fin + 1
sfin = 'Both_transmitters_adjacency_matrix' + str(int(h.prunenet)) \
+ str(tstop) + str(h.plastic) + str(int(h.ff)) + str(numcell) \
+ str(int(fin)) + '.png'
fig.savefig(sfin)

fin = fin + 1
fig = p1.figure()
fig.clf()
im = p1.imshow(Matrix, interpolation='nearest')
p1.colorbar(im)
p1.xlabel('columns = targets')
p1.ylabel('rows = sources')
p1.title('Degree matrix Excitatory and Inhibitory Connections')
p1.autoscale(True)
p1.grid(True)
sfin = 'Both_transmitters_degree_matrix' + str(int(h.prunenet)) \
+ str(tstop) + str(h.plastic) + str(int(h.ff)) + str(numcell) \
+ str(int(fin)) + '.png'
fig.savefig(sfin)

# p1 .grid(b='on',which='minor')   ##, which='major', axis='both', **kwargs)....
# p1.ax.grid(color='r',linewidth=1)

# im=p1.imshow(Matrix)

fin = fin + 1
fig = p1.figure()
fig.clf()

# im=p1.imshow(MatrixG)

im = p1.imshow(MatrixG, interpolation='nearest')
p1.autoscale(True)
p1.colorbar(im)
p1.xlabel('columns = targets')
p1.ylabel('rows = sources')
p1.title('Degree matrix GABA')
p1.grid(True)
sfin = 'GABA_degree_matrix' + str(int(h.prunenet)) + str(tstop) \
+ str(h.plastic) + str(int(h.ff)) + str(numcell) + str(int(fin)) \
+ '.png'
fig.savefig(sfin)

fin = fin + 1
fig = p1.figure()
fig.clf()

# im=p1.imshow(MatrixA)

im = p1.imshow(MatrixA, interpolation='nearest')
p1.colorbar(im)
p1.xlabel('columns = targets')
p1.ylabel('rows = sources')
p1.title('Degree matrix AMPA')
p1.grid(True)
p1.autoscale(True)
sfin = 'AMPA_degree_matrix' + str(int(h.prunenet)) + str(tstop) \
+ str(h.plastic) + str(int(h.ff)) + str(numcell) + str(int(fin)) \
+ '.png'
fig.savefig(sfin)

# j=len(colors)-1

numcell = numcell
numcell = numcell
roi_names = [0 for x in xrange(0, int(numcell))]
for i in xrange(0, int(numcell - 1)):
roi_names[i] = str(i) + ' ' + str(h.cells.o(i).nametype)
print roi_names[i]

# In the future it would be nice to label every xtick with
# the elements of roi_names

inimp0 = []
inimp20 = []
inimp100 = []
h('for i=0,cells.count-1{ py.inimp0.append(cells.o(i).inimp0) }')
h('for i=0,cells.count-1{ py.inimp20.append(cells.o(i).inimp20) }')
h('for i=0,cells.count-1{ py.inimp100.append(cells.o(i).inimp100) }')

across = np.arange(0, int(numcell), 1)

# p1.plot(tvec,intervec,'bs')#p1.plot(tvec,pyra,'g^')

fin += 1

p1.figure()
p1.hold(True)

p1.plot(across, np.array(inimp0), 'ro', linewidth=1,
   label='input impedence at DC')
p1.plot(across, np.array(inimp20), 'bo', linewidth=1,
   label='input impedence at 20Hz')  # ,c=colors[j], markeredgecolor = 'none')
j -= 1
p1.plot(across, np.array(inimp100), 'go', linewidth=1,
   label='input impedence at 100Hz')
p1.title('Input Impedence Versus cell')
p1.xlabel('cell number')
p1.ylabel('MOhm')
p1.hold(False)

# p1.legend(bbox_to_anchor=(0., -0.106, 1., .102), loc=8, ncol=2, mode="expand", borderaxespad=0.), label='input impedence at 20Hz')

sfin = 'Input Impedence Versus cell' + str(int(h.prunenet)) \
+ str(int(h.ff)) + str(int(fin))
p1.savefig(sfin)


# list(nx.simple_cycles(Ga)) These lines lead to memory errors.
# list(nx.simple_cycles(Gg))

## In the future I should free memory by destroying objects
#Ga = nx.to_networkx_graph()  # directed graph.
#Gg = nx.to_networkx_graph()  # from_numpy_matrix(mg)
#structuralg = nx.to_networkx_graph()
"""
'''
