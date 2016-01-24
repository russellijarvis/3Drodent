#%run #ie this line keeps ipython open for debugging after running this program.

#saving figures is rather impotent god knows why. Its not very rewarding to debug the problem. So I will just save all the traces in run specific files and then open them and plot them all at once.



##
# To avoid name conflict Matrix is titled Matrix_here.
# 
##


#from numpy import *
#import scipy
import numpy as np
import pylab as plt
#from pylab import *
#from pylab import *
#from matplotlib.pyplot import step, show, xlim
import re
#from mpi4py import MPI

#import neuron
#from neuron import h

#If you redeclare the hoc object 'h', it means the old hoc objects contents are deleted. Very bad!

#num_cells=26#int(h.cells.count)

#print num_cells, "this is num_cells"

#num_cells = int(len(num_cells))
#num_cells = int(len(nfs))
#for i in sources
mbin=np.array(mbin)

column_index=0
row_index=0
old_row=0
old_column=0
for i in xrange(0,int(mbin.shape[1])):
 for j in xrange(0,int(mbin.shape[0])):
 
  if sum(mbin[j,:]) > old_row: # old :
   old_row = sum(mbin[j,:])
   row_index=j
   print old_row
  if sum(mbin[:,i]) > old_column: # old :
   old_column = sum(mbin[:,i])
   column_index=i
   print old_column


outdegree=column_index
indegree=row_index

print outdegree
print indegree  



listlargest=[]
print sources[:]#, max(sources)
print targets[:]
listlargest.append(max(sources))
listlargest.append(max(targets))
largest_num_cells=int(max(listlargest)) 


listsmallest=[]
 
listsmallest.append(len(sources))
listsmallest.append(len(targets))
min_num_cells=int(min(listsmallest)) 

#for i in range(0,len(sources)-1): 
# print sources[i], " ", targets[i]
 
 
# requirements. Sort list of sources.
# THere are 98 cells but only 74 of these cells are added to the list. So have to conclude that 24 of these cells are cell fragments.
# artifacts from the reconstruction process.
#num_cells=100

     

#def conn_matrix(binary,plot):

print num_cells

old = 0
binary_plot = [[0 for x in xrange(0,int(num_cells)+1)] for x in xrange(0,int(num_cells)+1)] 
Matrix_here = [[0 for x in xrange(0,int(num_cells)+1)] for x in xrange(0,int(num_cells)+1)] 
for i in xrange(0,len(sources)-1): 
  
  Matrix_here[int(sources[i])][int(targets[i])] =Matrix_here[int(sources[i])][int(targets[i])] + 1 
  binary_plot[int(sources[i])][int(targets[i])] = 1 
  #The degree of connections will be graphed as 
     
  # print  "sources ", indexs_sources, "targets ", indexs_targets
  #print "source cell number is: ", i, " target at this source is simply: ", targets[i], " source is simply: ", sources[i]
  #Matrix_here[0][j]=targets[i] #rows are sources
  #Matrix_here[int(sources[i])][int(targets[i])] = 1 #The degree of connections will be graphed as 
 # if binary==1:
#   Matrix_here[int(sources[i])][int(targets[i])] = 1
 # else:
  #[rows=sources][columns=targets]


os.chdir('results')
Matrix=np.array(Matrix)
Matrix_here=np.array(Matrix_here)#was3
binary_plot=np.array(binary_plot)#Numpy matrix is very flexible. ie a=Matrix2[295,0:400]  
f = open('Matrix_here', 'w')
np.savetxt('Matrix_here',Matrix_here)
f = open('binary_plot', 'w')
np.savetxt('binary_plot',binary_plot)
f = open('Matrix', 'w')
np.savetxt('Matrix',Matrix)
os.chdir('../')

Weights2=np.array(Weights)
Weights=Weights2
 # a colour
   #A = rand(5,5)
for j in xrange(0,int(Matrix_here.shape[0])):
 print sum(Matrix_here[j,:]), j 
for i in xrange(0,int(Matrix_here.shape[1])):
 print sum(Matrix_here[:,i]), i

# print sum(Matrix_here2[:,24]), i

column_index=0
row_index=0
old_row=0
old_column=0
for i in xrange(0,int(Matrix_here.shape[1])):
 for j in xrange(0,int(Matrix_here.shape[0])):
  if sum(Matrix_here[j,:]) > old_row: # old :
   old_row = sum(Matrix_here[j,:])
   row_index=j
  if sum(Matrix_here[:,i]) > old_column: # old :
   old_column = sum(Matrix_here[:,i])
   column_index=i



outdegree=column_index
indegree=row_index     

os.chdir('../')

#Matrix_here and conn_matrix are the same except one is the transpose of the other. Which is correct?
   
    
print column_index, "high in degree cell"
Matrix_here[:,column_index]
print row_index, "high out degree cell"
Matrix_here[row_index,:]

#nrnpython("print sum(Matrix_here2[:,column_index])")
#nrnpython("print sum(Matrix_here2[row_index,:])")

##
# Call this program from NEURON using: 
# nrnpython("execfile('connection_matrix.py')")
# Save these results to a file each.
# scanf these results into NEURON.
##

h.chdir(h.workingdir)

# 
#if plot==1:

#plt2, ax = plt.subplots(2)
from matplotlib.colors import LogNorm


plt.figure(1)

X=np.ones(num_cells)
Y=np.ones(num_cells)
Z=Matrix_here
X1=np.linspace(0,num_cells)
#,len(),len())
im=plt.imshow(Matrix_here, norm=LogNorm())
#plt.pcolor(X, Y, Matrix_here,edgecolor='white',linewidth=5)   
#im=plt.imshow(Matrix_here, norm=LogNorm())
#plt.colorbar(im)
plt.xlabel('columns = targets')
plt.ylabel('rows = sources')
plt.title('Adjacency matrix')

#plt.savefig("connmat.svg")
#plt.draw()
#plt.close()
#plt.xticks(X, [1,2,3,4,5,6,7])
#plt.autoscale(True)
#xticks = plt.get_xticks()
#plt.xaxis.set_ticklabels(scale_xaxis(xticks))
#X=np.arrange(len(num_cells))
#Y=np.arrange(len(num_cells))
#plt2.show()
#plt.savefig('connection_matrix3.png')  
#im=plt.imshow(Weights2)

plt.figure(2)

#Matrix_here[:,column_index]=100
#Matrix_here[row_index,:]=100
im=plt.imshow(binary_plot)

plt.colorbar(im)
plt.xlabel('columns = targets')
plt.ylabel('rows = sources')
plt.title('Weight matrix')
#plt.savefig('connection_matrix3.png')  

#plt.savefig("weightmat.svg")
#plt.draw()
#plt.close()
plt.figure(3)

X=np.ones(num_cells)
Y=np.ones(num_cells)
Z=Matrix_here
#,len(),len())
im=plt.imshow(Matrix_here, norm=LogNorm())
#plt.pcolor(X, Y, Matrix_here,edgecolor='white',linewidth=2)   
#im=plt.imshow(Matrix_here, norm=LogNorm())
plt.colorbar(im)
plt.xlabel('columns = targets')
plt.ylabel('rows = sources')
plt.title('Adjacency matrix')

plt.autoscale(True)

#plt.savefig("connmatverification.svg")
#plt.draw()
#plt.close()
#xticks = plt.get_xticks()
#plt.xaxis.set_ticklabels(scale_xaxis(xticks))
#X=np.arrange(len(num_cells))
#Y=np.arrange(len(num_cells))
#plt2.show()
#plt.savefig('connection_matrix3.png')  
#im=plt.imshow(Weights2)




#plt.show()

 #return (row_index,column_index)

#result1=conn_matrix(0,1) #binary, plot
#result2=conn_matrix(1,1) #This time binary is chosen.


