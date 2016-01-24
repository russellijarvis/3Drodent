
#Now that I have limited the number of connections these colours for the degree matrices are realistic.

import pyhoc    
#import statsmodels as statsmodels
from matplotlib.colors import LogNorm 
import matplotlib 
#matplotlib.use('Agg') 
import pylab as p2
import pylab as p1
#from scipy import signal 
from bsmart import granger # Load the Granger calculation tool
#and user make sure that param6.hoc has run.
# sfin=str(h.prunenet)+str(fin)+'.png' 
#fig.savefig(sfin) 

import matplotlib
matplotlib.use('Agg') 


#nrnpython("medr=np.arange(0,int(numcell),0)")
#nrnpython("np.divide(input_marke,rates,medr)")
#nrnpython("meanedr=np.mean(medr)")
#nrnpython("print meanedr, 'mean entropy divided by firing ratem, this should be low for FF, high for FB')


fin=0
h('print nclist.count, "nclist.count", minr, "minr"')

fin+=1
fig=p1.figure(fin)
fig.clf()
if int(h.get_dist)==1:
  #h('py.the_length=py.int(h.adjnqs.v[8].size)'); across=np.arange(0,the_length,0)
  h('py.dist=adjnqs.v[8].to_python()')
  h('py.transfer=adjnqs.v[7].to_python()')
  #p1.plot(across,transfer)
  p1.scatter(transfer,dist)
  p1.title('arc length versus transfer impedence at synapse positions 100Hz')
  p1.xlabel("um")
  p1.ylabel("Mohm")
  sfin='arc length '+str(int(h.prunenet))+str(tstop)+str (h.plastic)+str(int(h.ff))+str(numcell)+str(int(fin))+'.png'
  p1.savefig(sfin, format='png') 


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




