from allensdk.model.biophys_sim.config import Config
from utils import Utils
#from allensdk.model.biophysical_perisomatic.utils import Utils

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import glob

config = Config().load('config.json')
#config = Config().load('473863510_fit.json')

# This configuration can be modified inside Python!
# Delete the contents here and make up my own.

# configure NEURON
utils = Utils(config)
h = utils.h

# configure model
manifest = config.manifest
utils.generate_cells()
utils.connect_cells()
utils.connect_ring()

#utils.generate_morphology()
def read_local_swc(utils):
    morphs=[]
    cells1=[]
    Utils.initialize_hoc(utils)
    swclist=glob.glob('*.swc')
    for swcf in swclist:
        #if swcf != "466664172.swc":
        #morphology = swc.read_swc(swcf)
        print swcf
        #morphology=
        print utils.generate_morphology(utils,swcf)
        
        #cells1.append(h.Cell[0])
        #print morphology
        #return morphs
        #cell1=Utils.load_cell_parameters(utils)       
        #cells1.append(cell1)
            #print type(cells1)
            #print type(cell1)
            #morphology.root
            #morphs.append(morphology)
            #eturn morphs,swclist#,cells1
            #morphs,swclist=
#morphs=
#utils.generate_cells()
#sread_local_swc(utils)

h('topology()')
h('forall psection()')
#
#generate_morphology(self, cell, morph_filename) unbound utils.Utils method


#utils.connect_cells()

# configure stimulus
utils.setup_iclamp_step(utils.cells[0], 0.27, 1020.0, 750.0)
h.dt = 0.025
h.tstop = 3000

# configure recording
vec = utils.record_values()

# run the model
h.finitialize()
h.run()

# save output voltage to text file
data = np.transpose(np.vstack((vec["t"],
                               vec["v"][0],
                               vec["v"][1],
                               vec["v"][2])))
np.savetxt('multicell.dat', data)

# use matplotlib to plot to png image
fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True)
for i in range(len(utils.cells)):
    axes[i].plot(vec["t"], vec["v"][i])
    axes[i].set_title(utils.cells_data[i]["type"])
plt.tight_layout()
plt.savefig('multicell.png')


#Never forget dir(h) works
#As does dir(nrn)

h('topology()')
h('forall psection()')
#nrn('topology()')
import mpi4py as MPI
import numpy as np
import csv
import pickle
import sys
import glob

#COMM = MPI.COMM_WORLD
#SIZE = COMM.Get_size()
#RANK = COMM.Get_rank()
pc = h.ParallelContext()    

srcs=[]
tgts=[]
pc.barrier()
itlist=range(98673-1)    
h('objref py')
h('py = new PythonObject()')
h('nrnpython("gidpop={}")') 
h('NCELL=5500000')



#This checks that all cells have synapses and they are not just APcounts.                                                                                                      
hocstring='''
objref strobj \n
strobj = new StringFunctions() \n
objref py \n
py = new PythonObject() \n
nrnpython("strlist=[]") \n
nrnpython("gidlist=[]") \n
strdef cell \n
for i=0, 50000000-1{ \n
    if (pc.gid_exists(i)) { \n
       sprint(cell,"%s",pc.gid2cell(i)) \n
       py.strlist.append(cell) \n
       py.gidlist.append(i) \n
    } \n
} \n
nrnpython("gidstr=zip(strlist,gidlist)") \n
'''
h(hocstring)
gidstr='gidstr'

gidpop={}
pickle.dump(gidstr,open("gidpop"+str(pc.id())+".p", "wb" ))
pc.barrier()

#make it such that each host gets a population dictionary of that contains entries pertaining to all hosts.
#gidpopls=glob.glob('gidpop*.*.p')
#srcs=[]
#tgts=[]
#gidpop=[]
#for s in gidpopls:
#    print s
#    new_dict = pickle.load( open( s, "rb" ) )
#    gidpop = gidpop + new_dict
#    print len(gidpop), len(new_dict)
#
#
#
##Manually convert the gidpop list to a dictionary, since automatic conversion has unintended results.    
#dic={}
#for l in gidpop:
#    dic[l[1]]=l[0]
#


pc.barrier()
lsoftup=[]
ncsize=len(h.NetCon)
#make a list of tuples where each list element contains (srcid,tgtid,srcpop,tgtpop)
for i in xrange(0,ncsize-1):
    srcs.append(int(h.NetCon[i].srcgid()))
    tgts.append(int(h.NetCon[i].postcell().gid))
    srcind=int(h.NetCon[i].srcgid())
    tgtind=int(h.NetCon[i].postcell().gid)
    #print strlist[tgtind]==dic[tgtind], ' sanity check '
    #add to list of tuples, netcon src, netcon tgt, src index, target index.
    lsoftup.append((int(h.NetCon[i].srcgid()),int(h.NetCon[i].postcell().gid),dic[srcind],dic[tgtind]))
    
#pickle dump the list of tuples.
pickle.dump(lsoftup, open( 'lsoftup'+str(pc.id())+'.p', 'wb' ))
