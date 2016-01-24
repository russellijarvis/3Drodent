#!/usr/bin/env python
'''
LFPs from a population of cells relying on MPI
'''
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
import os
import urllib2
import zipfile
#import mpi4py
from mpi4py import MPI
#GPIO pins on microcontoller are both TX and RX
#import neuron
from neuron import h
import LFPy#Importing LFPy imports neuron. NEURON has to be imported after mpi4py


#initialize the MPI interface
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()



pc = h.ParallelContext()
h('objref pc')
h.pc=pc
#h('pc = new ParallelContext()')
#h('pc = py.pc')
s = "mpi4py thinks I am %d of %d,\
 NEURON thinks I am %d of %d\n"
cw = MPI.COMM_WORLD
print s % (cw.rank, cw.size, \
           pc.id(),pc.nhost())
#pc.done()

h.dt = 0.025
tstop = 1025
h('dc=-1')  # delete cell. This retains the cell, but makes all of its section diameter 0us, and it makes all projection weights to it 0.
h('default_channels=0')
h('plastic=1')
h('get_dist=0')
h('dlc=0')  # dilate centres.
structural_plotting=1
h('n_cell=50')
n_cell=h.n_cell
h('numcell=n_cell')  # should really use sed. Such as to reduce duplicate variable to one.
n_cell=h.n_cell

#I normaly use alphabet size of two symbols (0,1). There is no reason why this could not be increased. Using multilevel quantisation.
prunenet=0
h('qd=0')#quick and dirty for debugging
h('lpy=1')#declare lpy here, the value is not important, its only important that the declaration exists.
h('if(qd==1){prunenet=900}')
h('objref py')
#h('objref py')
#h('py = new PythonObject()')
#h('if(qd==1){ py.structural_plotting=0')
h('minr=10')#minumum r large, for some reason which I do not yet 
h('minrs=10')#minumum r small, understand
#the number of connections made, depends somewhat on the last run, this may be 
# caused by morph3.hoc, deleting inappropriate cells occasionally.
# h("offset=0")
##OLD
h('offset=0')
##
h('parallels=0')
# h("tstop=1000//1 second =1000 ms")
h('delay=3000')
h('iw=1000 // divisor/ or scaler, so if its bigger weight should be smaller. Was 1000 but this should be stronger. delay equals internal weight.'
  )
h('iwp=1')
h('ew=0.0185// External weight.This value very appropriate. But perhaps smaller how to automate the tuning of synaptic weights?'
  )
h('run_iter=1')
h('msdp=0 //basically prune some connections distal from the soma.')
h('fast_wire=0')
h('ff=0')
h('if(ff==1){ get_dist=0 }')
h('if(ff==1){ n_cell=40 }')
h('if(ff==1){ prunenet=0 }')
h('no_b=0')
h('if (numcell>99){large=1}')
h('if (numcell>99){large_scale=1}')
h('{load_file("stdgui.hoc")}')
h('{load_file("nrngui.hoc")}')

h('load_file("import3d.hoc")')
h('strdef workingdir')
h('workingdir = getcwd()')
h('strdef neurons_home')
h('neurons_home = neuronhome()')
h('strdef machinename')
h('machine_name(machinename)')
h('strdef worm')
h('strdef morphrat')
h('strdef morphhuman')
h('strdef graph_dir')
h('strdef results_dir')
h('strdef matrix_dir')
h('sprint(matrix_dir,"%s%s",workingdir,"matrix_dir")')

h('sprint(results_dir,"%s%s",workingdir,"results")')
h('sprint(graph_dir,"%s%s",workingdir,"graphs")')
h('sprint(morphrat,"%s%s",workingdir,"main")')
h('sprint(morphhuman,"%s%s",workingdir,"SWC-2013human")')
h('sprint(worm,"%s%s",workingdir,"/openworm/CNG version/")')


def stationary_poisson(nsyn,lambd,tstart,tstop):
    ''' Generates nsyn stationary possion processes with rate lambda between tstart and tstop'''
    interval_s = (tstop-tstart)*.001
    spiketimes = []
    for i in xrange(nsyn):
        spikecount = np.random.poisson(interval_s*lambd)
        spikevec = np.empty(spikecount)
        if spikecount==0:
            spiketimes.append(spikevec)
        else:
            spikevec = tstart + (tstop-tstart)*np.random.random(spikecount)
            spiketimes.append(np.sort(spikevec)) #sort them too!

    return spiketimes


#Fetch Mainen&Sejnowski 1996 model files
if not os.path.isfile('/home/zaza3/trunk/examples/patdemo/cells/j4a.hoc') and RANK==0:
    #get the model files:
    u = urllib2.urlopen('http://senselab.med.yale.edu/ModelDB/eavBinDown.asp?o=2488&a=23&mime=application/zip')
    localFile = open('patdemo.zip', 'w')
    localFile.write(u.read())
    localFile.close()
    #unzip:
    myzip = zipfile.ZipFile('patdemo.zip', 'r')
    myzip.extractall('.')
    myzip.close()

#resync MPI threads
COMM.Barrier()

# Define cell parameters

# Define synapse parameters
synapse_parameters = {
    'idx' : 0, # to be set later
    'e' : 0.,                   # reversal potential
    'syntype' : 'AMPA',       # synapse type
    'tau' : 5.,                 # syn. time constant
    'weight' : .001,            # syn. weight
    'record_current' : True,
}

# Define electrode parameters
point_electrode_parameters = {
    'sigma' : 0.3,      # extracellular conductivity
    'x' : 0.,  # electrode requires 1d vector of positions
    'y' : 0.,
    'z' : 0.,
}


#number of units



n_cells = int(n_cell)#SIZE
cell_id = RANK
#execfile('/home/zaza3/Downloads/trunk/examples/expericomp20140421/atlas_reader3.py')

allrows=pickle.load(open('allrows.p','rb'))
h.xopen("/home/zaza3/Downloads/trunk/examples/expericomp20140421/morph4.hoc")
h.xopen("/home/zaza3/Downloads/trunk/examples/expericomp20140421/nqs.hoc")

os.system('pwd')
import re


os.chdir('/home/zaza3/Downloads/trunk/examples/expericomp20140421/')
#h.xopen("find_distances_.hoc")

#set the numpy random seeds
global_seed = 1234
np.random.seed(global_seed)


#assign cell positions
x_cell_pos = np.linspace(-250., 250., n_cells)

z_rotation = np.random.permutation(np.arange(0., np.pi, np.pi / n_cells))

#synaptic spike times
n_pre_syn = 1000
pre_syn_sptimes = stationary_poisson(nsyn=n_pre_syn, lambd=5., tstart=0, tstop=1000)

#re-seed the random number generator
cell_seed = global_seed + cell_id
np.random.seed(cell_seed)
cnt1=pc.id
#fn = re.compile('*.swc')
os.chdir('/home/zaza3/Downloads/trunk/examples/expericomp20140421/main/')
h('objref cell')
h('objref cells')
h('cells = new List()')
h('objref gidvec')
h('gidvec =new Vector()')


#cells=0
#cells=[cells]*40
#cells=[ LFPy.Cell(**cell_parameters) for i in range(int(n_cell+1))];
cells=[ 0 for i in range(int(n_cell+1))];


color_vec = [plt.cm.rainbow(int(x*256./n_cells)) for x in xrange(int(len(allrows)))]

#cells=[for x in xrange(0,50)]
#This is to make every cell a different color.


h('objref py')
h('py = new PythonObject()')
cnt=0 #This value cnt, will be different on different hosts and thats what is required.
for i in range(RANK, int(len(allrows)), SIZE): #20 was int(len(allrows))
    if((cnt<n_cell-1) and (0<i)):#While there are less than 100 cells on every node.
        if(cnt<1):    
            if(int(len(s))>9):#//This condition is counter intuitive many cells are         
                s=allrows[i]    
                storename= str(s[3]) #//simply being inside a loop, may be the main pro
                if(re.search('.swc',storename)):
                                        
                    h.cell=h.mkcell(storename)
                    h('cell.pyr()')
                    h('strdef cellposition')
                    h('sprint(cellposition,"%s%d%s%d%s%d%s","cell.position(",py.int(py.s[0]),",",py.int(py.s[1]),",",py.int(py.s[2]),")")')
                    h('execute(cellposition)')
	
                    h('pc.set_gid2node(py.i, pc.id)')# // associate gid i with this host
                    h('objref nc')
                    h('nc = cell.connect2target(nil)')# // attach spike detector to cell
                    h('pc.cell(py.i, nc)')#//



                    h('cells.append(cell)')
                    h('gidvec.append(py.i)')
                    h('print py.i')

                    cnt+=1
                
             
h('forall{ for(x,0){ insert xtra }}') 
h('forall{ for(x,0){ insert extracellular}}') 

h('chdir("../")')
h('xopen("interpxyz.hoc")')
h('grindaway()')
#h('forall{ for(x,0){ print x_xtra }}')


h('objref coords')
h('coords = new Vector(5)')


h('objref coords2')
h('coords2 = new Vector(3)')

#print coords, RANK

"""
if RANK == 0:
    #data = np.arange(1000, dtype='f')
    COMM.Send([coords, MPI.FLOAT], dest=1, tag=77)
elif RANK == 1:
    #data = np.empty(1000, dtype='f')
    COMM.Recv([coords, MPI.FLOAT], source=0, tag=77)
    print "Message Received, data is: ", coords
if RANK == 0:
    coords=np.array(h.coords.to_python(),dtype=np.float32)

else:
    coords = None
data = COMM.bcast(coords, root=0)
print "bcast finished and data on RANK %d is: "%RANK, coords
"""
h('objref syn_')
h('objref synlist')
h('synlist= new List()')

COMM.Barrier()#wait until all hosts get to this point

h('objref cell2')
h('cell2=pc.gid2cell(gidvec.x[0])')
#h('forsec cell2.all{psection()}')
#h('for i=0,gidvec.size-1{ \
"""
h('forsec cell.all{ \
    for(x){\
		if((x>0)&&(x<1)){ \
                   print x_xtra(x), "got here"\
//coords.x[0]=x_xtra(x) \
//                   coords.x[1]=y_xtra(y) \
//                   coords.x[2]=z_xtra(z) \
//                   coords.x[3]=i \
//                   coords.x[4]=x \
//                   coords.printf \
} \
}\
                  } \
}')
"""
#print 'got here'
#for i in xrange(0,SIZE):
for i in xrange(0,SIZE):
    if i==RANK:
        h('coords.x[0]=cells.o(0).soma[0].x_xtra()')
        h('coords.x[1]=cells.o(0).soma[0].y_xtra()')
        h('coords.x[2]=cells.o(0).soma[0].z_xtra()')
        #h('coords.x[4]=x')
        h('coords.x[3]=gidvec.x[0]')

        coords=np.array(h.coords.to_python(),dtype=np.float32)
    else: 
        coords= None
    data = COMM.bcast(coords, root=RANK)
    #COMM.Barrier() 
    if i==RANK:
        h('strdef secnames')
        h('cells.o(0).soma[0]{ secnames=secname() }')
        secnames=h.secnames
        #h('py.secnames=py.str(secnames) }')

    else: 
        secnames= None
    data = COMM.bcast(secnames, root=RANK)
    #COMM.Barrier() 
    #The secnames of the post synapses, and coordinates have been broadcast.
    #Because its a broadcast. One statement sends and recieves. Its implied that all hosts will recieve.
    if secnames!=None:
        if coords!=None:
            print "bcast finished and data on RANK %d is: "%RANK, coords, secnames
            h('coords2.x[0]=cells.o(0).soma[0].x_xtra()')
            h('coords2.x[1]=cells.o(0).soma[0].y_xtra()')
            h('coords2.x[2]=cells.o(0).soma[0].z_xtra()')
            h('print coords2.x[2]')
            h('coordsx=0.0')
            h.coordsx=coords[0]
            h('coordsy=0.0')
            h.coordsy=coords[1]
            h('coordsz=0.0')
            h.coordsz=coords[2]

            r=0.
            h('r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)')
            r=float(h.r)
            print r
            h('print r, "this is r"')
            #Assign a post synapse.
            post_syn=secnames+' '+'syn_ = new AMPA(0.5)'
            h(post_syn)
            #place the post synapse.
            h('synlist.append (syn_)')
            h('nc = pc.gid_connect(coords.x[3], syn_)')#connect the pre synapse to the post synapse.
            h('print coords.x[3], "this is the global identifier"')
            h('nc.delay = 1')
            h('nc.weight = 0.01')

        


"""
for i in range(RANK, int(len(allrows)), SIZE):
    if rank == i:
        for i_proc in xrange(1, SIZE):
            comm.Send([coords, MPI.FLOAT], dest=i_proc, tag=77)
    elif rank == !i:
    
        comm.Recv([coords, MPI.FLOAT], source=0, tag=77)
        print "Message Received, data is: ", coords


#    if RANK==i:#If the host 0, send stimulation from electrode to all other hosts.
#        for i_proc in xrange(1, SIZE):
#            COMM.recv(   ,source=i_proc)#COMM.r
#        COMM.send(,dest=0) #if not on the 0th host convey all the information LFP to the first host.

zips = []
i=0
for i in xrange(0,int(len(cells)-1)):
    print i
    print cells[i], RANK
    #if(i>1):
    if((i<cnt-1) and (0<i)):#While there are less than 100 cells on every node.
        
        for x, z in cells[1].get_idx_polygons():
            #print len(cells), 'length of cells', np.shape(cells)
        #print i
            zips.append(zip(x, z))
        #This is not a 3d plot as there is no y.



linecol = LineCollection(zips,
                         edgecolor = 'none',
                         facecolor = color_vec[i],
                         rasterized=False,
                         )            


fig = plt.figure()#figsize=(12, 8))
# Morphologies axes:
plt.axes([.175, .0, .65, 1], aspect='equal')
plt.axis('off')

ax = plt.gca()
ax.add_collection(linecol)
axis = ax.axis(ax.axis('equal'))
ax.axis(np.array(axis) / 1.15)
axis = ax.axis(ax.axis('equal'))
ax.axis(np.array(axis) / 1.15)
save_str='example_me'+str(RANK)+'.png'
os.chdir('/home/zaza3/trunk/examples/expericomp20140421/')
fig.savefig(save_str, dpi=300)
#crash here please

#exit()    
#quit()    
# Create cell


#Have to position and rotate the cells!

#assign spike times to different units
n_synapses = 100

# Create synapse and set time of synaptic input
pre_syn_pick = np.random.permutation(np.arange(n_pre_syn))[0:n_synapses]

for i_syn in xrange(n_synapses):
    syn_idx = int(cell.get_rand_idx_area_norm())
    synapse_parameters.update({'idx' : syn_idx})
    synapse = LFPy.Synapse(cell, **synapse_parameters)
    synapse.set_spike_times(pre_syn_sptimes[pre_syn_pick[i_syn]])

#run the cell simulation
cell.simulate(rec_imem=True,rec_isyn=True)

#set up the extracellular device
point_electrode = LFPy.RecExtElectrode(cell, **point_electrode_parameters)
point_electrode.calc_lfp()

if RANK==0:#If the host 0, send stimulation from electrode to all other hosts.
    single_LFPs = [point_electrode.LFP[0]]
    for i_proc in xrange(1, SIZE):
        single_LFPs = np.r_['0,2', single_LFPs, COMM.recv(source=i_proc)]#COMM.recv 
        #Facilitates point to point communication.
        #The results of LFPs are conveyed from all other hosts to the 0th host, where they are collated.
        
else:
    COMM.send(point_electrode.LFP[0], dest=0) #if not on the 0th host convey all the information LFP to the first host.

#MPI_Reduce
#Reduces values on all processes to a single value (by summing them in parallel) 
# we can also use MPI to sum arrays directly:
summed_LFP = COMM.reduce(point_electrode.LFP[0])


if RANK==0:
    #assign color to each unit
    color_vec = [plt.cm.rainbow(int(x*256./n_cells)) for x in xrange(n_cells)]

    #figure
    fig = plt.figure(figsize=(12, 8))
    
    # Morphologies axes:
    plt.axes([.175, .0, .65, 1], aspect='equal')
    plt.axis('off')

    for i_cell in xrange(n_cells):
        cell = LFPy.Cell(cell_parameters['morphology'],
                         nsegs_method='lambda_f',
                         lambda_f=5)
        cell.set_rotation(x=4.99, y=-4.33, z=z_rotation[i_cell])
        cell.set_pos(xpos=x_cell_pos[i_cell])

        zips = []
        for x, z in cell.get_idx_polygons():
            zips.append(zip(x, z))
            #These look like plot arguments.
            linecol = LineCollection(zips,
                    edgecolor = 'none',
                    facecolor = color_vec[i_cell],
                    rasterized=False,
                    )            

        ax = plt.gca()
        ax.add_collection(linecol)
    
    axis = ax.axis(ax.axis('equal'))
    ax.axis(np.array(axis) / 1.15)


    #adding a blue dot:
    ax.plot(point_electrode.x, point_electrode.z, 'o',
            markeredgecolor='none', markerfacecolor='b', markersize=3,
            zorder=10, clip_on=False)
    plt.annotate("Electrode",\
            xy=(0., 0.), xycoords='data',\
            xytext=(-100., 1000.),
            arrowprops=dict(arrowstyle='wedge',
                            shrinkA=1,
                            shrinkB=1,
                            #lw=0.5,
                            mutation_scale=20,
                            fc="0.6", ec="none",
                            edgecolor='k', facecolor='w'))

    plt.xlim([-700., 700.])

    ax.plot([100, 200], [-250, -250], 'k', lw=1, clip_on=False)
    ax.text(150, -300, r'100$\mu$m', va='center', ha='center')

    #presynaptic spike trains axes
    plt.axes([.05, .35, .25, .55])

    pop_sptimes = []
    for i_pre in xrange(n_pre_syn):
        sp = pre_syn_sptimes[i_pre]
        for i_sp in xrange(len(sp)):
            pop_sptimes.append(sp[i_sp])
               
    for i_pre in xrange(n_pre_syn):
        plt.scatter(pre_syn_sptimes[i_pre],
                    i_pre*np.ones(len(pre_syn_sptimes[i_pre])),
                    s=1, edgecolors='none', facecolors='k')

    plt.ylim([0,n_pre_syn])
    plt.xlim([0,cell_parameters['tstopms']])
    plt.ylabel('train #', ha='left')
    plt.title('Presynaptic spike times')
    
    ax = plt.gca()
    for loc, spine in ax.spines.iteritems():
        if loc in ['right', 'top']:
            spine.set_color('none')            
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
    ax.set_xticklabels([])

    #spike rate axes
    plt.axes([.05,.12,.25,.2])

    binsize = 5
    bins=np.arange(0, cell_parameters['tstopms']+1., binsize)
    count,b = np.histogram(pop_sptimes, bins=bins)
    rate = count*(1000./binsize)*(1./n_pre_syn)
    plt.plot(b[0:-1],rate,color='black',lw=1)

    plt.xlim([0,cell_parameters['tstopms']])
    plt.ylim([0,10.])
    
    tvec = np.arange(point_electrode.LFP.shape[1])*cell.timeres_python 

    plt.xlabel('$t$ (ms)')
    plt.ylabel('Rate (spike/s)')
    
    ax = plt.gca()
    for loc, spine in ax.spines.iteritems():
        if loc in ['right', 'top']:
            spine.set_color('none')            
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    #single neuron EPs axes
    plt.axes([.7,.35,.25,.55])

    plt.title('Single neuron extracellular potentials')
    plt.axis('off')

    for i_cell in xrange(n_cells):
        plt.plot(tvec,
                        i_cell+2.e3*single_LFPs[i_cell],
                        color=color_vec[i_cell], lw=1,
                        )

    plt.ylim([-1,n_cells-.5])

    #Summed LFPs axes
    plt.axes([.7,.12,.25,.2])
    plt.plot(tvec, 1E3*summed_LFP, color='black', lw=1)
    plt.ylim([-5.e-1,5e-1])

    plt.title('Summed extracellular potentials')
    plt.xlabel(r'$t$ (ms)')
    plt.ylabel(r'$\mu$V',ha='left',rotation='horizontal')

    ax = plt.gca()
    for loc, spine in ax.spines.iteritems():
        if loc in ['right', 'top']:
            spine.set_color('none')            
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    save_str='example3.png'+str(RANK)
    fig.savefig(save_str, dpi=300)
objref cell
objref gidvec  // to associate gid and position in cells List
proc mbsp(){local i,j,k, cell_diam, cell_length, r  localobj nc, nil, rs, ranlist, morph, import
//, gidvecgidvec = 
   gidvec = new Vector()


   chdir(workingdir)
   nrnpython("from neuron import h")
	
   nrnpython("execfile('atlas_reader3.py')") 
   nrnpython("import cPickle")
   nrnpython("names = cPickle.load(open('names.p', 'rb'));")
   chdir(morphrat)
   //nrnpython("fg = open('goodindexfile', 'a')")
   //nrnpython("fb = open('badindexfile.sh', 'a')") //open badindexfile for both reading and writing.

   halfn_cell=$1/2
   n_cell = $1
   strdef storename
  
   cnt+=pc.id
   i=0
   
   for (i=pc.id; i < py.int(py.len(py.allrows)) ; i += pc.nhost) {   
   //for (i=pc.id; i <  py.int(py.len(py.allrows)) ; i += pc.nhost) {
    if(cnt<n_cell+1){

 
    //print i
    py.s=py.allrows[i]
    
    storename= py.str(py.s[3]) //simply being inside a loop, may be the main problem. 
    
   
   if(int(py.len(py.s))>9){ //This condition is counter intuitive many cells are getting created, which subsequentely cannot pass 
     //n_cell=10
     //print a, b, c  
     //if(cells.count<n_cell){
  
       //If there are enough columns in the file. Ie there should be nine not 7. Updated from 5
        //print i

        //print storename

        chdir(morphrat)
     
        nrnpython("i=int(h.i)")
        nrnpython("i=int(i)")
        nrnpython("print type(i)")

        nrnpython("sx=str(names[int(i)]+'.CNG.swc')") 
        //If the line below does not execute the bad index is stored in bif.
        //cell = mkcell(py.sx)
        print i
        
        cell=mkcell(storename)  
        cell.name=storename 
        //cell.nametype=py.str(py.s[5]) //assign the char value of cell type
        //cell.num_type=py.int(py.s[6])
        //cell.population=py.str(py.s[7])
      
        cell.pyr()
        //access cell.soma[0]
        //psection()
        //forsec cell.all{ insert hh }
        //cell = new Cellbasic()
         //forsec cell.basal{ print secname() }
         //object_id(cell.basal)!=0 //cell has been cast to a double, but how.
        
        //The expecting a double message, genuinely means that the cell has been cast to a double some how.
        
   
        //The append error is something to do with the cell template, but something after p3dadd remember debugging messages are misleading.
        print cell
        cells.append(cell)
        
        pc.set_gid2node(cnt, pc.id) // associate gid i with this host

        nc=cell.connect2target(nil)
        pc.cell(cnt, nc)

        
        gidvec.append(cnt)
        cnt+=pc.nhost

        

        chdir(morphrat)
 
      //}
      }
           
       //File names required not repo names. 
    }   
  }

}
pc.barrier()


"""
