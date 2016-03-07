'''
Author information. This is an extension to the Utils class from Allen Brain API. 
The parallel wiring related functions are written by Russell Jarvis russell_jarvis@riseup.net
https://www.linkedin.com/in/russell-jarvis-02433a30?trk=hp-identity-name
'''


from allensdk.model.biophys_sim.neuron.hoc_utils import HocUtils
import logging
import glob
from mpi4py import MPI
import numpy as np
import logging
import networkx
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(funcName)s - %(lineno)d')
fh = logging.FileHandler('wiring.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)
import unittest
import pdb as pdb

class Utils(HocUtils):#search multiple inheritance unittest.
    _log = logging.getLogger(__name__)
    def __init__(self, description, NCELL=20,readin=0):
        super(Utils, self).__init__(description)
        h=self.h  
        h('objref pc, py')
        h('pc = new ParallelContext()')
        h('py = new PythonObject()')
        self.readin=readin
        self.synapse_list=[]
        self.stim = None
        self.stim_curr = None
        self.sampling_rate = None
        self.cells = []
        self.gidlist=[]
        self.NCELL=NCELL
        self.celldict={}
        self.COMM = MPI.COMM_WORLD
        self.SIZE = self.COMM.Get_size()
        self.RANK = self.COMM.Get_rank()
        self.allsecs=None #global list containing all NEURON sections, initialized via mkallsecs
        self.coordict=None
        self.celldict={}
        self.cellmorphdict={}
        self.nclist = []
        self.seclists=[]
        self.tvec=self.h.Vector()    
        self.idvec=self.h.Vector() 
        self.icm = np.zeros((self.NCELL, self.NCELL))
        self.ecm = np.zeros((self.NCELL, self.NCELL))
        self.visited = np.zeros((self.NCELL, self.NCELL))
        self.ecg = networkx.DiGraph()
        self.icg = networkx.DiGraph()
        self.my_visited = np.zeros_like(self.icm)
        self.my_icm = np.zeros_like(self.icm)
        self.my_ecm = np.zeros_like(self.ecm)
        self.my_ecg = networkx.DiGraph()
        self.my_icg = networkx.DiGraph()
        self.debugdata=[]

    def prep_list(self):                    
        '''
        find which list has the shortest length.
        and construct a new list with 1 in 3 inhibitory neurons and 2 out of 3 
         excitatory neurons. 
        It would be preferable to make an exhaustive list of all neurons
        however this is not practical for debugging small models, composed
        of a balance between excitation and inhibition.
        '''
        import pickle
        allrows = pickle.load(open('allrows.p', 'rb'))
        allrows.remove(allrows[0])#The first list element are the column titles. 
        allrows = [i for i in allrows if int(len(i))>9 ]
        excitatory = [i for i in allrows if i[5]!="interneuron" ]        
        interneurons = [i for i in allrows if i[5]=="interneuron" ]     
        bothtrans=[]
        if len(excitatory) > len(interneurons):
            length=len(interneurons)
        else:
            length=len(excitatory)
        for i in xrange(0,length):
            #If change conditions, then must re-run wiring not from file
            #if ((i%2)==0):
            if ((i%3)==0): #2/3 excitatory to reflect cortical balance of transmitters.
                bothtrans.append(interneurons[i]) 
            else:
                bothtrans.append(excitatory[i])
        return bothtrans
    
    def read_local_swc(self):
        h=self.h    
        NCELL=self.NCELL
        SIZE=self.SIZE
        RANK=self.RANK
        pc=h.ParallelContext()
        morphs=[]
        cells1=[]
        self.initialize_hoc()
        swclist=glob.glob('*.swc')
        itergids = iter( i for i in range(RANK, len(swclist), SIZE) )
        for i in itergids:
            cell = h.mkcell(swclist[i])
            self.generate_morphology(cell,swclist[i])
            self.load_cell_parameters(cell, fit_ids[utils.cells_data[i]['type']])
            cells1.append(cell1)
            morphology.root
            morphs.append(morphology)
        return morphs,swclist,cells1

    #TODO use neuro electro to test cortical pyramidal cells, and baskett cells before including
    #them in the network.
    #Call a method test_cell inside the make_cells function.
    #def test_cell(self,type):
    #    from neuron.unit.neuroelectro import NeuroElectroSummary
    #    if type=='hippocampus':
    #        summary = NeuroElectroSummary(neuron={'name':'Hippocampus CA1 Pyramidal Cell'},
    #                                    ephysprop={'name':'spike width'})
    #        observation = summary.get_observation(show=True)
    #        from neuronunit.tests import import SpikeWidthTest
    #        ca1_pyramdical_spike_width_test=SPikeWidthTest(observation=observation)
            
            
    def make_cells(self,polarity):
        h=self.h    
        NCELL=self.NCELL
        SIZE=self.SIZE
        RANK=self.RANK
        tvec=h.Vector()
        gidvec=h.Vector()
        h('objref tvec, gidvec')
        h('gidvec = new Vector()')
        h('tvec = new Vector()')
        #self.testcell(type)
        pc=h.ParallelContext()
        d = { x: y for x,y in enumerate(polarity)}         
        itergids = iter( i for i in range(RANK, NCELL, SIZE) )        
        #TODO keep rank0 free of cells, such that all the memory associated with that CPU is free for graph theory related objects.
        #itergids = iter( i for i in range(RANK+1, NCELL, SIZE-1) )        
        fit_ids = self.description.data['fit_ids'][0] #excitatory        
        for i in itergids:
            cell = h.mkcell(d[i][3])            
            print cell, d[i][3]
            cell.geom_nseg()
            cell.gid1=i #itergids.next()
            
            #excitatory cell.
            if 'pyramid' in d[i]:                
                self.load_cell_parameters(cell, fit_ids[self.cells_data[0]['type']])
                cell.polarity=1
                #TODO use neuroelectro here, to unit test each cell and to check if it will fire.
            else:            
                #inhibitory cell.
                #TODO use neuroelectro here, to unit test each cell and to check if it will fire.
                self.load_cell_parameters(cell, fit_ids[self.cells_data[2]['type']])
                cell.polarity=0           
            h('Cell[0].soma[0] nc =  new NetCon(&v(0.5), nil)')                        
            pc.set_gid2node(i,RANK)
            h('pc.cell('+str(i)+', nc)')
            hocstring='pc.spike_record('+str(i)+',tvec,gidvec)'
            h(hocstring)
            pc.spike_record(i,tvec,gidvec)
            
            cell1=pc.gid2cell(i)
            self.celldict[i]=cell
            self.cells.append(cell)
    
    def gcs(self,NCELL):
        """Instantiate NEURON cell Objects in the Python variable space, such that cell
        that all cells have unique identifiers."""
        NCELL=self.NCELL
        SIZE=self.SIZE
        RANK=self.RANK
        h=self.h    
        pc=h.ParallelContext()     
        h('objref nc, cells')
        swcdict={}
        NFILE = 3175
        fit_ids = self.description.data['fit_ids'][0] #excitatory        
        self.cells_data = self.description.data['biophys'][0]['cells']
        bothtrans =self.prep_list()    
        import os
        os.chdir(os.getcwd() + '/main')  
        self.make_cells(bothtrans)
        len(h.List('NetCon'))                        
        pol=[ a.polarity for a in self.cells ]       
        os.chdir(os.getcwd() + '/../')               
        self.h.define_shape()        
        self.h('forall{ for(x,0){ insert xtra }}')
        self.h('forall{ for(x,0){ insert extracellular}}')    
        self.h('xopen("interpxyz.hoc")')
        self.h('grindaway()')    
    
    def matrix_reduce(self, matrix=None):
        # TODO make this method argument based so it can handle arbitary input matrices not one in particular
        #
        import networkx as nx
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        import numpy as np

        ##Argument matrix
        if matrix!=None:
            self.my_matrix = np.zeros_like(self.matrix)
            COMM.Reduce([self.matrix, MPI.DOUBLE], [self.my_matrix, MPI.DOUBLE], op=MPI.SUM,
                        root=0)
            
        ##Standard matrices that will always need to be reduced.    
        COMM.Barrier()
        self.my_visited = np.zeros_like(self.visited)
        COMM.Reduce([self.visited, MPI.DOUBLE], [self.my_visited, MPI.DOUBLE], op=MPI.SUM,
                    root=0)

        self.my_icm = np.zeros_like(self.icm)
        COMM.Reduce([self.icm, MPI.DOUBLE], [self.my_icm, MPI.DOUBLE], op=MPI.SUM,
                    root=0)
        self.my_ecm = np.zeros_like(self.ecm)
        COMM.Reduce([self.ecm, MPI.DOUBLE], [self.my_ecm, MPI.DOUBLE], op=MPI.SUM,
                    root=0)
        self.my_icg = nx.DiGraph(self.my_icm)
        self.my_ecg = nx.DiGraph(self.my_ecg)
        if RANK==0:
            assert np.sum(self.my_ecm)!=0
        if RANK==0:
            assert np.sum(self.my_icm)!=0
        #TODO check if visited is empty.
        #if RANK==0:
        #    assert np.sum(self.my_visited)!=0
        

    def vec_reduce(self):       
        print type(np.array(self.idvec.to_python()))
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        COMM.Barrier()
        self.my_tvec = np.zeros_like(self.idvec.to_python())
        COMM.Reduce([np.array(self.tvec.to_python()), MPI.DOUBLE], [self.my_tvec, MPI.DOUBLE], op=MPI.SUM,
                    root=0)
        self.my_idvec = np.zeros_like(self.tvec.to_python())
        COMM.Reduce([np.array(self.idvec.to_python()), MPI.DOUBLE], [self.my_idvec, MPI.DOUBLE], op=MPI.SUM,
                    root=0
        )


    def prun(self,tstop):
        h=self.h    
        pc=h.ParallelContext()
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        checkpoint_interval = 50000.

        #The following code is from the open source code at:
        #http://senselab.med.yale.edu/ModelDB/ShowModel.asp?model=151681
        #with some minor modifications
        cvode = h.CVode()
        cvode.cache_efficient(1)
        # pc.spike_compress(0,0,1)
    
        pc.setup_transfer()
        mindelay = pc.set_maxstep(10)
        if RANK == 0:
            print 'mindelay = %g' % mindelay
        runtime = h.startsw()
        exchtime = pc.wait_time()
    
        inittime = h.startsw()
        h.stdinit()
        inittime = h.startsw() - inittime
        if RANK == 0:
            print 'init time = %g' % inittime
    
        while h.t < tstop:
            told = h.t
            tnext = h.t + checkpoint_interval
            if tnext > tstop:
                tnext = tstop
            pc.psolve(tnext)
            if h.t == told:
                if RANK == 0:
                    print 'psolve did not advance time from t=%.20g to tnext=%.20g\n' \
                        % (h.t, tnext)
                break    
            print 'working', h.t
        runtime = h.startsw() - runtime
        comptime = pc.step_time()
        splittime = pc.vtransfer_time(1)
        gaptime = pc.vtransfer_time()
        exchtime = pc.wait_time() - exchtime
        if RANK == 0:
            print 'runtime = %g' % runtime
        print comptime, exchtime, splittime, gaptime


    def nestedpre(self,j):

        h=self.h   
        pc=h.ParallelContext()
        shiplist=[]
        h('objref coords') 
        h('coords = new Vector(3)')
        #self.celldict.items()[0]
        if j in self.celldict.keys():
            seglist= iter( (seg, sec, self.celldict[j]) for sec in self.celldict[j].spk_trig_ls for seg in sec )     
            for (seg,sec, cellc) in seglist:
                sec.push()
                get_cox = str('coords.x[0]=x_xtra('
                              + str(seg.x) + ')')
                h(get_cox)                   
                get_coy = str('coords.x[1]=y_xtra('
                              + str(seg.x) + ')')
                h(get_coy)
                get_coz = str('coords.x[2]=z_xtra('
                              + str(seg.x) + ')')
                h(get_coz)
                coordict={} #Destroy dictionary some-how.
                            #This is a pointer problem.
                coordict['hostfrom'] = pc.id()
                coordict['coords'] = np.array(h.coords.to_python(),
                                          dtype=np.float64)
                coordict['gid']= int(j)
                coordict['seg']= seg.x                    
                secnames = self.h.cas().name()#sec.name()  
                coordict['secnames'] = str(secnames)
                shiplist.append(coordict)
                self.h.pop_section()
               
        return shiplist
        
    def alloc_synapse_ff(self,r,post_syn,cellind,k,gidn,i):

        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        #from neuron import h
        h=self.h  
        
        pc=h.ParallelContext()
        polarity = 0        
        polarity=int(h.Cell[int(cellind)].polarity)
        if polarity==1:
            #TODO pickle load the graphs here instead of making them manually.
            self.ecg.add_edge(i,gidn,weight=r/0.4)
            self.ecm[i][gidn] = self.ecm[i][gidn] + 1
            assert np.sum(self.ecm)!=0
        else:
            self.icg.add_edge(i,gidn,weight=r/0.4)
            self.icm[i][gidn] = self.icm[i][gidn] + 1
            self.icg.add_edge(i,gidn,weight=r/0.4)
            assert np.sum(self.icm)!=0                
            #TODO Add other edge attributes like secnames etc.
        h('objref syn_')   
        h(post_syn)
        syn_=h.syn_
        h.syn_.cid=i
        h.Cell[cellind].ampalist.append(h.syn_)
        h.Cell[cellind].div.append(k['gid'])
        h.Cell[cellind].gvpre.append(k['gid'])
        nc=pc.gid_connect(k['gid'],syn_)                                        
        nc.threshold = -20
        nc.delay=1+r/0.4
        nc.weight[0]=r/0.4    
        self.nclist.append(nc)
        

    def alloc_synapse(self,r,h,sec,seg,cellind,secnames,k,i,gidn):
        
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        from neuron import h
        pc=h.ParallelContext()
        h=self.h
        self.visited[i][gidn] = self.visited[i][gidn] + 1              
        if r < 2.5: #2.5 micro metres.
            polarity = 0        
            polarity=int(h.Cell[int(cellind)].polarity)
            h('objref syn_')        
            if int(polarity) == int(0):
                post_syn = secnames + ' ' + 'syn_ = new FastInhib(' + str(seg.x) + ')'
                #post_syn = secnames + ' ' + 'syn_ = new GABAa(' + str(seg.x) + ')'
                self.icm[i][gidn] = self.icm[i][gidn] + 1
                self.icg.add_edge(i,gidn,weight=r/0.4)
                self.icg[i][gidn]['post_loc']=secnames
                self.icg[i][gidn]['pre_loc']=k['secnames']

                assert np.sum(self.icm)!=0
        
            else:
                if (k['gid']%2==0):
                    post_syn = secnames + ' ' + 'syn_ = new ExpSid(' + str(seg.x) + ')'                       
                    self.ecm[i][gidn] = self.ecm[i][gidn] + 1
                    self.ecg.add_edge(i,gidn,weight=r/0.4)
                    self.ecg[i][gidn]['post_loc']=secnames
                    self.ecg[i][gidn]['pre_loc']=k['secnames']
                    self.seclists.append(secnames)
                    assert np.sum(self.ecm)!=0
                else:
                    #TODO Find standard open source brain affiliated code for NMDA synapse
                    post_syn = secnames + ' ' + 'syn_ = new ExpSid(' + str(seg.x) + ')'                       
                    self.ecm[i][gidn] = self.ecm[i][gidn] + 1
                    self.ecg.add_edge(i,gidn,weight=r/0.4)
                    self.ecg[i][gidn]['post_loc']=secnames
                    self.ecg[i][gidn]['pre_loc']=k['secnames']
                    self.seclists.append(secnames)
                    assert np.sum(self.ecm)!=0
        
            h(post_syn)
            self.synapse_list.append((r,post_syn,cellind,k,gidn,i))
            syn_=h.syn_
            h.syn_.cid=i
            h.Cell[cellind].ampalist.append(h.syn_)
            h.Cell[cellind].div.append(k['gid'])
            h.Cell[cellind].gvpre.append(k['gid'])
            nc=pc.gid_connect(k['gid'],syn_)                                        
            nc.threshold = -20
            nc.delay=1+r/0.4
            nc.weight[0]=r/0.4    
            self.nclist.append(nc)

        
       
    def nestedpost(self,data):
        """
        
        This is the inner most loop of the parallel wiring algorithm.
        13For ever GID
        For every coordinate thats received from a broadcast.
        for i,t in self.celldict.iteritems():
        For ever GID thats on this host (in the dictionary)
        if i in self.celldict.keys():
        for k,i,t in iterdata :
        if the gids are not the same.
        Rule out self synapsing neurons (autopses), with the condition
        pre GID != post GID
        If the putative post synaptic gid exists on this CPU, the referen
        tree.
        My wiring algorithm uses HOC variables that interpolate the middl
        some C libraries to achieve collision detection for synapse alloc
        from neuromac.segment_distance import dist3D_segment_to_segment
        from segment_distance import dist3D_segment_to_segment
        and I am considering using them here also
        """
        #from segment_distance import dist3D_segment_to_segment
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        from neuron import h
        pc=h.ParallelContext()
        h=self.h    
        secnames = ''
        cellind =0 
        polarity = 0
        h('objref coords')
        h('coords = new Vector(3)')
        h('objref pc')
        h('pc = new ParallelContext()')        
        h('objref coords2') 
        h('coords2 = new Vector(3)')
        #self.debugdata=data
        for q,s in enumerate(data):
            
            def test_wiring(q,s,data):
                if q<len(data)-1:
                    left=(str(s[q]['secnames'])+str(s[q]['seg']))
                    right=(str(s[q+1]['secnames'])+str(s[q+1]['seg']))
                    assert left!=right
                         
            test_wiring(q,s,data)
            for t in s:
                k={} #The only point of this redundantvariable switching is to force the dictionary k to be redclared 
                k=t #such that it is not prevented from updating
                itercell= ( (i,t) for i,t in self.celldict.iteritems() if i in self.celldict.keys() if int(t.gid1) != int(k['gid']) )       
                for i,t in itercell :                          
                    iterseg=iter( (seg,sec) for sec in t.spk_rx_ls for seg in sec)               
                    for (seg,sec) in iterseg:
                        segxold=seg.x
                        h('objref cell1')
                        h('cell1=pc.gid2cell('+str(i)+')')
                        secnames = sec.name()
                        cellind = int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.
                        h(str('coords2.x[2]=') + str('z_xtra(')
                          + str(seg.x) + ')')
                        h(str('coords2.x[1]=') + str('y_xtra(')
                          + str(seg.x) + ')')
                        h(str('coords2.x[0]=') + str('x_xtra(')
                        + str(seg.x) + ')')
            
                        h('coordsx=0.0')
                        h.coordsx = k['coords'][0]
                        h('coordsy=0.0')
                        h.coordsy = k['coords'][1]
                        h('coordsz=0.0')
                        h.coordsz = k['coords'][2]  
            #h('coordsx') and coordsx are not tautolous. 
            #One is a variable in the HOC space, the other is in
            #coordsx from the Python space has been broadcast ov      
                        coordsx = float(k['coords'][0])
                        coordsy = float(k['coords'][1])
                        coordsz = float(k['coords'][2])
              
            #Find the euclidian distance between putative presynaptic segments, 
            #and putative post synaptic segments.    
            #If the euclidian distance is below an allowable threshold in micro 
            #meters, continue on with code responsible for assigning a 
            #synapse, and a netcon. Neurons parallel context class can handle the actual message passing associated with sending and receiving action potentials on different hosts.                               
              
              
                        r = 0.
                        import math
                        r=math.sqrt((h.coords2.x[0] - coordsx)**2+(h.coords2.x[1] - coordsy)**2+(h.coords2.x[2] - coordsz)**2)
                        gidn=k['gid']    
                        r = float(r)                      
                        self.alloc_synapse(r,h,sec,seg,cellind,secnames,k,i,gidn)



    def wirecells(self):
        """This function constitutes the outermost loop of the parallel wiring algor
        The function returns two adjacency matrices. One matrix whose elements are excitatory connections and another matrix of inhibitory connections"""
        #from segment_distance import dist3D_segment_to_segment
        import pickle
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        h=self.h    
        pc=h.ParallelContext()
        secnames = ''
        cellind =0 
        polarity = 0
        h('objref coords')
        h('coords = new Vector(3)')
        h('objref pc')
        h('pc = new ParallelContext()')        
        coordict=None
        coordictlist=[]

        #Iterate over all CPU ranks, iterate through all GIDs (global 
        #identifiers, stored in the python dictionary).
        if self.readin!=1:    
            for s in xrange(0, SIZE):
                #Synchronise processes here, all ranks must have finished receiving 
                #transmitted material before another transmission of the coordictlis begins, potentially
                #overwritting a coordictlist before it has been properly exhausted.
                COMM.barrier() #New line could save CPU but functional? Can be removed
                coordictlist=[]
                if COMM.rank==s: #New line could save CPU but functional? Can be removed
                    print 'begin creating message for transmision on rank ', COMM.rank
                    celliter= iter(i for i in self.celldict.keys())  
                    for i in celliter:  
                        cell1=pc.gid2cell(i)
                        coordictlist.append(self.nestedpre(i))
                    print 'end tx on rank ', COMM.rank

                #All CPUs wait here so that recieving only occurs when the rank==s has finished
                #tranmsitting              
                data = COMM.bcast(coordictlist, root=s)  # ie root = rank
                print 'checking for rx on rank ', COMM.rank
                if len(data) != 0:
                    print 'receieved rx on rank ', COMM.rank
                    self.nestedpost(data)
                    print 'using received message on rank ', COMM.rank
                    print len(data)

            print('finished wiring of connectivity\n')
            fname='synapse_list'+str(RANK)+'.p'
            assert len(self.synapse_list)!=0
            with open(fname, 'wb') as handle:
                pickle.dump(self.synapse_list, handle)
            fname='visited'+str(RANK)+'.p'
            with open(fname, 'wb') as handle:
                pickle.dump(self.visited,handle)
        else:
            fname='synapse_list'+str(RANK)+'.p'
            with open(fname, 'rb') as handle:
                self.synapse_list=pickle.load(handle)
                #for s in self.synapse_list:
                for (r,post_syn,cellind,k,gidn,i) in self.synapse_list:
                    self.alloc_synapse_ff(r,post_syn,cellind,k,gidn,i)
  
    def tracenet(self):
        '''
        This method does two things.
        1 Send a matrix to rank0.
        2 Do a local hub node computation.
        Ideally there should be too many neurons to properly visualise them in a network graph.
        Future design decision. Keep rank0 free of cells, such that it has the RAM to store
        big matrices.
        # Then destroy after sending to rank 0.
        # Maybe stimulate one cell per CPU.
        #

        '''
        ncsize=len(self.h.NetCon)
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        self.matrix_reduce()
        lsoftup=[]
        for i, j in enumerate(self.h.NetCon):
            if type(j)!=None:
                assert type(j)!=None
                srcind=int(j.srcgid())
                tgtind=int(j.postcell().gid1)
                print int(j.srcgid()),int(j.postcell().gid1),self.celldict[srcind],self.celldict[tgtind]
                lsoftup.append((int(utils.h.NetCon[i].srcgid()),int(utils.h.NetCon[i].postcell().gid1),utils.celldict[srcind],utils.celldict[tgtind]))
        return lsoftup
      
    def plotgraph(self):
        assert self.COMM.rank==0        
        import json
        import networkx as nx
        from networkx.readwrite import json_graph
        #import http_server  
        #G=nx.from_numpy_matrix(self.my_ecm)
        G=self.my_ecg
        d = json_graph.node_link_data(G)     
        #t = json_graph.tree_graph(G)       
        json.dump(d, open('js/excitatory_network.json','w'))
        G=self.my_icg
        #G=nx.from_numpy_matrix(self.my_icm)
        d = json_graph.node_link_data(G)     
        #t = json_graph.tree_graph(G)       
        json.dump(d, open('js/inhibitory_network.json','w'))

        print('Wrote node-link JSON data to js/network.json')
        # open URL in running web browser
        #http_server.load_url('force/force.html')


        

    def generate_morphology(self, cell, morph_filename):
        #This code is from the Allen Brain API examples.

        h = self.h
        swc = self.h.Import3d_SWC_read()
        swc.input(morph_filename)
        imprt = self.h.Import3d_GUI(swc, 0)
        h('execute("forall delete_section()",cell)')
        imprt.instantiate(cell)
        
      
        
        for seg in cell.soma[0]:
            seg.area()
        for sec in cell.allsec():
            sec.nseg = 1 + 2 * int(sec.L / 40)
        h.define_shape()
        
        #cell.simplify_axon()
        #for sec in cell.axonal:
        #    sec.L = 30
        #    sec.diam = 1
        #    sec.nseg = 1 + 2 * int(sec.L / 40)
        #cell.axon[0].connect(cell.soma[0], 0.5, 0)
        #cell.axon[1].connect(cell.axon[0], 1, 0)
    
    def load_cell_parameters(self, cell, type_index):
        #This code is from the Allen Brain API examples.
        h=self.h
        passive = self.description.data['fit'][type_index]['passive'][0]
        conditions = self.description.data['fit'][type_index]['conditions'][0]
        genome = self.description.data['fit'][type_index]['genome']

        # Set passive properties
        cm_dict = dict([(c['section'], c['cm']) for c in passive['cm']])
        for sec in cell.all:
            sec.Ra = passive['ra']
            sec.cm = cm_dict[sec.name().split(".")[1][:4]]
            sec.insert('pas')
            for seg in sec:
                seg.pas.e = passive["e_pas"]

        # Insert channels and set parameters
        for p in genome:
            sections = [s for s in cell.all if s.name().split(".")[1][:4] == p["section"]]
            for sec in sections:
                sec.push()
                if p["mechanism"] != "":
                    print p["mechanism"]
                    sec.insert(p["mechanism"])
                    h('print psection()')
                setattr(sec, p["name"], p["value"])
                self.h.pop_section()
        # Set reversal potentials
        for erev in conditions['erev']:
            sections = [s for s in cell.all if s.name().split(".")[1][:4] == erev["section"]]
            for sec in sections:
                sec.ena = erev["ena"]
                sec.ek = erev["ek"]


    def setup_iclamp_step(self, target_cell, amp, delay, dur):
        self.stim = self.h.IClamp(target_cell.soma[0](0.5))
        self.stim.amp = amp
        self.stim.delay = delay
        self.stim.dur = dur

    def record_values(self):
        vec = { "v": {}, #define a dictionary.
                "t": self.h.Vector() }
        for i, cell in enumerate(self.cells):
            vec["v"][int(cell.gid1)]=self.h.Vector()
            vec["v"][int(cell.gid1)].record(cell.soma[0](0.5)._ref_v)
        vec["t"].record(self.h._ref_t)
    
        return vec
  
    def spikerecord(self):   
        for cell in self.cells:
            self.h.pc.spike_record(int(cell.gid1), self.tvec, self.idvec)



