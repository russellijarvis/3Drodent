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
import neuroelectro


#import neuroelectro, and test each cell to see if it conforms to acceptable dynamics for the allen brain 
#ontology it is supposed to represent.
#


class Utils(HocUtils):#search multiple inheritance unittest.
    _log = logging.getLogger(__name__)
    def __init__(self, description, NCELL=20,readin=0):
        super(Utils, self).__init__(description)
        # Logically the pc, and py object attributes of the HOC object should be initialized here too.
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
        self.my_tvec = np.zeros(0)
        self.my_idvec = np.zeros(0)

        self.icm = np.zeros((self.NCELL, self.NCELL))
        self.ecm = np.zeros((self.NCELL, self.NCELL))
        self.ecg = networkx.Graph()
        self.icg = networkx.Graph()
        if self.RANK==0:        
            self.my_icm = np.zeros((self.NCELL, self.NCELL))
            self.my_ecm = np.zeros((self.NCELL, self.NCELL))

        #self.pc=h.pc

    def __del__(self):
        """
        AUTHORS:
        - THOMAS MCTAVISH (2010-11-04): initial version. Modified version of the
                project by Hines and Carnevale. (Hines M.L. and Carnevale N.T, 
                Translating network models to parallel hardware in NEURON,
                Journal of Neuroscience Methods 169 (2008) 425-455).
        In the case of multiple runs where NEURON is not quit and reloaded,
        we need to clear NEURON so that we can run a new network.
        There is a particular order that things need to be deleted:
        1) Any recording vectors
        2) pc.gid_clear(0)
        3) Destroy NetCons
        4) Destroy Cells
        """
        self.t_vec = [] # Must come before gid_clear
        self.id_vec = [] # Must come before gid_clear
        self.pc.gid_clear(0)
        self.nclist = []  # Synaptic NetCon list on this host
        self.stim = None
        self.cells = []          # Cells on this host
        # I do not know how to refer to relative paths in Python, 
        # the below emulates a call to a relative path.


    
   



    def register_gid(self, gid, source, section=None):
        """Register a global ID with the global `ParallelContext` instance."""
        ###print "registering gid %s to %s (section=%s)" % (gid, source, section)
        self.parallel_context.set_gid2node(gid, self.mpi_rank) # assign the gid to this node
        if is_point_process(source):
            nc = h.NetCon(source, None)                          # } associate the cell spike source
        else:
            nc = h.NetCon(source, None, sec=section)
        self.parallel_context.cell(gid, nc)                     # } with the gid (using a temporary NetCon)
        self.gid_sources.append(source) # gid_clear (in _State.reset()) will cause a
                                        # segmentation fault if any of the sources
                                        # registered using pc.cell() no longer exist, so
                                        # we keep a reference to all sources in the
                                        # global gid_sources list. It would be nicer to
                                        # be able to unregister a gid and have a __del__
                                        # method in ID, but this will do for now.


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
            if ((i%3)==0):
                bothtrans.append(interneurons[i]) 
            else:
                bothtrans.append(excitatory[i])
        return bothtrans#(excitatory, interneurons)      
        
    def read_local_swc(self):
        h=self.h    
        NCELL=self.NCELL
        SIZE=self.SIZE
        RANK=self.RANK
        #from neuron import h
        pc=h.ParallelContext()
        
        morphs=[]
        cells1=[]
        self.initialize_hoc()
        swclist=glob.glob('*.swc')
        itergids = iter( i for i in range(RANK, len(swclist), SIZE) )

        #for swcf, i in enumerate(swclist):
        for i in itergids:
            #morphology = swc.read_swc(swcf)
            cell = h.mkcell(swclist[i])
                #self.generate_morphology(cell, d[i][3])
            self.generate_morphology(cell,swclist[i])
            self.load_cell_parameters(cell, fit_ids[utils.cells_data[i]['type']])
            
            #cell1=self.load_cell_parameters()       
            cells1.append(cell1)
            print type(cells1)
            print type(cell1)
            morphology.root
            morphs.append(morphology)
        return morphs,swclist,cells1


    def make_cells(self,polarity):
        h=self.h    
        NCELL=self.NCELL
        SIZE=self.SIZE
        RANK=self.RANK
        #from neuron import h
        pc=h.ParallelContext()
        d = { x: y for x,y in enumerate(polarity)}         
        itergids = iter( i for i in range(RANK, NCELL, SIZE) )        
        #TODO keep rank0 free of cells, such that all the memory associated with that CPU is free for graph theory related objects.
        itergids = iter( i for i in range(RANK+1, NCELL, SIZE-1) )        
        fit_ids = self.description.data['fit_ids'][0] #excitatory        
        for i in itergids:
            cell = h.mkcell(d[i][3])
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
        print np.sum(pol)
        os.chdir(os.getcwd() + '/../')               
        self.h.define_shape()        
        self.h('forall{ for(x,0){ insert xtra }}')
        self.h('forall{ for(x,0){ insert extracellular}}')    
        self.h('xopen("interpxyz.hoc")')
        self.h('grindaway()')    
         

 # no .clear() command
        
    def htype (obj): st=obj.hname(); sv=st.split('['); return sv[0]
    def secname (obj): obj.push(); print self.h.secname() ; self.h.pop_section()
    def psection (obj): obj.push(); print self.h.psection() ; self.h.pop_section()
    
    # still need to generate a full allsecs
    def mkallsecs():
        """ mkallsecs - make the global allsecs variable, containing
        all the NEURON sections.
        """
        #global allsecs
        allsecs=self.h.SectionList()
        return allsecs

    def matrix_reduce(self):       
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        COMM.Barrier()
        self.my_icm = np.zeros_like(self.icm)
        COMM.Reduce([self.icm, MPI.DOUBLE], [self.my_icm, MPI.DOUBLE], op=MPI.SUM,
                    root=0)
        self.my_ecm = np.zeros_like(self.ecm)
        COMM.Reduce([self.ecm, MPI.DOUBLE], [self.my_ecm, MPI.DOUBLE], op=MPI.SUM,
                    root=0)

    def graph_reduce(self):       
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        COMM.Barrier()
        self.my_icg = self.icg
        COMM.Reduce([self.icm, MPI.DOUBLE], [self.my_icm, MPI.DOUBLE], op=MPI.SUM,
                    root=0)
        self.my_ecg = self.ecg
        COMM.Reduce([self.ecm, MPI.DOUBLE], [self.my_ecm, MPI.DOUBLE], op=MPI.SUM,
                    root=0)
        

    def vec_reduce(self,idvec,tvec):       
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        COMM.Barrier()
        my_tvec = np.zeros_like(idvec)
        COMM.Reduce([tvec, MPI.DOUBLE], [self.my_tvec, MPI.DOUBLE], op=MPI.SUM,
                    root=0)
        my_idvec = np.zeros_like(tvec)
        COMM.Reduce([idvec, MPI.DOUBLE], [self.my_idvec, MPI.DOUBLE], op=MPI.SUM,
                    root=0
        )
        return my_idvec, my_tvec


    def prun(self,tstop):
        h=self.h    
        pc=h.ParallelContext()
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        checkpoint_interval = 50000.

        #This code is from:
        #http://senselab.med.yale.edu/ModelDB/ShowModel.asp?model=151681
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


    def nestedpre_test(self,j):
        #from neuron import h
        h=self.h   
        pc=h.ParallelContext()
        shiplist=[]
        #coordict={}
        #pc=h.pc
        h('objref coords') 
        h('coords = new Vector(3)')
        #self.celldict.items()[0]
        if j in self.celldict.keys():
            seglist= iter( (seg, sec, self.celldict[j]) for sec in self.celldict[j].spk_trig_ls for seg in sec )     
            for (seg,sec, cellc) in seglist:
                sec.push()
                #sec.nseg
                #for (seg,sec, cellc) in seglist:
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
                #print secnames, seg.x                       
                #shiplist.append((secnames,seg.x))
                coordict['secnames'] = str(secnames)
                shiplist.append(coordict)
                self.h.pop_section()
               
        return shiplist
        
    def alloc_synapse_ff(self,r,post_syn,cellind,k,gidn,i):
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        from neuron import h
        pc=h.ParallelContext()
        h=self.h  
        polarity = 0        
        polarity=int(h.Cell[int(cellind)].polarity)
        if polarity==1:
            self.ecg.add_edge(i,gidn,weight=r/0.4)
            #Add other edge attributes like secnames etc.
        else:
            self.icg.add_edge(i,gidn,weight=r/0.4)
            #Add other edge attributes like secnames etc.

                  
        h('objref syn_')   
        h(post_syn)
        print post_syn
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
        if r < 5:
            polarity = 0        
            polarity=int(h.Cell[int(cellind)].polarity)
            h('objref syn_')        
            if int(polarity) == int(0):
                post_syn = secnames + ' ' + 'syn_ = new GABAa(' + str(seg.x) + ')'
                self.icm[i][gidn] = self.icm[i][gidn] + 1
                self.icg.add_edge(i,gidn)
                self.icg.add_edge(i,gidn,weight=r/0.4)
                self.icg[i][gidn]['post_loc']=secnames
                self.icg[i][gidn]['pre_loc']=str(sec.name())

                assert np.sum(self.icm)!=0
        
            else:
                #TODO modify mod file for exp2syn such that cid exists.
                #Also because exp2syn is an inbuilt mechanism need to refactor explicitly such that custom file 
                #myexp2syn is used instead.
                #post_syn = secnames + ' ' + 'syn_ = new exp2syn(' + str(seg.x) + ')'
                #post_syn='syn_ = self.h.Exp2Syn('+str(seg.x)+',sec='+secnames+')'                        
                #syn_.e = connection["erev"]
                if (k['gid']%2==0):
                    post_syn = secnames + ' ' + 'syn_ = new ExpSid(' + str(seg.x) + ')'                       
                    self.ecm[i][gidn] = self.ecm[i][gidn] + 1
                    
                    self.ecg.add_edge(i,gidn,weight=r/0.4)
                    self.ecg[i][gidn]['post_loc']=secnames
                    self.ecg[i][gidn]['pre_loc']=str(sec.name())

                    self.seclists.append(secnames)
                    assert np.sum(self.ecm)!=0
                else:
                    post_syn = secnames + ' ' + 'syn_ = new NMDA(' + str(seg.x) + ')'                       
                    self.ecm[i][gidn] = self.ecm[i][gidn] + 1
                    self.ecg.add_edge(i,gidn,weight=r/0.4)
                    self.ecg[i][gidn]['post_loc']=secnames
                    self.ecg[i][gidn]['pre_loc']=str(sec.name())
                    self.seclists.append(secnames)
                    assert np.sum(self.ecm)!=0
        
            h(post_syn)
            self.synapse_list.append((r,post_syn,cellind,k,gidn,i))
            print post_syn
            h('print syn_')
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

        
       
    def nestedpost_test(self,data):
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
        from segment_distance import dist3D_segment_to_segment
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
        for s in data: 
            k={} #The only point of this redundantvariable switching is to force the dictionary k to be redclared 
            k=s #such that it is not prevented from updating
            itercell= ( (i,t) for i,t in self.celldict.iteritems() if i in self.celldict.keys() if int(t.gid1) != int(k['gid']) )       
            for i,t in itercell :                          
                iterseg=iter( (seg,sec) for sec in t.spk_rx_ls for seg in sec)               
                for (seg,sec) in iterseg:
                    secold=sec.name()
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
        #h('coordsx') and coordsx are not tautolous they are
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



    def wirecells_test(self):
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
        coordictlist=None
        #Iterate over all CPU ranks, iterate through all GIDs (global 
        #identifiers, stored in the python dictionary).
        if self.readin!=1:    
            for s in xrange(0, SIZE):
                celliter= iter( (i, j) for i,j in self.celldict.iteritems() )  
                for (i,j) in celliter:  
                    cell1=pc.gid2cell(i)
                    coordictlist=self.nestedpre_test(i)
                data = COMM.bcast(coordictlist, root=s)  # ie root = rank
                if len(data) != 0:
                    self.nestedpost_test(data)
            print('finished wiring of connectivity\n')
            fname='synapse_list'+str(RANK)+'.p'
            assert len(self.synapse_list)!=0
            with open(fname, 'wb') as handle:
                pickle.dump(self.synapse_list, handle)
        else:
            fname='synapse_list'+str(RANK)+'.p'
            with open(fname, 'rb') as handle:
                #print pickle.load(handle)
                self.synapse_list=pickle.load(handle)
                for s in self.synapse_list:
                    print s
                for (r,post_syn,cellind,k,gidn,i) in self.synapse_list:
                    self.alloc_synapse_ff(r,post_syn,cellind,k,gidn,i)
        #return (self.nclist, self.ecm, self.icm)

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
                print 'this is evidently not printing, some of these objects have been destroyed!'
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
        #print('Or copy all files in force/ to webserver and load force/force.html')



    def generate_morphology(self, cell, morph_filename):
        h = self.h
        
        swc = self.h.Import3d_SWC_read()
        swc.input(morph_filename)
        imprt = self.h.Import3d_GUI(swc, 0)
        imprt.instantiate(cell)
        
        for seg in cell.soma[0]:
            seg.area()

        for sec in cell.all():
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

    def connect_cells(self):
        self.synlist = []
        self.nclist = []
        connections = self.description.data["biophys"][0]["connections"]

        for connection in connections:
            for target in connection["targets"]:
                source_cell = self.cells[connection["source"]]
                target_cell = self.cells[target]

                syn = self.h.Exp2Syn(0.5, sec=target_cell.dend[0])
                syn.e = connection["erev"]
                source_section = source_cell.soma[0]
                nc = self.h.NetCon(source_section(0.5)._ref_v, syn, sec=source_section)
                nc.weight[0] = connection["weight"]
                nc.threshold = -20
                nc.delay = 2.0

                self.synlist.append(syn)
                self.q.append(nc)


    def connect_ring(self):
        self.synlist = []
        self.nclist = []
        connections = self.description.data["biophys"][0]["connections"]

        #for connection in connections:
        #    for target
            #in connection["targets"]:
        for i,discard in enumerate(self.cells):
            for j, discard in enumerate(self.cells):
                if i!=j:
                    print type(i), type(j)
                    
                    source_cell = self.cells[i]
                    target_cell = self.cells[j]
                    syn = self.h.Exp2Syn(0.5, sec=target_cell.dend[0])
                    self.nclist[0].syn.Section
                    #syn.e = connection["erev"]
                    #syn.e = -0.6
                    source_section = source_cell.soma[0]
                    nc = self.h.NetCon(source_section(0.5)._ref_v, syn, sec=source_section)
                    nc.weight[0]=0.005                   
                    #nc.weight[0] = connection["weight"]
                    nc.threshold = -20
                    nc.delay = 2.0        
                    self.synlist.append(syn)
                    self.nclist.append(nc)


    #from neuromac.segment_distance import dist3D_segment_to_segment
        

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


    '''
    def proc spikeout(self):   
        self.COMM.barrier()# // wait for all hosts to get to this point
        if (self.COMM.rank==0):
            print("\ntime\t cell\n") #// print header once
    
        for rank in xrange(0,self.COMM.size): #{ // host 0 first, then 1, 2, etc.
            if (rank==self.COMM.rank):
                print(fno,"%s_spt.dat", fstem)
                fo = new File(fno)  
                fo.aopen()
                for i in len(tvec)#.size-1:   
                    printf ("%g\t %d\n", tvec.x[i], idvec.x[i])
                    fo.printf("%g\t %d\n", tvec.x[i], idvec.x[i])
                fo.close()
        self.COMM.barrier()
    '''                   
  
        
