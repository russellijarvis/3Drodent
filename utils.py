from allensdk.model.biophys_sim.neuron.hoc_utils import HocUtils
import logging
import glob
from mpi4py import MPI
import btmorph
import numpy

import numpy as np
class Utils(HocUtils):

    _log = logging.getLogger(__name__)
    
    def __init__(self, description):
      
        super(Utils, self).__init__(description)
        self.stim = None
        self.stim_curr = None
        self.sampling_rate = None
        self.cells = []
        self.gidlist=[]
        self.NCELL=0
        self.celldict={}
        self.COMM = MPI.COMM_WORLD
        self.SIZE = self.COMM.Get_size()
        self.RANK = self.COMM.Get_rank()
        self.allsecs=None #global list containing all NEURON sections, initialized via mkallsecs
        self.coordict=None
        self.celldict={}
        self.cellmorphdict={}
        self.nclist = []

        icm = np.zeros((self.NCELL, self.NCELL))
        ecm = np.zeros((self.NCELL, self.NCELL))
        if self.RANK==0:        
            my_icm = np.zeros((self.NCELL, self.NCELL))
            my_ecm = np.zeros((self.NCELL, self.NCELL))

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


    
    def prep_list(self):                
        import pickle
        allrows = pickle.load(open('allrows.p', 'rb'))
        allrows.remove(allrows[0])#The first list element are the column titles. 
        allrows2 = [i for i in allrows if int(len(i))>9 ]
        return allrows2        



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
    
    def gcs(self,NCELL):
        """Instantiate NEURON cell Objects in the Python variable space, such that cell
        that all cells have unique identifiers."""
        NCELL=self.NCELL
        SIZE=self.SIZE
        RANK=self.RANK
        h=self.h    
        pc=h.ParallelContext()        
        h('objref pc, py, nc, cells')
        h('pc = new ParallelContext()')
        h('py = new PythonObject()')
        swcdict={}
        NFILE = 3175
        fit_ids = self.description.data['fit_ids'][0] #excitatory        
        self.cells_data = self.description.data['biophys'][0]['cells']
        info_swc=self.prep_list()
        d = { x: y for x,y in enumerate(info_swc)}
        import os
        os.chdir(os.getcwd() + '/main')   
        itergids = iter( i for i in range(RANK, NCELL, SIZE) )
        for i in itergids:
            cell = h.mkcell(d[i][3])
            cell.geom_nseg()
            cell.gid1=i #itergids.next()
            #excitatory cell.
            if 'pyramid' in d[i]:            
                self.load_cell_parameters(cell, fit_ids[self.cells_data[0]['type']])
                cell.polarity=1
            else:            
            #inhibitory cell.
                self.load_cell_parameters(cell, fit_ids[self.cells_data[2]['type']])
                cell.polarity=0
            h('Cell[0].soma[0] nc =  new NetCon(&v(0.5), nil)')                        
            pc.set_gid2node(i,RANK)
            h('pc.cell('+str(i)+', nc)')
            cell1=pc.gid2cell(i)
            self.celldict[i]=cell
            self.cells.append(cell)
        len(h.List('NetCon'))                        
        pol=[ a.polarity for a in self.cells ]       
        import numpy as np
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

    def matrix_reduce(self,ecm,icm):       
        import numpy as np
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        icm = np.zeros((NCELL, NCELL))
        ecm = np.zeros((NCELL, NCELL))
        COMM.Barrier()
        my_icm = np.zeros_like(icm)
        COMM.Reduce([icm, MPI.DOUBLE], [my_icm, MPI.DOUBLE], op=MPI.SUM,
                    root=0)
        my_ecm = np.zeros_like(ecm)
        COMM.Reduce([ecm, MPI.DOUBLE], [my_ecm, MPI.DOUBLE], op=MPI.SUM,
                    root=0)
        return ecm, icm
        
       
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
    
        # if h.t%2==0: The problem is h.t is float multiple not integer multiple
    
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
        """
        
        This function searches for pre synaptic sources in the parallel wiring algorim
        This algorithm includes the soma as proxy for other regions which can exocytosis
        TODO: this algorithm is not currently limited to leaf nodes (such as dendrite sp
        """
        from neuron import h
        pc=h.ParallelContext()
        h=self.h   
        import numpy as np
        coordictlist=[]
        coordict={}
        pc=h.pc
        h('objref coords') 
        h('coords = new Vector(3)')
        if j in self.celldict.keys():
            seglist= iter( (seg, sec, self.celldict[j]) for sec in self.celldict[j].spk_trig_ls for seg in sec )              
            for (seg,sec, cellc) in seglist:
                    get_cox = str('coords.x[0]=x_xtra('
                                  + str(seg.x) + ')')
                    h(get_cox)                   
                    get_coy = str('coords.x[1]=y_xtra('
                                  + str(seg.x) + ')')
                    h(get_coy)
                    get_coz = str('coords.x[2]=z_xtra('
                                  + str(seg.x) + ')')
                    h(get_coz)
                    coordict['hostfrom']=pc.id()
                    coordict['coords'] = np.array(h.coords.to_python(),
                                              dtype=np.float64)
                    coordict['gid']= int(j)
                    coordict['seg']= seg.x                    
                    secnames = sec.name()  # h.secnames                            
                    coordict['secnames'] = str(secnames)
                    h('print x_xtra('+ str(seg.x) +')')
                    coordictlist.append(coordict)               
        print len(coordictlist)                                   
        return coordictlist
             

        '''
        Long term goal is to make tracenet a class method 
        ecm and icm should be class attributes.
    class NetStat(self):
        #They advantage of creating a seperate object for stats of abstract facts about the model
        #Versus the biological model is because
        __init__(self):
        self.ecm
        self.icm
        '''
       
    def tracenet(self):
        ncsize=len(self.h.NetCon)
        import numpy as np
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        icm = np.zeros((NCELL, NCELL))
        ecm = np.zeros((NCELL, NCELL))
        COMM.Barrier()
        my_icm = np.zeros_like(icm)
        COMM.Reduce([icm, MPI.DOUBLE], [my_icm, MPI.DOUBLE], op=MPI.SUM,
                    root=0)
        my_ecm = np.zeros_like(ecm)
        COMM.Reduce([ecm, MPI.DOUBLE], [my_ecm, MPI.DOUBLE], op=MPI.SUM,
                    root=0)
        lsoftup=[]
        #return ecm, icm
        #make a list of tuples where each list element contains (srcid,tgtid,srcpop,tgtpop)
        #for s in xrange(0,SIZE):            
        for i in xrange(0,ncsize-1):
            #srcs.append(int(self.h.NetCon[i].srcgid()))
            #tgts.append(int(self.h.NetCon[i].postcell().gid1))
            srcind=int(self.h.NetCon[i].srcgid())
            tgtind=int(self.h.NetCon[i].postcell().gid1)
            print int(utils.h.NetCon[i].srcgid()),int(utils.h.NetCon[i].postcell().gid1),utils.celldict[srcind],utils.celldict[tgtind]
                #print strlist[tgtind]==dic[tgtind], ' sanity check '
                #add to list of tuples, netcon src, netcon tgt, src index, target index.
            lsoftup.append((int(utils.h.NetCon[i].srcgid()),int(utils.h.NetCon[i].postcell().gid1),utils.celldict[srcind],utils.celldict[tgtind]))
        return lsoftup
        #The broadcasting and gathering should happen on a different host.
        #Actually this should all be reduced to rank0
        #    data = COMM.bcast(lsoftup, root=s)  # ie root = rank
        #return data        
  
       
 
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
        from segment_distance import dist3D_segment_to_segment
        import numpy as np
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        icm = np.zeros((NCELL, NCELL))
        ecm = np.zeros((NCELL, NCELL))
        #self.nclist
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
        iterdata=iter( (k,i,t) for k in data for i,t in self.celldict.iteritems() if i in self.celldict.keys() if int(t.gid1) != int(k['gid']))
                 
        for k,i,t in iterdata :                          
            iterseg=iter( (seg,sec) for sec in t.spk_rx_ls for seg in sec)                    
            for seg,sec in iterseg:
                print seg.x, sec.name(), k['secnames']
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
                if r < 1:  
    
                    print r,# 'this is not hopefuly wiring everything to everything'
                    polarity = 0        
                    polarity=int(h.Cell[int(cellind)].polarity)
                    print seg.x, k['seg'], k['secnames'], sec.name(), RANK, k['hostfrom'], k['gid'], int(h.Cell[int(cellind)].gid1)
                    
                    print polarity
                    h('objref syn_')        
                    if int(polarity) == int(0):
                        post_syn = secnames + ' ' + 'syn_ = new GABAa(' + str(seg.x) + ')'
                        icm[i][gidn] = icm[i][gidn] + 1
                    else:
    
                        post_syn = secnames + ' ' + 'syn_ = new AMPA(' + str(seg.x) + ')'
                        ecm[i][gidn] = ecm[i][gidn] + 1
    
                    h(post_syn)
                    print post_syn
                    h('print syn_')
                    syn_=h.syn_
                    h.syn_.cid=i
                    h.Cell[cellind].ampalist.append(h.syn_)
                    h.Cell[cellind].div.append(k['gid'])
                    h.Cell[cellind].gvpre.append(k['gid'])
                    nc=pc.gid_connect(k['gid'],syn_)                                        
                    nc.delay=1+r/0.4
                    nc.weight[0]=r/0.4    
                    self.nclist.append(nc)
            h('uninsert xtra')                          
        return self.nclist, ecm, icm
    
    def wirecells(self):
        """This function constitutes the outermost loop of the parallel wiring algor
        The function returns two adjacency matrices. One matrix whose elements are excitatory connections and another matrix of inhibitory connections"""
        from segment_distance import dist3D_segment_to_segment
        import numpy as np
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        icm = np.zeros((NCELL, NCELL))
        ecm = np.zeros((NCELL, NCELL))
        #self.nclist
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
        for s in xrange(0, SIZE):
            celliter= iter( (i, j) for i,j in self.celldict.iteritems() )  
            for (i,j) in celliter:  
                #print i==j           
                cell1=pc.gid2cell(i)
                coordictlist=self.nestedpre(i)
                print 'entered parallel wiring now', s, i, j
            data = COMM.bcast(coordictlist, root=s)  # ie root = rank
            
            if len(data) != 0:
                self.nclist, ecm, icm = self.nestedpost(data)

   
        data=None                        
        ecm,icm = self.matrix_reduce(ecm,icm)
        return (self.nclist, ecm, icm)


    def generate_morphology(self, cell, morph_filename):
        h = self.h
        
        swc = self.h.Import3d_SWC_read()
        swc.input(morph_filename)
        imprt = self.h.Import3d_GUI(swc, 0)
        imprt.instantiate(cell)
        
        for seg in cell.soma[0]:
            seg.area()

        for sec in cell.all:
            sec.nseg = 1 + 2 * int(sec.L / 40)
        
        #cell.simplify_axon()
        #for sec in cell.axonal:
        #    sec.L = 30
        #    sec.diam = 1
        #    sec.nseg = 1 + 2 * int(sec.L / 40)
        #cell.axon[0].connect(cell.soma[0], 0.5, 0)
        #cell.axon[1].connect(cell.axon[0], 1, 0)
        h.define_shape()
    
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
                if p["mechanism"] != "":
                    sec.insert(p["mechanism"])
                    h('print psection()')
                setattr(sec, p["name"], p["value"])
        
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
        vec = { "v": [],
                "t": self.h.Vector() }
    
        for i, cell in enumerate(self.cells):
            vec["v"].append(self.h.Vector())
            vec["v"][i].record(cell.soma[0](0.5)._ref_v)
        vec["t"].record(self.h._ref_t)
    
        return vec
        
        
    def wirecells_s(self):
        '''wire cells on the same hosts'''
        import numpy as np
        NCELL=self.NCELL
        SIZE=self.SIZE
        COMM = self.COMM
        RANK=self.RANK
        icm = np.zeros((NCELL, NCELL))
        ecm = np.zeros((NCELL, NCELL))
        #self.nclist
        from neuron import h
        pc=h.ParallelContext()
        h=self.h
        celliter= iter( (j, i) for j,l in self.celldict.iteritems() for i,t in self.celldict.iteritems() if i!=j )  
        for (j,i) in celliter:  
            print i==j           
            cell1=pc.gid2cell(i)            
            coordictlist=self.inner1(j)
            for coordict in coordictlist:    
                seglist= iter( (seg, sec, self.celldict[j]) for sec in self.celldict[j].spk_rx_ls for seg in sec )                          
                for (seg, sec, cellc) in seglist:
                    secnames = sec.name()
                    cellind = int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.
                    h('objref coords2') 
                    h('coords2 = new Vector(3)')
                    h(str('coords2.x[2]=') + str('z_xtra(')
                      + str(seg.x) + ')')
                    h(str('coords2.x[1]=') + str('y_xtra(')
                      + str(seg.x) + ')')
                    h(str('coords2.x[0]=') + str('x_xtra(')
                    + str(seg.x) + ')')
                    h('coordsx=0.0')
                    h.coordsx = coordict['coords'][0]
                    h('coordsy=0.0')
                    h.coordsy = coordict['coords'][1]
                    h('coordsz=0.0')
                    h.coordsz = coordict['coords'][2]
      
                    coordsx = float(coordict['coords'][0])
                    coordsy = float(coordict['coords'][1])
                    coordsz = float(coordict['coords'][2])
                    r = 0.
                    import math
                    r=math.sqrt((h.coords2.x[0] - coordsx)**2+(h.coords2.x[1] - coordsy)**2+(h.coords2.x[2] - coordsz)**2)
                    gidn=coordict['gid']    
                    r = float(r)
                    if r < 1:  
                        print r,# 'this is not hopefuly wiring everything to everything'
                        gidcompare = ''
                        polarity = 0
                        #cellind is a cell index, that is relative to the host. So the identifier repeats on different hosts.
                        #gidn is a global identifier. These numbers are not repeated on different hosts.
                        polarity=int(h.Cell[int(cellind)].polarity)
                        #print seg.x, coordict['seg'], coordict['secnames'], sec.name(), RANK, coordict['hostfrom'], coordict['gid'], int(h.Cell[int(cellind)].gid1)                        
                        #print polarity
                        h('objref syn_')        
                        if int(polarity) == int(0):
                            post_syn = secnames + ' ' + 'syn_ = new GABAa(' + str(seg.x) + ')'
                            icm[i][gidn] = icm[i][gidn] + 1
                        else:

                            post_syn = secnames + ' ' + 'syn_ = new AMPA(' + str(seg.x) + ')'
                            ecm[i][gidn] = ecm[i][gidn] + 1

                        h(post_syn)
                        h('print syn_')
                        syn_=h.syn_
                        h.syn_.cid=i
                        h.Cell[cellind].ampalist.append(h.syn_)
                        h.Cell[cellind].div.append(coordict['gid'])
                        h.Cell[cellind].gvpre.append(coordict['gid'])
                        h('objref nc')            
                        ls=str(coordict['secnames'])+' nc =  new NetCon(&v('+str(coordict['seg'])+'),'+str(sec.name())+')'
                        print str(sec.name()), coordict['secnames'] , coordict['seg']                       
                        h(ls)
                        nc=h.nc                                                    
                        nc.delay=1+r/0.4
                        nc.weight[0]=r/0.4   
                        self.nclist.append(nc)

        ecm,icm = self.matrix_reduce(ecm,icm)
        return (self.nclist, ecm, icm)

                   

        
