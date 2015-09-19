from allensdk.model.biophys_sim.neuron.hoc_utils import HocUtils
import logging
import glob
from mpi4py import MPI

# initialize the MPI interface
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

#pc = h.ParallelContext()
#h('objref pc')
#h.pc = pc
#s = "mpi4py thinks I am %d of %d,\
# NEURON thinks I am %d of %d\n"
#cw = MPI.COMM_WORLD
#print s % (cw.rank, cw.size, pc.id(), pc.nhost())
#h('time_start=pc.time()')
#PROJECTROOT=os.getcwd()

class Utils(HocUtils):
    _log = logging.getLogger(__name__)
    
    def __init__(self, description):
        super(Utils, self).__init__(description)
        self.stim = None
        self.stim_curr = None
        self.sampling_rate = None
        self.cells = []
        self.RANK=0
        self.NCELL=0
        self.SIZE=0
    
# I do not know how to refer to relative paths in Python,
# the below emulates a call to a relative path.


    
    def prep_list(self):                
        import pickle
        allrows = pickle.load(open('allrows.p', 'rb'))
        allrows.remove(allrows[0])#The first list element are the column titles. 
        allrows2 = [i for i in allrows if int(len(i))>9 ]
        return allrows2        







    def gcs(self,NCELL):
        swcdict={}
        NFILE = 3175
        fit_ids = self.description.data['fit_ids'][0] #excitatory        
        self.cells_data = self.description.data['biophys'][0]['cells']
        info_swc=self.prep_list()
        d = { x: y for x,y in enumerate(info_swc)}
        #d = {a[3]: a for a in info_swc}
        import os
        os.chdir(os.getcwd() + '/main')        
        #swclist=glob.glob('*.swc')
        # Filter out the bad morphology, using a list comprehension.
        #swclist.remove("Scnn1a-Tg3-Cre_Ai14_IVSCC_-177300.01.02.01_473845048_m.swc")            
        #swclist.remove("466664172.swc") 
           
        gids = [ i for i in range(RANK, NCELL, SIZE) ]
        itergids=iter(gids)
        iterd=iter(d)
        
        for i in itergids:
        #for i,j in enumerate(d):
            cell = self.h.cell()
            #print len(itergids)
            #print len(iterd)
            cell.gid=i #itergids.next()

            print cell.gid    
            self.generate_morphology(cell, d[i][3])#iterd.next())#iterswc.next())
            if 'pyramid' in d[i]:            
                self.load_cell_parameters(cell, fit_ids[self.cells_data[0]['type']])
                cell.polarity=1
            else:            
            #inhibitory type stained cell.
                self.load_cell_parameters(cell, fit_ids[self.cells_data[2]['type']])
                cell.polarity=0

            
            self.cells.append(cell)
        from neuron import h
        pol=[ a.polarity for a in self.cells ]       
        import numpy as np
        print np.sum(pol)
        os.chdir(os.getcwd() + '/../')               
        self.h('forall{ for(x,0){ insert xtra }}')
        self.h('forall{ for(x,0){ insert extracellular}}')    
        self.h('xopen("interpxyz.hoc")')
        self.h('grindaway()')    
        self.h('xopen("seclists.hoc")')

    
    def wirecells3(self):
    #def wirecells(RANK,NCELL,SIZE,h,icm,ecm):
        from segment_distance import dist3D_segment_to_segment
        import numpy as np
        icm = np.zeros((self.NCELL, self.NCELL))
        ecm = np.zeros((self.NCELL, self.NCELL))
        h=self.h    
        pc=h.ParallelContext()
        #global s,j,i,test,test2,r
        self.s=0
        self.j=0
        self.i=0
        self.gidcompare = ''
     
        secnames = ''# sec.name()
        cellind =0 #int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.
    
        polarity = 0
        
        self.h.py.j=0
    
        for s in xrange(0, SIZE):#Was full SIZE, not SIZE-1
            if RANK == s:
    
                for j in xrange(0,NCELL):
                    self.py.j=int(j)
                #print s, j
                    coordlist = [[0 for x in xrange(0, 2)] for x in xrange(0,
                             2)]
    
                    test = 0
                    #h('print py.j, py.test')
                    test=int(pc.gid_exists(j))
                    print 'test ', test
                    #h('py.test=py.int(pc.gid_exists(py.int(py.j)))')
                    if test != 0:
                        h.source=pc.gid2cell(j)
                        print h.source
                        #h('source=pc.gid2cell(py.int(py.j))')
                        #h('print source')
                        h.source=h.pc.gid2cell(int(j))
                        h('print source')
    
                        for sec in h.source.spk_trig_ls:
                            for seg in sec:
                                #print sec, seg
                                get_cox = str('coords.x[0]=x_xtra('
                                              + str(seg.x) + ')')
                                h(get_cox)
    
    
                                get_coy = str('coords.x[1]=y_xtra('
                                              + str(seg.x) + ')')
                                h(get_coy)
                                get_coz = str('coords.x[2]=z_xtra('
                                              + str(seg.x) + ')')
                                h(get_coz)
                                h.coords.x[3] = int(j)
                                h.coords.x[4] = seg.x # ie the gidvec
                                #h('coords.x[4]=py.seg.x')  # ie the gidvec.
    
                                coords = np.array(h.coords.to_python(),
                                                  dtype=np.float64)
    
                                coordlist[0][1] = coords
    
    
                                secnames = sec.name()  # h.secnames
                                coordlist[1][0] = str(secnames)
    
    
                    if test == 0:
                        coordlist = None
                else:
                    coordlist = None
                data = COMM.bcast(coordlist, root=s)  # ie root = rank
    
                if data != None: 
                    #pdb.set_trace()
     
                    #h.py.i=0
                    for i in xrange(0, NCELL-1):
                        h.py.i=i
                        gidn = int(data[0][1][3])
                        if i != int(gidn):  # if the gids are not the same
                            test2 = 0
                            test2=pc.gid_exists(i)
                            #h('py.test2=py.int(pc.gid_exists(py.int(py.i)))')
                            print test2
                            if test2 != 0:
                                h.target=pc.gid2cell(i)
                                print h.target
                                h('target=pc.gid2cell(py.int(py.i))')
                                h('print target')
                                for sec in h.target.spk_rx_ls:
                                    for seg in sec:
                                        h(str('coords2.x[2]=') + str('z_xtra(')
                                          + str(seg.x) + ')')
                                        h(str('coords2.x[1]=') + str('y_xtra(')
                                          + str(seg.x) + ')')
                                        h(str('coords2.x[0]=') + str('x_xtra(')
                                        + str(seg.x) + ')')
    
                                        h('coordsx=0.0')
                                        h.coordsx = data[0][1][0]
                                        h('coordsy=0.0')
                                        h.coordsy = data[0][1][1]
                                        h('coordsz=0.0')
                                        h.coordsz = data[0][1][2]
                          
                                        coordsx = float(data[0][1][0])
                                        coordsy = float(data[0][1][1])
                                        coordsz = float(data[0][1][2])
                                        r = 0.
                                        r=math.sqrt((h.coords2.x[0] - coordsx)**2+(h.coords2.x[1] - coordsy)**2+(h.coords2.x[2] - coordsz)**2)
                                        if front.parent == None or o_front.parent == None:
                                            D = np.sqrt(np.sum((front.xyz-o_front.xyz)**2))
                                        else:
                                            D = dist3D_segment_to_segment (front.xyz,front.parent.xyz,o_front.parent.xyz,o_front.xyz)
    
                                        #h('py.r = sqrt((coords2.x[0] - coordsx)^2 + (coords2.x[1] - coordsy)^2 + (coords2.x[2] - coordsz)^2)')
                                        r = float(r)
                                        if r < 10:  
                                            print r, 'this is not hopefuly wiring everything to everything'
                                            gidcompare = ''
    
                                            secnames = sec.name()
                                            cellind = int(secnames[secnames.find('Cell[') + 5:secnames.find('].')])  # This is the index of the post synaptic cell.
    
                                            polarity = 0
                                       
    
                                            #cellind is a cell index, that is relative to the host. So the identifier repeats on different hosts.
                                            #gidn is a global identifier. These numbers are not repeated on different hosts.
                                            polarity=int(h.Cell[int(cellind)].polarity)
                                            
                                            print polarity
    
                                            if int(polarity) == int(1):
                                                post_syn = secnames + ' ' + 'syn_ = new GABAa(' + str(seg.x) + ')'
                                                icm[i][gidn] = icm[i][gidn] + 1
                                            else:
    
                                                post_syn = secnames + ' ' + 'syn_ = new AMPA(' + str(seg.x) + ')'
                                                ecm[i][gidn] = ecm[i][gidn] + 1
    
                                            h('gidn=0')
                                            h.gidn = int(data[0][1][3])
    
                                            h(post_syn)
                                            print post_syn
                                            h('print syn_')
                                            #h.syn_.cid=i
                                            h('syn_.cid=py.int(py.i)')  # This makes the gid a property of the synapse.
                                            h.Cell[cellind].ampalist.append(h.syn_)
    
    
    
                                            h('synlist.append(syn_)')
                                            h('gidn=0')
                                            h.gidn = int(data[0][1][3])
                                            h.Cell[cellind].div.append(h.gidn)
                                            h.Cell[cellind].gvpre.append(h.gidn)
    
    
                                            h('nc = pc.gid_connect(gidn, syn_)')
    
                 
    
                                            h('print nc," ", nc.srcgid()," ", gidn')
                                            print 'the gids of the connected cells are: ', i, ' ', gidn, '\n'
    
    
                                            h.py.r=r
                                            h('nc.delay = 1+((py.r)/delay)')
                                            h('nc.weight = py.r/iw')
                                            h('nclist.append(nc)')
    
    
        COMM.Barrier()
    
        return (h.nclist, ecm, icm)

       

 
    
    def mb(RANK, NCELL, SIZE, allrows2, gidvec, h,s1):
        #make soff and boff optional parameters in python function.
        #pdb.set_trace()
        #h.load_file("stdlib.hoc")
        #global s1
        #global i, cnt, cnti
        #h('chdir(workingdir)')
        #os.chdir(os.getcwd() + '/main')
        h('py = new PythonObject()')
        h('cnt=0')
        ie0 = np.zeros((NCELL, 2))
        ie1 = np.zeros((NCELL, 2))
    
        i=0
        cnti = 0
        cnt=0
        s1=''
        storename=''
        h('cnt=pc.id')
           
        for i in range(RANK, NCELL,SIZE):#NCELL-1, SIZE):  # 20 was int(len
            h.py.i=i
            s1 = allrows2[i]
            h.py.s1=s1
    
            storename = str(s1[3])  # //simply being inside a loop, may be the main problem
            if re.search('.swc', storename):
                h.cell = h.mkcell(storename)
    
    
                h('cell.geom_nseg()')
                #h.cell.geom_nseg()
                h('cell.gid1=py.i')
                h('cell.gvpre.printf')
                h('cell.gvpost.printf')
                
    
    
    
                h.cell.nametype=str(s1[5])
                #h('cell.nametype=py.str(py.s1[5])')
                h.cell.num_type=int(s1[6])
                #h('cell.num_type=py.int(py.s1[6])')
                h('cell.population=py.str(py.s1[7])')
                #h('cell.reponame=py.str(py.storename)')
                h.cell.reponame=str(storename)
                h('cell.div.resize(py.int(py.NCELL))')
                h('cell.conv.resize(py.int(py.NCELL))')
                #h('cell.nametype=py.str(py.s1[5])')
    
                h('if(strcmp("pyramid",py.s1[5])==0){pyr_list.append(cell)}')
    
                h('if(strcmp(cell.population,"neocortex")==0){ if(strcmp("pyramid",py.s[5])==0){ cell.pyr() }}')
    
    
    
                h('if(strcmp(cell.population,"hippocampus")==0){ if(strcmp("pyramid",py.s[5])==0){ cell.pyr2() }}')
    
                #if 'hippocampus' in s1:
                #    h('cell.pyr2()')
    
    
                if 'interneuron' in s1:
                    ie0[cnti] = i
                    ie1[cnti] = 1
                    cnti += 1
                    h.cell.basket()
                    h('cell.polarity=1')
                if 'pyramid' in s1:
                    ie0[cnti] = i
                    ie1[cnti] = 0
                    cnti += 1
                    h.cell.pyr()
                    h('cell.polarity=0')
                #h.cell.basket()
                #h('cell.polarity=1')
    
    
                h('if(strcmp("interneuron",py.s1[5])==0){ inter_list.append(cell) }')
                h('if(strcmp("interneuron",py.s1[5])==0){ cell.basket()}')
                h('if(strcmp("interneuron",py.s1[5])==0){ print "inter neuron"}')
    
                h('if(strcmp("aspiny",py.s1[5])==0){ aspiny_list.append(cell) }')
    
                h('if(strcmp(cell.population,"hippocampus")==0){ hipp.append(cell)  }')
    
                h('if(strcmp(cell.population,"neocortex")==0){ neoc.append(cell) }')
    
                h('strdef cellposition')
                h('sprint(cellposition,"%s%d%s%d%s%d%s","cell.position(",py.float(py.s1[0]),",",py.float(py.s1[1]),",",py.float(py.s1[2]),")")')
                h('print cellposition')
                h('execute(cellposition)')
                
                #h.pc.set_gid2node(int(i),RANK)
                h('pc.set_gid2node(int(py.i), pc.id)')  # // associate gid i with this host
    
                h('cell.soma[0] nc =  new NetCon(&v(0.5), nil)')
                
                h('pc.cell(int(py.i), nc)')  # //
                #h.pc.cell(int(i),h.nc)
    
                h('cells.append(cell)')
                h('gidvec.append(py.i)')
                #h.gidvec.append(i)
                gidvec.append(i)
                print i
                h('print py.i, " py.i", cnt, " cnt"')
                h('cnt+=pc.nhost')
    
                cnt += 1
                #strangley i can update in this environment but not in wirecells.
        #destroy the allrows list to free up memory, then return the empty list.
        #I can't seem to destroy allrows here without wrecking something in rigp later.
        #allrows=[]
        #pdb.set_trace()
        return (h.cells, allrows2, ie0, ie1)

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
        
        cell.simplify_axon()
        for sec in cell.axonal:
            sec.L = 30
            sec.diam = 1
            sec.nseg = 1 + 2 * int(sec.L / 40)
        cell.axon[0].connect(cell.soma[0], 0.5, 0)
        cell.axon[1].connect(cell.axon[0], 1, 0)
        h.define_shape()
    
    def load_cell_parameters(self, cell, type_index):
        from neuron import h
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
                self.nclist.append(nc)


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
        
