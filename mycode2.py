
#h.xopen("/home/zaza3/Downloads/trunk/examples/expericomp20140421/mycode3.hoc")

from mpi4py import MPI
from neuron import h
import numpy as np
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
h('objref py')
h('py = new PythonObject()')
data1 = np.arange(4, dtype='f')
h('objref coords')
h('coords = new Vector(3)')
h('coords.x[0]=1')
h('coords.x[1]=2')
h('coords.x[2]=3')
coords=np.array(h.coords.to_python(),dtype=np.float32)
print coords
#if rank == 0:
#    data = {'a': 7, 'b': 3.14}
#    comm.send(data, dest=1, tag=11)
#    print "Message sent, data is: ", data
#elif RANK == 1:
#    data = comm.recv(source=0, tag=11)
#    print "Message Received, data is: ", data
if RANK == 0:
    #data = np.arange(1000, dtype='f')
    comm.Send([coords, MPI.FLOAT], dest=1, tag=77)
elif RANK == 1:
    #data = np.empty(1000, dtype='f')
    comm.Recv([coords, MPI.FLOAT], source=0, tag=77)
    print "Message Received, data is: ", coords

# automatic MPI datatype discovery
if RANK == 0:
    data = np.arange(100, dtype=np.float32)
    comm.Send(data, dest=1, tag=13)
elif RANK == 1:
    data = np.empty(100, dtype=np.float32)
    comm.Recv(data, source=0, tag=13)
    print "Message Received, data is: ", data


if RANK == 0:
    data = {'key1' : [7, 2.72, 2+3j],
            'key2' : ( 'abc', 'xyz')}
else:
    data = None
    data = COMM.bcast(data, root=0)
    print "bcast finished and data \
    on RANK %d is: "%COMM.RANK, data

if RANK == 0:
    data = {'key1' : [7, 2.72, 2+3j],
            'key2' : ( 'abc', 'xyz')}
else:
    data = None
    data = COMM.bcast(data, root=0)
    print "bcast finished and data \
    on RANK %d is: "%COMM.RANK, data


quit()


for i in range(RANK, int(len(allrows)), SIZE):

    for y in xrange(0,int(ncell)):
        h('''forsec cells.o(py.y).spk_trig_ls{ for(x){ if((x>0)&&(x<1)){ py.data1[0]=x_xtra   py.data1[1]=y_xtra   py.data1[2]=z_xtra   py.data1[3]=x''')
		   

        if RANK == 0:
            data = data1
        else:
            data = None
            data = COMM.scatter(data1, root=0)
        #assert data == (RANK+1)**2
            data = COMM.gather(data1, root=0)
            if RANK == 0:
                data = data1#np.arange(1000, dtype='i')
        #COMM.Send([data, MPI.INT], dest=1, tag=77)
                COMM.Send(data, dest=1, tag=13)

            elif RANK == 1:
                data = data1
                COMM.Recv(data, source=0, tag=13)

        #COMM.Recv([data, MPI.INT], source=0, tag=77)
"""
	    for z=0,cells.count-1{
		if(y!=z){
		    ////print y, z, "is this even executing?"
		    wire_distances(cells.o(z),cells.o(y))
		}
		
		
	    }//OLD y spk_rx_ls is done with and will never be iterated.
	    //cells.o(y).spk_rx_ls=new SectionList()
	    
	    
	}

proc find_distances(){local sx1, sy1, sz1, sx0, sy0, sz0, polarity, tri, dist 
    //Now objects such as conn_matrix object is only available in this method because procedures protect them.
    
    
    
    nrnpython("import csv")
    chdir(workingdir)
    syn_//, gabablist, ampalist//, prng//, cnt_r  // now expects xyc coords as arguments
    sx1=sy1=sz1=sx0=sy0=sz0=0
    
    strdef sec_string1 //needs more meaningful names. Secstring1, should be  precell
    strdef sec_string2 //secstring2 should be postcell
    strdef cell_here
    strdef executable //I know bad names for future debugging
    strdef executable2
    strdef executable3
    strdef executable4
    strdef tail_src
    strdef gabalists, ampalists
    strdef tail_target
    strdef filename
    
    sprint(filename,"%s%d","synapse_file.hoc",large_scale)
    strdef pre_src
    
    
    forsec $o7.spk_rx_ls{
	
	////print secname() 
	//NEW DUPLICATE CODE PLUS CONDITION
	sec_string1=$s4 
	num_seg=$5
	sec_string2=secname() 
	strobj.head(sec_string2, "].", cell_here) 
	strobj.head(sec_string1, "].", pre_src) 
	strobj.tail(pre_src,"Cell",tail_src)
	strobj.tail(cell_here,"Cell",tail_target) //tgt	
	strobj.right(tail_target, 1)
	strobj.right(tail_src, 1)
	if(strcmp(tail_src,tail_target)==0){ break }  //Ie if the source cell and target are
	if(Cell[py.int(tail_src)].div.x[py.int(tail_target)]>3){ break }
	
	
	if (ismembrane("xtra")) { //its something about xtra that only has segment indexs at 0.5
	    //for (x,0) {
	    for(x){
		if((x>0)&&(x<1)){
		    
		    
		    
		    ////
		    // Cleft Distance. Ie is r sufficiently small such that the presence of vesicles is warranted.
		    //// wedssaq 
		    
		    r = (sqrt((x_xtra(x) - $1)^2 + (y_xtra(x) - $2)^2 + (z_xtra(x) - $3)^2))
		    py.cleft_hist.append(r) 
		    sum_r=sum_r + r  
		    //checkif=(r<minr) //check distances
		    ////print checkif
		    ////print "is this even executing?"
		    if(large==1){
			checkif=(r<minr) //check distances
			
		    }else{
			checkif=(r<minrs)
		    }
		    
		    if(checkif){  
			
			
			// Since the threshold for a synaptic cleft has been satisfied.
			// Code below calculates the Somatic Distance
			//// 
			
			/* Remove duplicate  
			sec_string1=$s4 
			num_seg=$5
			sec_string2=secname() 
			
			strobj.head(sec_string2, "].", cell_here) 
			strobj.head(sec_string1, "].", pre_src) 
			strobj.tail(pre_src,"Cell",tail_src)
			strobj.tail(cell_here,"Cell",tail_target) //tgt	
			strobj.right(tail_target, 1)
			strobj.right(tail_src, 1)
			*/
			
			sx0= Cell[py.int(tail_src)].soma[0].x_xtra()
			sy0= Cell[py.int(tail_src)].soma[0].y_xtra()
			sz0= Cell[py.int(tail_src)].soma[0].z_xtra()
			sx1= Cell[py.int(tail_target)].soma[0].x_xtra()
			sy1= Cell[py.int(tail_target)].soma[0].y_xtra()
			sz1= Cell[py.int(tail_target)].soma[0].z_xtra() 
			rs = (sqrt(( sx0 - sx1)^2 + ( sy0 - sy1)^2 + ( sz0 - sz1)^2))
			////
			// 2/3rds of values are above 6.0272221
			//(sum_rs/py.len(py.dist_hist))-(((sum_rs/py.len(py.dist_hist))*2)/6)
			//	6.0272221 
			//
			//9.0408332
			//  The pre synaptic transmitter polarity is determined by rs.
			////
			//Only cells who have clefts close to each other shall obtain a type. If cells never have clefts that surpass minimum threshold 
			//The cells will never even obtain polarity.
			if(Cell[py.int(tail_src)].type==-1){ //If cell type not already classified.
			    if(rs>=5){//if(rs>=0.10408076){
				Cell[py.int(tail_src)].type=1 //cell type is assigned excitatory for the furtherest 2/3 of somatic distances.
			    }else if(rs<5){//       if(rs<0.10408076){	
				Cell[py.int(tail_src)].type=0 //cell type is assigned inhibitory for the closest 1/3 of somatic distances.
			    }
			    
			}
			// //print rs,"This is the somatic distance"
			py.rs=rs
			nrnpython("dist_hist.append(rs)")
			sum_rs=sum_rs+rs
			
			//OLD may wish to restore.
			//if(strcmp(tail_src,tail_target)!=0){  //Ie if the source cell and target are different. No neural short circuits.
			
			
			storex=x
			//print num_seg, x
			
			
			
			// At this point cells which are neither type 3, 4 or 5 may pass.
			
			//print Cell[py.int(tail_src)].div.x[py.int(tail_target)]
			//write to sparse Matrix's.
			
			//if(Cell[py.int(tail_src)].type==1){ //If excitatory cell type
			if(Cell[py.int(tail_src)].num_type==5||(Cell[py.int(tail_src)].num_type==4)){//Here I am bypassing if excitatory cell type was 
			    //deterined by somatic distance. Be aware that this method of assigning synapse polarity leaves many neurons unaccounted for.
			    
			    
			    
			    //If excitatory cell type
			    //For varying STDP synapses use the below. which combines AMPA and NMDA
			    //sprint(synapse_post, "%s%s%g%s", sec_string2, " syn_ =  new AMPA(",storex,")")
			    // For fixed gain synapses use the below.
			    if(plastic==1){
				
				if(cnt_f%2==0){
				    
				    //sprint(synapse_post, "%s%s%g%s", sec_string2, " syn_ =  new STDPE2(",storex,")")
				    sprint(synapse_post, "%s%s%g%s", sec_string2, " syn_ =  new NMDA(",storex,")")
				    execute(synapse_post) //put a post synapse at the target 
				    sprint(ampalists,"%s%s", cell_here,"].nmdalist.append(syn_)")
				    execute(ampalists) //Now the synapse belongs to two lists, harder to destroy but thats okay.
				    
				}else{
				    sprint(synapse_post, "%s%s%g%s", sec_string2, " syn_ =  new AMPA(",storex,")")
				    ////print "fid executes okay!"
				    
				}
				
			    }else{
				sprint(synapse_post, "%s%s%g%s", sec_string2, " syn_ =  new AMPA(",storex,")")
				////print fid(Cell[py.int(tail_target)],sec_string2), "fid executes okay!"
				
				
			    }
			    
			    execute(synapse_post) //put a post synapse at the target 
			    
			    sprint(ampalists,"%s%s", cell_here,"].ampalist.append(syn_)")
			    execute(ampalists)
			    
			    //plastic synapses are added to two lists so they are harder to destroy.
			    
			    
			    sprint(synapse_pre, "%s%s%g%s", sec_string1," nc = new NetCon (&v(",num_seg,"),syn_)")
			    // sprint(synapse_pre, "%s%s", sec_string1,"     nc = new NetCon (&v (0.5), syn_)")
			    
			    
			    
			    execute(synapse_pre)
			    ////print synapse_pre
			    nc.threshold = 30
			    if(plastic==1){  
				if(cnt_f%5==0){
				    nc.weight = 0.0035
				    ptb.append(nc.weight,cnt_f)
				}else{
				    nc.weight = (r+rs)/iw// //was only divisor 100
				}
				
			    }else{   
				nc.weight = (r+rs)/iw// //was only divisor 100
			    }
			    nc.delay = (r + rs)/delay// Ie the delay is proportional to the cleft + soma distance 
			    nclist.append(nc)
			    //print py.int(tail_target), "tail_target"
			    //print py.int(tail_src), "tail_src"
			    
			    // This code has to be inside these two if conditions because If the cell is of type -1 
			    // tail_target, and tail src strings won't obtain a real value, and will crash code.
			    py.targets=(py.int(tail_target)) // the post synapse.
			    py.sources=(py.int(tail_src)) //the pre synapse.
			    //Why am I am making a list of these values.? I should be just using the values.
			    
			    //print py.int(tail_target)
			    //print py.int(tail_src)
			    
			    py.weights=py.int(Cell[py.int(tail_src)].type)
			    
			    //py.pyf.write(py.int(tail_target),'\t',py.int(tail_src),'\t',nc.weight)
			    
			    
			    //for key, value in orfdict.iteritems():
			    chdir(workingdir)
			    //f//print("%s%s%s%s%s%g%s%s%g%s",synapse_post,"\n",synapse_pre,"\n","nc.weight=",nc.weight,"\n","nc.delay=",nc.delay,"\n") //dont bother with    //print py.int(tail_target)
			    //print py.int(tail_src)
			    //f//print("%s%s%s%s%g%s%g%s",synapse_post,"\n",synapse_pre,"\n",nc.weight,"\n",nc.delay,"\n") //dont bother with new lines. The string start and stop is encoded.
			    
			    //adjnqs.append(synapse_post,synapse_pre,nc.weight)//fill spread sheet row by row
			    //NQS won't accept string arguments unless you declare them in advance.
			    polarity=1
			    
			    if(get_dist!=1){
				adjnqs.append(py.int(tail_target),py.int(tail_src),nc.weight[4],nc.delay,polarity,sec_string1,secname())
			    }
			    if(get_dist==1){
				
				tri=transfer1(sec_string2,Cell[py.int(tail_target)],20)
				//print tri, "tri"
				dist=fid(Cell[py.int(tail_target)],sec_string2)
				//print dist, "dist" 
				
				adjnqs.append(py.int(tail_target),py.int(tail_src),nc.weight[4],nc.delay,polarity,sec_string1,secname(),dist,tri)
				
				//  
				
				//
				
				//("presynapse_cell","postsynapse_cell", "weight","delay","polarity","presection","postsection","distance","input_impedence") //The titles 
				py.presec=py.str(sec_string2)
				py.postsec=py.str(sec_string1)
				
				nrnpython("sparsem[0].append(presec)")
				nrnpython("sparsem[1].append(postsec)") 
				
			    }
			    
			    //nrnpython("MatrixG[int(targets[len(targets)-1])][int(sources[len(sources)-1])]=Matrix[int(targets[len(targets)-1])][int(sources[len(sources)-1])]+1")
			    
			    
			    
			    
			    
			    
			    
			    nrnpython("mbin[int(sources)][int(targets)] = 1")
			    
			    nrnpython("Matrix[int(sources)][int(targets)] =Matrix[int(sources)][int(targets)] + 1")// #The degree of connections will be graphed as 
			    nrnpython("MatrixA[int(sources)][int(targets)] =MatrixA[int(sources)][int(targets)] + 1")// #The degree of connections will be graphed as 
			    // nrnpython("WeightsA[int(sources)][int(targets)] =WeightsA[int(sources)][int(targets)] + 1")// #The degree of connections will be graphed as 
			}
			
			//        nrnpython("Weights[int(targets[len(targets)-1])][int(sources[len(sources)-1])]=Weights[int(targets[len(targets)-1])][int(sources[len(sources)-1])]+weights[len(weights)-1]")
			
			conn_matrix.x[py.int(tail_target)][py.int(tail_src)] = conn_matrix.x[py.int(tail_target)][py.int(tail_src)]+1
			conn_matrixA.x[py.int(tail_target)][py.int(tail_src)] = conn_matrixA.x[py.int(tail_target)][py.int(tail_src)]+1
			///
			Cell[py.int(tail_src)].div.x[py.int(tail_target)]=Cell[py.int(tail_src)].div.x[py.int(tail_target)]+1
			Cell[py.int(tail_target)].conv.x[py.int(tail_src)]=Cell[py.int(tail_target)].conv.x[py.int(tail_src)]+1
			
			
			
			
			
			//       if(Cell[py.int(tail_src)].type==0){ //If excitatory cell type
			if((Cell[py.int(tail_src)].num_type==3)){ //If the neuron type is interneuron.
			    //I should be aware that this method of assigning synapse polarity leaves many neurons unaccounted for.
			    
			    sprint(synapse_post, "%s%s%g%s", sec_string2, " syn_ = new  GABAa(",storex,")") //MyExp2Syn(",storex,")")
			    
			    
			    
			    
			    execute(synapse_post) //put a post synapse at the target
			    sprint(gabalists,"%s%s", cell_here,"].gabalist.append(syn_)")
			    execute(gabalists)
			    
			    
			    sprint(synapse_pre, "%s%s%g%s", sec_string1," nc = new NetCon (&v(",num_seg,"),syn_)")
			    execute(synapse_pre)
			    nc.threshold = 5 //To much for GABAergic synapse.
			    nc.weight = (r+rs)/iw// kills network. with such low value caused by divisor: 4000 //was (r)/100 when everything worked well.
			    
			    //weights should be loaded from the final state of the previous network.     
			    nc.delay = (r+rs)/delay// Ie the delay is proportional to the cleft + soma distance 
			    //(note not arc length so only approximation)
			    
			    nclist.append(nc)
			    
			    //print py.int(tail_target), "tail_target"
			    //print py.int(tail_src), "tail_src"
			    //print py.int(tail_target)
			    //print py.int(tail_src)
			    // This code has to be inside these two if conditions because If the cell is of type -1 
			    // tail_target, and tail src strings won't obtain a real value, and will crash code.
			    py.targets=py.int(tail_target) // the post synapse.
			    
			    py.sources=(py.int(tail_src)) //the pre synapse.
			    py.weights=(py.int(Cell[py.int(tail_src)].type))
			    chdir(workingdir)
			    
			    //x_xtra pertains to post.
			    //The arguments pertain to pre.
			    //f//print("%s%s%s%s%s%g%s%s%g%s",synapse_post,"\n",synapse_pre,"\n","nc.weight=",nc.weight,"\n","nc.delay=",nc.delay,"\n") //dont bother with new lines. The string
			    //     f//print("%s%s%s%s",synapse_post,"\n",synapse_pre,"\n") //dont bother with new lines. The string start and stop is encoded.
			    //	           f//print("%s%s%s%s%d%s",synapse_post,"\n",synapse_pre,"\n",nc.weight,"\n") //dont bother with new lines. The string start and stop is
			    //  py.weights.append(nc.weight*-1)
			    //  nrnpython("ldict = {'postsecname': h('sec_string2'), 'xpst': h('x_xtra(x)'), 'ypst': h('y_xtra(x)'), 'zpst': h('z_xtra(x)'),'pressecname': h('sec_string1'), 'xpre': h('$1'), 'ypre': h('$2'), 'zpre': h('$3') }")
			    
			    
			    
			    
			    polarity=-1
			    
			    if(get_dist!=1){
				adjnqs.append(py.int(tail_target),py.int(tail_src),nc.weight[4],nc.delay,polarity,sec_string1,secname())
			    }
			    if(get_dist==1){
				////print "got here", fid(Cell[py.int(tail_target)],sec_string2)
				
				//
				
				tri=transfer1(sec_string2,Cell[py.int(tail_target)],20)
				//print tri, "tri"
				dist=fid(Cell[py.int(tail_target)],sec_string2)
				//print dist, "dist" 
				
				adjnqs.append(py.int(tail_target),py.int(tail_src),nc.weight[4],nc.delay,polarity,sec_string1,secname(),dist,tri)
				
				
				
				//  
				//transfer1(secname(),cells.o(py.int(tail_target)),100)  
				//
				py.presec=py.str(sec_string2)
				py.postsec=py.str(sec_string1)
				py.cnt_f=py.int(cnt_f)
				nrnpython("sparsem[0].append(presec)")
				nrnpython("sparsem[1].append(postsec)")
			    }
 			    //   nrnpython("import neuron as h")
			    // nrnpython("Matrix[int(targets[len(targets)-1])][int(sources[len(sources)-1])]=Matrix[int(targets[len(targets)-1])][int(sources[len(sources)-1])]+1")
			    
			    nrnpython("mbin[int(sources)][int(targets)] = 1")
			    nrnpython("Matrix[int(sources)][int(targets)] =Matrix[int(sources)][int(targets)] + 1")
			    
			    nrnpython("MatrixG[int(sources)][int(targets)] =MatrixG[int(sources)][int(targets)] + 1")
			    // It is not actually these matrices that are causing the memory errors. A massive amount of NEURON is written in Python. Error messages of Python running out of memory, are the internal NEURON functions/methods running out of memory 
			    //Its probable that no Gabargic connections are being made because of the sparesness of neurons used. Distance between neurons too large.
			    
			    
			    //Makes sense to store weights and delays as sparse matrix, but values get lost in a conventional matrix.
			    
			    //    nrnpython("WeightsG[int(sources)][int(targets)] =WeightsG[int(sources)][int(targets)] + 1")
			    
			    
			    //These matrices are the wrong way around compared to the pythonic matrices.
			    
			    
			    conn_matrix.x[py.int(tail_target)][py.int(tail_src)] = conn_matrix.x[py.int(tail_target)][py.int(tail_src)]+1
			    
			    conn_matrixG.x[py.int(tail_target)][py.int(tail_src)] = conn_matrixG.x[py.int(tail_target)][py.int(tail_src)]+1
			    
			    ///
    			    Cell[py.int(tail_src)].div.x[py.int(tail_target)]=Cell[py.int(tail_src)].div.x[py.int(tail_target)]+1
			    Cell[py.int(tail_target)].conv.x[py.int(tail_src)]=Cell[py.int(tail_target)].conv.x[py.int(tail_src)]+1
			    
      			}
			strdef execution
			strdef execution1
			strdef prelists
			strdef postlists
			
 			
			sprint(executable,"%s%s", cell_here,"].synlist.append(syn_)")
			//sprint(prelists,"%s%s%s%s", cell_here,"].prelist.append(",sec_string1,")")
			//sprint(postlists,"%s%s%s%s", cell_here,"].postlist.append(",secname(),")")
			
			//sprint(postlists,"%s%s", secname()," a = new SectionRef()")
			execute(executable)
			//execute(prelists)
			//execute(postlists)
			// The following is performed in run_stats now
			
			//The reason that previously the sum of all (valid) network connections has been less than the sum	of pre and post synapse pairs. Is the same reason why post and presynaptic synapse numbers have not always been equal, and that is because synapses where forcibly put in the centre of each segment and that is not always a reasonable place to put a synapse. It is also not faithful to the wiring distance algorithm which more often allocates synapses to segments of sections which are not in the centre of the section.
			//Later I can check for consistancy between this way of finding sources and targets, and NEURON built in methods that I systematically exploit in param.hoc	
			//Faster to construct this matrix here all in one pass.
			//Get the coordinates of the spike initation zone of the soma, or the nearest dendrite.
    			strdef check_soma0, check_soma1 
			////
			//py.somatic_histogram_py.append(rs)
			//the somatic_histogram_py has become null for reasons its hard to explain
			//py.len(py.somatic_histogram_py)
			////
 			cnt_conn+=1
			cnt_f+=1
			py.cnt_f=py.cnt_f+1
			
			
		    }
		}
	    }
	}
    }
    //If this fails it means there were no connections
    //Its better that this fails, than runs.
    
    // return conn_matrix
} 

//quit()

proc wire_distances(){
    cnt_segments=0
    
    forsec $o1.spk_trig_ls{ //localobj sec_string
	
	////print secname() 
	strdef sec_string1
	
	for (x){
	    if((x>0)&&(x<1)){
		if(ismembrane("xtra")){
		    cnt_segments+=1 
		    sec_string1=secname()
		    seg_num=x
		    
		    
		    
		    
		    
		    //This is hugely erroneous. There needs to be only one external conn_matrix not one that is remade everytime this procedure is enterered;
		    find_distances(x_xtra(x), y_xtra(x), z_xtra(x), sec_string1, seg_num, average_distance,$o2)//
		}
	    } 
	} 
    }
    
    //nrnpython("file=open('connection_matrix.txt',w)")
    //nrnpython("np.save('connection_matrix.txt',Matrix)")
}

proc reapply(){
    for indexs=0,cells.count-1{
	cells.o(indexs).div=conn_matrix.getrow(indexs) //rows are sources
	cells.o(indexs).conv=conn_matrix.getcol(indexs) //columns are targets
    }
}
"""
