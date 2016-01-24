objref ghist
ghist = new Graph()
objref vo1, vo2, vb2[2], vb[2]
nrnpython("from bsmart import granger # Load the Granger calculation tool")
//proc ent_table_sqsbf(){local cnt_size, i, j, k, largest_size, cnt localobj vo1, vo2, vb2[2], vb[2]
nrnpython("execfile('pyhoc.py')") 

  for i=0,cells.count-1{//count up
  py.ii=py.int(i)
    j=0
    //Lets make just one column 0,60 versus 0.
 
       h("maxt=spk[j].max")
        h("printf("maxt=%g\n",maxt)")
        
        if(spk[i].max>maxt)maxt=spk[i].max
        h("printf("maxt=%g\n",maxt)")
        h("binsz=10")
        if(spk[i].size>0){
           if(maxt>0){
          
           h("for k=0,1 vb[k] = new Vector() //vb stands for vector binned.")
           h("vb[0].hist(spk[i],0,(maxt+binsz-1)/binsz,binsz)")
           h("vb[1].hist(spk[j],0,(maxt+binsz-1)/binsz,binsz)")
           
           h("for k=0,1 vb2[k] = new Vector()") //vb stands for vector binned.
           h("vb2[0].hist(times[i],0,(maxt+binsz-1)/binsz,binsz)")
           h("vb2[1].hist(times[j],0,(maxt+binsz-1)/binsz,binsz)")
                
           oi=downsample(time_courses[int(ii)],oldrate=40000,newrate=200)
           oj=downsample(time_courses[0],oldrate=40000,newrate=200)
           F,pp,cohe,Fx2y,Fy2x,Fxy=granger(np.array(oi.to_python()),np.array(oj.to_python()),20)
           
           oid=downsample(h.vb[0].to_python(),oldrate=40000,newrate=200)
           ojd=downsample(h.vb[1].to_python(),oldrate=40000,newrate=200)
           F,pp,cohe,Fx2yd,Fy2xd,Fxy=granger(np.array(oid),np.array(ojd),20)
           
           //vb[0].graph(ghist,0.1)
           //vb[1].graph(ghist,0.1)

          vo1=normte(vb[0],vb[1],20)
          vo2=normte(vb[1],vb[0],20)
          te01=vo1.x(2)
          te10=vo2.x(2)
          if(vo1.x(4)<=0)vo1.x(4)=1
          if(vo2.x(4)<=0)vo2.x(4)=1
          if(te01>0 || te10>0) {
          pf01=(te01-te10)/(te01+te10) //difference. normalised
          pf10=(te10-te01)/(te01+te10) //difference. normalised 
          } else {
          pf01=pf10=0
          }


          print "pref dir nent", pf01
          nrnpython("print 'direction SGC', np.mean(Fx2y)-np.mean(Fy2x)")
          
          print "transfered amount", vo1.x(2)
          nrnpython("print 'direction SGC', np.mean(Fx2y)")
          //nrnpython("print 'direction SGC', Fy2x")
          //vo1=normte(vb[0],vb[1],20)
          //d=vo1.x(2)
          
          
          
          //storval=GetTENQ(vb[0],vb[1],20,i,j)

          py.dirin=py.float(storval.v[5].x[0])
          nrnpython("msgcv[ii][0]=np.mean(Fx2y)")
          nrnpython("dir[ii][0] = float(dirin)")


     
          print d, "double"
          py.i=py.int(i)
     
          py.j=py.int(j)
          if(storval.v[8].x[0]>0.5){ 

            ntefa.append(i,j,vo1.x(2)) 
          }

        

          }
       }
    }
   // }
  }
  //nrnpython("plt.show()")
  chdir(results_dir)
  nqte.sv(svnqte)
  chdir(workingdir)
  //nqte.sv(svnqte)
}
//ent_table_sqsbf()
