from bsmart import timefreq, pwcausalr
from scipy import array, size
import pyhoc  
def get_sgc(indexi,indexj):
    vec1=downsample(time_courses[int(indexi)],oldrate=40000,newrate=200)
    vec2=downsample(time_courses[int(indexj)],oldrate=40000,newrate=200)
        
    """
    order=10
    rate=200
    maxfreq=0
    
    if maxfreq==0: F=timefreq(outj,rate) # Define the frequency points
    else: F=array(range(0,maxfreq+1)) # Or just pick them
    npts=size(F,0)

   
    h("for k=0,1 vb[k] = new Vector() //vb stands for vector binned.")
    h("vb[0].hist(times[py.indexi],0,(py.maxt+py.binsz-1)/py.binsz,py.binsz)")
    h("vb[1].hist(times[py.indexj],0,(py.maxt+py.binsz-1)/py.binsz,py.binsz)")

    h("vo1=normte(vb[0],vb[1],20)")
    #return F,pp[0,:],cohe[0,:],Fx2y[0,:],Fy2x[0,:],Fxy[0,:]
    
    """
    order=15
    npts=len(vec1)
    f_s=200
    n=len(vec1)
    max_freq=(f_s)/2
    p=15
    ntrls=1
                          #(x,ntrls,npts,p,fs,freq)
    x=array([vec1,vec2])
    F,pp,cohe,Fx2y,Fy2x,Fxy=pwcausalr(x,1,npts,order,rate,max_freq)
    print "GSC components", np.mean(Fx2y), np.mean(Fy2x) 
    print "GSC difference", np.mean(Fx2y)-np.mean(Fy2x)
    print "GSC difference", np.mean(Fy2x)-np.mean(Fx2y)
    print indexi
    return np.mean(Fx2y), np.mean(Fy2x), np.mean(Fx2y) -np.mean(Fy2x)

def get_sgc_wa(vec1,vec2):
    #vec1=downsample(indexi),oldrate=40000,newrate=200)
    #vec2=downsample(indexj),oldrate=40000,newrate=200)
        
    """
    order=10
    rate=200
    maxfreq=0
    
    if maxfreq==0: F=timefreq(outj,rate) # Define the frequency points
    else: F=array(range(0,maxfreq+1)) # Or just pick them
    npts=size(F,0)

   
    h("for k=0,1 vb[k] = new Vector() //vb stands for vector binned.")
    h("vb[0].hist(times[py.indexi],0,(py.maxt+py.binsz-1)/py.binsz,py.binsz)")
    h("vb[1].hist(times[py.indexj],0,(py.maxt+py.binsz-1)/py.binsz,py.binsz)")

    h("vo1=normte(vb[0],vb[1],20)")
    #return F,pp[0,:],cohe[0,:],Fx2y[0,:],Fy2x[0,:],Fxy[0,:]
    
    """
    order=15
    npts=len(vec1)
    f_s=100# for time bins of 10ms
    n=len(vec1)
    max_freq=(f_s)/2
    p=15
    ntrls=1
                          #(x,ntrls,npts,p,fs,freq)
    #x=array([vec1,vec2])
    #F,pp,cohe,Fx2y,Fy2x,Fxy=granger(vec1,vec2,ntrls,15)#,p,fs,freq)#,200,100)# nfs2[0], nfs2[1] ,
    F,pp,cohe,Fx2y,Fy2x,Fxy=pwcausalr(x,1,npts,order,rate,max_freq)
    print "GSC components", np.mean(Fx2y), np.mean(Fy2x) 
    print "GSC difference", np.mean(Fx2y)-np.mean(Fy2x)
    print "GSC difference", np.mean(Fy2x)-np.mean(Fx2y)
    print indexi
    return np.mean(Fx2y), np.mean(Fy2x), np.mean(Fx2y) -np.mean(Fy2x)

#a=np.zeros([int(h.ntefa.v[2].size()-1)])
#b=np.zeros([int(h.ntefa.v[2].size()-1)])
"""
a=[];b=[]
for i in xrange(0,int(h.ntefa.v[2].size()-1)):
  if(h.ntefa.v[2].x[i]>0.1):
    print 'entropy', h.ntefa.v[2].x[i]
    #a.append(
    #x,y,z=
    get_sgc(h.ntefa.v[0].x[i],h.ntefa.v[1].x[i])
    #) 

    b.append(h.ntefa.v[2].x[i])    
    
p1.plot(a)
p1.plot(b)    
"""    
