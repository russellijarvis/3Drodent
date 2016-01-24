#Transfer entropy estimation
#Greg Ver Steeg and Aram Galstyan 
#{gregv,galstyan}@isi.edu
#August 22,2012
#This code is released under GPLv3. 
#The associated paper "Inferring predictive links in social media using content transfer"
#is available on arxiv
#The paper uses this code along with the package "gensim" (available free online)
#To estimate "content transfer" among users on Twitter
#For request the data used in the paper, contact Sofus Macskassy, sofmac@isi.edu
#Example usage of the function is shown in an example in the body of the program
#(Also requires scipy/numpy, freely available scientific computing packages for python)

from math import log,pi
from scipy.special import digamma,gamma
import numpy.random as nr
import scipy.spatial as ss

def zip2(*args):
  return [sum(sublist,[]) for sublist in zip(*args)]

inf = float("inf")
def palus(points,dvec):
  #use palus,vejmelka 2008 method 
  N = len(points)
  tree = ss.KDTree(points)
  entropy = 0.
  for i in range(N):
    dist = dvec[i]
    num_points = len(tree.query_ball_point(points[i],dist-1e-15,p=inf)) 
    # -1 for center point +1 for est.
    #-1e-15, strictly less
    entropy += -digamma(num_points)/N
  return entropy

def cmi(x,y,z,k=3):
  #Mutual information of x and y, conditioned on z
  #add low intensity noise
  intens = 1e-10 #1./len(x)
  x = [list(p + intens*nr.rand(len(x[0]))) for p in x]
  y = [list(p + intens*nr.rand(len(y[0]))) for p in y]
  z = [list(p + intens*nr.rand(len(z[0]))) for p in z]
  points = zip2(x,y,z)
  #Find nearest neighbors in joint space
  tree = ss.KDTree(points)
  #indice and distance to each k-nn
  lvec = [tree.query(point,k+1,p=inf)[1][k] for point in points]
  dvec = [ss.distance.chebyshev(points[i],points[lvec[i]]) for i in range(len(points))]
  a,b,c,d = palus(zip2(x,z),dvec), palus(zip2(y,z),dvec), palus(z,dvec), -digamma(k)
  return a+b-c-d


if __name__ == "__main__":
  from numpy.linalg import det
  import numpy as np
  import pylab
  import random

  #CMI testing suite
  #generate samples with known CMI
  #and uses Palus estimator
  
  knn = "3"
  nsamples = 400 #Number of estimates for C.I.
  samplo = int(0.025*nsamples) #95% confidence intervals
  samphi = int(0.975*nsamples)
  Ntry = [10,25,50,100,200,400,800] #,1000,2000] #Number of samples to use in estimate

  #Generate points and true entropies for a multivariate Gaussian defined below
  d1 = [1,1,0]
  d2 = [1,0,1]
  d3 = [0,1,1]
  mat = [d1,d2,d3]
  tmat = np.transpose(mat)
  diag = [[3,0,0],[0,1,0],[0,0,1]]
  mean = np.array([0,0,0])
  cov = np.dot(tmat,np.dot(diag,mat)) 
  trueent = -0.5*(3+log(8.*pi*pi*pi*det(cov))) 
  trueent += -0.5*(1+log(2.*pi*cov[2][2])) #z sub
  trueent += 0.5*(2+log(4.*pi*pi*det([[cov[0][0],cov[0][2]],[cov[2][0],cov[2][2]]] ))) #xz sub
  trueent += 0.5*(2+log(4.*pi*pi*det([[cov[1][1],cov[1][2]],[cov[2][1],cov[2][2]]] ))) #yz sub
  print 'true CMI', trueent
  trueent = 0.5*(1+log(2.*pi*cov[0][0])) #x sub
  trueent += 0.5*(1+log(2.*pi*cov[1][1])) #y sub
  trueent += -0.5*(2+log(4.*pi*pi*det([[cov[0][0],cov[0][1]],[cov[1][0],cov[1][1]]] ))) #xz sub
  print 'true MI', trueent
  allpoints = nr.multivariate_normal(mean,cov,5*max(Ntry))

  #ESTIMATION should produce points in Fig 2(a) of paper
  for NN in Ntry:
    tempent = []
    for j in range(nsamples):
      points = random.sample(allpoints,NN)
      #We have N samples of a d-dimensional variable.
      #So x is a list of N elements. Each element is a list of d coordinates
      x = [point[:1] for point in points] 
      y = [point[1:2] for point in points] 
      z = [point[2:] for point in points] 
      #THIS IS THE CALL TO ESTIMATE CMI
      ###############
      tempent.append(cmi(x,y,z,k=int(knn)))
      ###############
    tempent.sort()
    tempmean = np.mean(tempent)
    print "For "+str(NN)+" samples, avg CMI:", tempmean
    print "CI:",tempent[samplo],tempent[samphi]
