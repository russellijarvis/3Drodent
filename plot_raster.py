
#import matplotlib
#matplotlib.use('Agg')

import pylab as plt

"""
#now the plotting code:
plt.figure(1) # to increase the size of the graph use: mp.figure(1, figsize=(20,20))
plt.scatter(tvec, idvec, c='b', marker='+')

plt.xlabel('spike times')
plt.ylabel('n=neuron')
plt.title('Raster Plot')
#plt.savefig("ff.svg") # works
#plt.draw()
plt.show() #displays the plot immediately after the simulation is over comment if this is not needed
#plt.close()
#times=[len(idvec)]
dt=0.025
print h.dt
print "need some way of scaling by dt"
"""
plt.figure(2)
for i in range(0,int(num_cells)):
  plt.plot(time_courses[i]) # The length of this vector is time / samples (tstop/dt)
  plt.hold(True)

plt.show()
plt.hold(False)



plt.figure(3)
for i in range(0,int(num_cells)):
  plt.plot(trains[i]) # The length of this vector is time / samples (tstop/dt)
  plt.hold(True)

plt.show()
plt.hold(False)


execfile('spike_distance.py')
ti = 0
tf = int(h.tstop)
#multivariate_spike_distance(timespython, 0,1000, 1000)
t, Sb = multivariate_spike_distance(timespython, ti, tf, 10)

plt.figure(4)
plt.title('Spikes')
#plt.subplot(212)
plt.plot(t,Sb,'k')
plt.xlim([0, tf])
plt.ylim([0, 1])
plt.xlabel("Time (ms)")
plt.title("Multivariate SPIKE distance")
plt.show()

#plt.plot(times[:]) 
#im=plt.imshow(times)
#imshow(times)

