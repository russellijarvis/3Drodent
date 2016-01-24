from mpi4py import MPI
import numpy as np

def psum(a):
    locsum = np.sum(a)
    rcvBuf = np.array(0.0,'d')
    MPI.COMM_WORLD.Allreduce([locsum, MPI.DOUBLE],
        [rcvBuf, MPI.DOUBLE],
        op=MPI.SUM)
    return rcvBuf


from IPython.parallel import Client

c = Client(profile='mpi')

view = c[:]

view.activate() # enable magics

# run the contents of the file on each engine:
view.run('psum.py')

view.scatter('a',np.arange(16,dtype='float'))

view['a']

#%px totalsum = psum(a)

view['totalsum']
