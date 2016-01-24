from IPython.parallel import Client
c = Client(profile='mpi')
view = c[:]
view.activate() # enable magics
view.run('initpar.py')
