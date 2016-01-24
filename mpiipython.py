

from IPython.parallel import Client
c = Client()
c.ids
#c[:].apply_sync(lambda: "Hello, world!")
view = c[:]
view.activate()
#view.run('memb_pot.py')
