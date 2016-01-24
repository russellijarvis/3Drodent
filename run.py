#%run run.py
#pdb.run('execfile("run.py")')

from neuron import h
#feed in parameters to the program somehow. Especially to the param file.
h.xopen("init.hoc")

#No this style of program needs to happen in BASH instead.

del h
from neuron import h
h.xopen("init.hoc")
# #Interpretor stay open for debugging.

#change input parameters to init.hoc here and open it again. 

# run a different simulation.

#keep opening and closing it with different parameters to emulate a batch effect.
