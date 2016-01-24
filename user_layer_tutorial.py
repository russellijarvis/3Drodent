import nineml.user_layer as nineml
    
exc_celltype = nineml.SpikingNodeType(
                     name="Excitatory neuron type",
                     definition="IaF_tau.xml",
                     )
    

parameters = {"membraneCapacitance": (1.0, "nF"),
               "membraneTimeConstant": (20.0, "ms"),
               "refractoryTime": (5.0, "ms"),
               "threshold": (-50.0, "mV"),
               "restingPotential": (-65.0, "mV"),
               "resetPotential": (-70.0, "mV")}

exc_celltype = nineml.SpikingNodeType(
                     name="Excitatory neuron type",
                     definition="http://svn.incf.org/svn/nineml/catalog/neurons/IaF_tau.xml",
                     parameters=parameters
                     )
    

exc_celltype
exc_celltype.name
exc_celltype.definition.url
exc_celltype.parameters    
ParameterSet({'membraneCapacitance': Parameter(name=membraneCapacitance, value=1.0, unit=nF),
    'refractoryTime': Parameter(name=refractoryTime, value=5.0, unit=ms), 'threshold':
    Parameter(name=threshold, value=-50.0, unit=mV), 'restingPotential':
    Parameter(name=restingPotential, value=-65.0, unit=mV), 'membraneTimeConstant':
    Parameter(name=membraneTimeConstant, value=20.0, unit=ms), 'resetPotential':
    Parameter(name=resetPotential, value=-70.0, unit=mV)})


from lxml import etree
print etree.tostring(exc_celltype.to_xml(), pretty_print=True)
    
#Suppose that the inhibitory neurons will use the same model, but with some
#randomness in their membrane time constants::
tau_distr = nineml.RandomDistribution(
                 name="normal(20.0,3.0)",                
                 definition="http://svn.incf.org/svn/nineml/catalog/randomdistributions/normal_distribution.xml",
                 parameters={'standardDeviation': (3.0, "dimensionless"),
                             'mean': (20.0, "dimensionless")})
'''
inh_celltype = nineml.SpikingNodeType(
    ...                 name="Inhibitory neuron type",
    ...                 reference=exc_celltype.name,
    ...                 parameters={'membraneTimeConstant': (tau_distr, "ms")})
    
Now we have specified that `inh_celltype` should use the same definition and
parameters as `exc_celltype`, except for the membrane time constant. For now,
`inh_celltype` is unresolved: we haven't matched the reference to the actual object::

    >>> inh_celltype.reference
    >>> inh_celltype
    SpikingNodeType(name="Inhibitory neuron type", UNRESOLVED)
    
Normally, you don't need to worry about resolving components -- it is done
automatically when you build a full model -- but we will do it manually for
illustration purposes::

    >>> inh_celltype.resolve(exc_celltype)
    >>> inh_celltype          
    SpikingNodeType(name="Inhibitory neuron type", definition="http://svn.incf.org/svn/nineml/catalog/neurons/IaF_tau.xml")
    >>> inh_celltype.parameters
    ParameterSet({'membraneCapacitance': Parameter(name=membraneCapacitance, value=1.0, unit=nF),
    'refractoryTime': Parameter(name=refractoryTime, value=5.0, unit=ms),
    'threshold': Parameter(name=threshold, value=-50.0, unit=mV),
    'restingPotential': Parameter(name=restingPotential, value=-65.0, unit=mV),
    'membraneTimeConstant': Parameter(name=membraneTimeConstant,
    value=RandomDistribution(name="normal(20.0,3.0)",
    definition="http://svn.incf.org/svn/nineml/catalog/randomdistributions/normal_distribution.xml"), unit=ms), 'resetPotential': Parameter(name=resetPotential, value=-70.0, unit=mV)})
    
OK, we've got some prototypes for the neuron models, now let's create some populations of neurons::

    >>> grid2D = nineml.Structure(
    ...                 name="2D grid",
    ...                 definition="http://svn.incf.org/svn/nineml/catalog/networkstructures/2Dgrid.xml",
    ...                 parameters={'fillOrder': ("sequential", None),
    ...                             'aspectRatioXY': (1.0, "dimensionless"),
    ...                             'dx': (10.0, u"um"), 'dy': (10.0, u"um"),
    ...                             'x0': (0.0, u"um"), 'y0': (0.0, u"um")})
    >>> exc_cells = nineml.Population(
    ...                 name="Excitatory cells",
    ...                 number=100,
    ...                 prototype=exc_celltype,
    ...                 positions=nineml.PositionList(structure=grid2D))

(If your terminal won't let you type "um" (Alt-m), "um" works just as well). Here
we have specified that the excitatory cells should be laid out on a square grid.
The actual positions won't be calculated until the model is simulated. We can
also specify the cell positions explicitly as a list of (x,y,z) coordinates::

    >>> positions = [(10.0*x,0.0,0.0) for x in range(25)]
    >>> inh_cells = nineml.Population(
    ...                 name="Inhibitory cells",
    ...                 number=25,
    ...                 prototype=inh_celltype,
    ...                 positions=nineml.PositionList(positions=positions))
    
A population cannot exist by itself: we have to add it to a `Group`::

    >>> network = nineml.Group("Network")
    >>> network.add(exc_cells, inh_cells)
    
To create a network, we need to create some more components: for the post-synaptic
response, for the connection type (management of weights, delays, etc.) and the
connectivity rule::

    >>> connection_rule = nineml.ConnectionRule(
    ...                 name="random connections",
    ...                 definition="http://svn.incf.org/svn/nineml/catalog/connectionrules/fixed_probability.xml",
    ...                 parameters={'p_connect': (0.1, "dimensionless")})
    >>> exc_psr = nineml.SynapseType(
    ...                 name="Excitatory post-synaptic response",
    ...                 definition="http://svn.incf.org/svn/nineml/catalog/postsynapticresponses/exp_g.xml",
    ...                 parameters={'decayTimeConstant': (5.0, "ms"), 'reversalPotential': (0.0, "mV"), 'unitAmplitude': (0.1, "nS")})
    >>> inh_psr = nineml.SynapseType(
    ...                 name="Inhibitory post-synaptic response",
    ...                 reference=exc_psr.name,
    ...                 parameters={'reversalPotential': (-70.0, "mV")})
    >>> connection_type = nineml.ConnectionType(
    ...                 name="Static connections",
    ...                 definition="http://svn.incf.org/svn/nineml/catalog/connectiontypes/static_connection.xml",
    ...                 parameters={'delay': (0.3, "ms")})
    
And now we can create our Projections::

    >>> exc2all = nineml.Projection(
    ...                 name="Excitatory connections",
    ...                 source=exc_cells,
    ...                 target=network,
    ...                 rule=connection_rule,
    ...                 synaptic_response=exc_psr,
    ...                 connection_type=connection_type)
    >>> inh2all = nineml.Projection(
    ...                 name="Inhibitory connections",
    ...                 source=inh_cells,
    ...                 target=network,
    ...                 rule=connection_rule,
    ...                 synaptic_response=inh_psr,
    ...                 connection_type=connection_type)
    
and add them to our Group::

    >>> network.add(exc2all, inh2all)
    
Finally, we create a container for the entire model and add the group to it::

    >>> model = nineml.Model("User-layer tutorial model")
    >>> model.add_group(network)

And now we can export it as XML::

    >>> model.write("user_layer_tutorial_model.xml")
'''
