
# --------brainx utils------------------------------------------------------
# These utils were copied over from brainx - needed for viz

def plot_channels(xdata,ydata, channel_names=None, fig=None, x_tick_rot=0,
                        size=None, cmap=plt.cm.RdBu_r, colorbar=True,
                        color_anchor=None, title=None):
    r"""Creates a lower-triangle of the matrix of an nxn set of values. This is
    the typical format to show a symmetrical bivariate quantity (such as
    correlation or coherence between two different ROIs).

    Parameters
    ----------

    in_m: nxn array with values of relationships between two sets of rois or
    channels

    channel_names (optional): list of strings with the labels to be applied to
    the channels in the input. Defaults to '0','1','2', etc.

    fig (optional): a matplotlib figure

    cmap (optional): a matplotlib colormap to be used for displaying the values
    of the connections on the graph

    title (optional): string to title the figure (can be like '$\alpha$')

    color_anchor (optional): determine the mapping from values to colormap
        if None, min and max of colormap correspond to min and max of in_m
        if 0, min and max of colormap correspond to max of abs(in_m)
        if (a,b), min and max of colormap correspond to (a,b)

    Returns
    -------

    fig: a figure object

    """
    N = in_m.shape[0]
    ind = np.arange(N)  # the evenly spaced plot indices

    def channel_formatter(x, pos=None):
        thisind = np.clip(int(x), 0, N - 1)
        return channel_names[thisind]

    if fig is None:
        fig = plt.figure()

    if size is not None:

        fig.set_figwidth(size[0])
        fig.set_figheight(size[1])

    w = fig.get_figwidth()
    h = fig.get_figheight()

    ax_im = fig.add_subplot(1, 1, 1)

    #If you want to draw the colorbar:
    if colorbar:
        divider = make_axes_locatable(ax_im)
        ax_cb = divider.new_vertical(size="10%", pad=0.1, pack_start=True)
        fig.add_axes(ax_cb)

    #Make a copy of the input, so that you don't make changes to the original
    #data provided
    m = in_m.copy()

    #Null the upper triangle, so that you don't get the redundant and the
    #diagonal values:
    #idx_null = triu_indices(m.shape[0])
    #m[idx_null] = np.nan

    #Extract the minimum and maximum values for scaling of the
    #colormap/colorbar:
    max_val = np.nanmax(m)
    min_val = np.nanmin(m)

    if color_anchor is None:
        color_min = min_val
        color_max = max_val
    elif color_anchor == 0:
        bound = max(abs(max_val), abs(min_val))
        color_min = -bound
        color_max = bound
    else:
        color_min = color_anchor[0]
        color_max = color_anchor[1]

    #The call to imshow produces the matrix plot:
    #im = ax_im.imshow(m, origin='upper', interpolation='nearest',
    #   vmin=color_min, vmax=color_max, cmap=cmap)
    im=ax_im.plot(xdata,ydata)
    #Formatting:
    ax = ax_im
    ax.grid(True)
    #Label each of the cells with the row and the column:
    if channel_names is not None:
        for i in xrange(0, m.shape[0]):
            if i < (m.shape[0] - 1):
                ax.text(i - 0.3, i, channel_names[i], rotation=x_tick_rot)
            if i > 0:
                ax.text(-1, i + 0.3, channel_names[i],
                        horizontalalignment='right')

        ax.set_axis_off()
        ax.set_xticks(np.arange(N))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(channel_formatter))
        fig.autofmt_xdate(rotation=x_tick_rot)
        ax.set_yticks(np.arange(N))
        ax.set_yticklabels(channel_names)
        ax.set_ybound([-0.5, N - 0.5])
        ax.set_xbound([-0.5, N - 1.5])

    #Make the tick-marks invisible:
    for line in ax.xaxis.get_ticklines():
        line.set_markeredgewidth(0)

    for line in ax.yaxis.get_ticklines():
        line.set_markeredgewidth(0)

    ax.set_axis_off()

    if title is not None:
        ax.set_title(title)
    """
    #The following produces the colorbar and sets the ticks
    if colorbar:
        #Set the ticks - if 0 is in the interval of values, set that, as well
        #as the maximal and minimal values:
        if min_val < 0:
            ticks = [color_min, min_val, 0, max_val, color_max]
        #Otherwise - only set the minimal and maximal value:
        else:
            ticks = [color_min, min_val, max_val, color_max]

        #This makes the colorbar:
        cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal',
                          cmap=cmap,
                          norm=im.norm,
                          boundaries=np.linspace(color_min, color_max, 256),
                          ticks=ticks,
                          format='%.2f')

    # Set the current figure active axis to be the top-one, which is the one
    # most likely to be operated on by users later on
    fig.sca(ax)
    """
    return fig
# --------brainx utils------------------------------------------------------
# These utils were copied over from brainx - needed for viz


#fig03 = draw_graph(dirfinals2)
#p1.show()
#sfin='grange_graph'+str(h.prunenet)+str(fin)+'.png' 
#fig02.savefig(sfin)
#sfin='grange_graphfinal'+str(h.prunenet)+str(fin)+'.png' 

fin+=1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(thmsgcv2)#, roi_names)
p1.title('outright SGC')
p1.draw() 

fin+=1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(thent)#, roi_names)
p1.title('outright nTE')
p1.draw() 



"""



class results(object):
    idvec0 = []
    tvec0 = []

    idvec1 = []
    tvec1 = []

    # The class "constructor" - It's actually an initializer 
    def __init__(self, name, age, major):
        self.idvec = idvec
        self.tvec = tvec

def make_student(tvec, idvec):
    r = results(idvec, tvec 2
    return r




pyra=[int(i) for i in pyra]
for j in xrange(0,max(pyra)):
rows12 = [i for i in range(0,len(pyra)) if pyra[i]==j]
if(len(rows12)!=0):
T12 = [tvec[i] for i in rows12]
T12=[int(i) for i in T12]

s=DiscreteSystem(np.array(T12), (1,max(T12)+1), np.array(T12), (1,max(T12)+1))
s.calculate_entropies(method='plugin', calc=['HX'])
print s.H



tvec2[:]=tvec[:]*1000
tvec3=[int(i) for i in tvec2]


T13 = [tvec[i] for i in rows13]
s=DiscreteSystem(np.array(tvec3), (1,max(tvec3)+1), np.array(tvec3), (1,max(tvec3)+1))


T13=[int(i) for i in T13]



for indexi in xrange(0,int(numcell-1)):
  for indexj in xrange(0,int(numcell-1)):
    #X=pe.quantise_discrete(np.array(time_courses[int(indexi)]), 2)
    #Y=pe.quantise_discrete(np.array(time_courses[int(indexj)]), 2)
    #X=trains[indexi]
    #Y=trains[indexj]
    #tvec,pyra
    X= [int(float(x)) for x in timespython[indexi]]
    Y= [int(float(x)) for x in timespython[indexj]]

    X= [int(float(x)) for x in trains[indexi]]
    Y= [int(float(x)) for x in trains[indexj]]


    #tvec[pyra[0]] # tvec[np.where(pyra==2)]
    #Y=#tvec[pyra[1]]
    
    out=np.where(np.array(tvec)!=idvec[i] np.array(tvec), 0)   
    # reshape to input, output vectors
    #http://code.google.com/p/pyentropy/wiki/SupplementalData
    #R = data.T.reshape(nt*ns,L).T.astype(int)
    #S = r_[0:ns].repeat(nt).astype(int)
    #Needs a spike time format similar to tvec. Only, it needs to be an integer version of tvec.
    #Without losing precision would involve multiplication.
    #sys = DiscreteSystem(R,(L,2),S,(1,13))

    X_dims=(len(X))
    Y_dims=(len(Y))
    s=DiscreteSystem(X, X_dims, Y, Y_dims)
    s.calculate_entropies(method='plugin', calc=['HX', 'HXY'])
    s.I()
    s.calculate_entropies(method='pt', calc=['HX', 'HXY'])
    s.I()
    s.calculate_entropies(method='qe', calc=['HX', 'HXY'])
    s.I()
    s.calculate_entropies(method='nsb', calc=['HX', 'HXY'])
    s.I()
fin=0
fig=p1.figure(fin)
fig.clf()
"""


#p1.show()
#numpy.where(condition[, x, y])
    

#p1.clf()
#nx.draw_networkx(dirfinal,nodelist=dicten,edgelist=dictee)#, roi_names)
#nx.draw(dirfinal)
#nx.draw(poot)
#fig = drawgraph_channels(dirfinal, roi_names)
#sfin='Transfer Entropy graph centrality 
#sfin='now 2'+str(h.prunenet)+str(fin)+'.png' 
#fig.savefig(sfin)

execfile('pyhoc.py')

""" 
dictc=nx.degree_centrality(dirfinals2)
#nx.draw(dirfinals2[dictc[0,10]])#, roi_names)
    
fin+=1
fig = p1.figure(fin)
p1.clf()
nx.draw_networkx(dirfinals,nodelist=dictc)#, roi_names)
sfin='subset'+str(h.prunenet)+str(fin)+'.png' 
fig.savefig(sfin)

dicten=nx.betweenness_centrality(dirfinal)
dictee=nx.edge_betweenness_centrality(dirfinal)
#nx.draw(dirfinals2[dictc[0,10]])#, roi_names)
listout=sorted(dicte.items(), key=itemgetter(1))
l2=listout[len(listout)-15:]


x= np.empty((len(vec1), vec1.shape[0]))
y = np.empty((len(vec2), vec2.shape[0]))

for i in xrange(N):
    frex, x[i], nu = alg.multi_taper_psd(vec1)
    frex, y[i], nu = alg.multi_taper_psd(vec2)
"""    
#from bsmart import timefreq, pwcausalr
#from scipy import array, size

#This is what the algorithm is actually doing.
          #Rxx = tsu.autocov_vector(np.vstack([time_courses[indexi],time_courses[indexj]]), nlags=lag)
          #coef, ecov = alg.lwr_recursion(np.array(Rxx).transpose(2, 0, 1))
          #w, f_x2y, f_y2x, f_xy, Sw = alg.granger_causality_xy(coef,
          #ecov,
          #n_freqs=n_freqs)
          
