
from mpl_toolkits.axes_grid import make_axes_locatable
import matplotlib.ticker as ticker


""" --------brainx utils------------------------------------------------------
 These utils were copied over from brainx - needed for viz
Copyright (c) 2006-2009, NIPY Developers
All rights reserved.

https://github.com/fperez/brainx/blob/master/LICENSE

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

    * Neither the name of the NIPY Developers nor the names of any
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
     "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
     LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
     A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
     OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
     SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
     LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
     DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
     THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
     (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


def mutual_information(d1, d2):
    """Mutual information between two graph partitions.<-In Brain structure.

Read in two dictionaries of sets (i.e. a graph partition) and assess how
similar they are using mutual information as in Danon, Diaz-Guilera, Duch &
Arenas, J Statistical Mechanics 2005.

Parameters
----------
d1 : dict
dictionary of 'real communities'
d2 : dict
dictionary of 'found communities'

Returns
-------
mi : float
Value of mutual information between the two partitions.
"""
    log = np.log
    nansum = np.nansum

    N = confusion_matrix(d1, d2)

    nsum_row = N.sum(0)[np.newaxis, :]
    nsum_col = N.sum(1)[:, np.newaxis]

    # Sanity checks: a zero in either of these can only happen if there was an
    # empty module in one of the input partitions. Rather than manually check
    # the entire partitions, we look for this problem at this stage, and bail
    # if there was an empty module.
## if (nsum_row==0).any():
## EmptyModuleError("Empty module in second partition.")
## if (nsum_col==0).any():
## EmptyModuleError("Empty module in first partition.")

    # nn is the total number of nodes
    nn = nsum_row.sum()
    num = nansum(N*log(N*nn/(nsum_row*nsum_col)))
    den = nansum(nsum_row*log(nsum_row/nn)) + nansum(nsum_col*log(nsum_col/nn))

    return -2*num/den



def find_unconnected_nodes(self):
    """ checks for nodes in graph with no edges """
    graph = nx.from_numpy_matrix(self.graph_adj_matrix)
    unconnected = [ n for n,d in graph.degree_iter() if d==0 ]
    return unconnected

def drawmatrix_channels(
    in_m,
    channel_names=None,
    fig=None,
    x_tick_rot=0,
    size=None,
    cmap=p1.cm.RdBu_r,
    colorbar=True,
    color_anchor=None,
    title=None,
    ):
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
        fig = p1.figure()

    if size is not None:

        fig.set_figwidth(size[0])
        fig.set_figheight(size[1])

    w = fig.get_figwidth()
    h = fig.get_figheight()

    ax_im = fig.add_subplot(1, 1, 1)

    # If you want to draw the colorbar:

    if colorbar:
        divider = make_axes_locatable(ax_im)
        ax_cb = divider.new_vertical(size='10%', pad=0.1,
                pack_start=True)
        fig.add_axes(ax_cb)

    # Make a copy of the input, so that you don't make changes to the original
    # data provided

    m = in_m.copy()

    # Null the upper triangle, so that you don't get the redundant and the
    # diagonal values:
    # idx_null = triu_indices(m.shape[0])
    # m[idx_null] = np.nan

    # Extract the minimum and maximum values for scaling of the
    # colormap/colorbar:

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

    # The call to imshow produces the matrix plot:

    im = ax_im.imshow(
        m,
        origin='upper',
        interpolation='nearest',
        vmin=color_min,
        vmax=color_max,
        cmap=cmap,
        )

    # Formatting:

    ax = ax_im
    ax.grid(True)

    # Label each of the cells with the row and the column:

    if channel_names is not None:
        for i in xrange(0, m.shape[0]):

            # if i < (m.shape[0] - 1):
            #    ax.text(i - 0.3, i, channel_names[i], rotation=x_tick_rot)
            # I suspect this is the cause of the x axis ticks.

            if i > 0:
                ax.text(-1, i + 0.3, channel_names[i],
                        horizontalalignment='right')

        ax.set_axis_off()
        ax.set_xticks(np.arange(N))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(channel_formatter))
        fig.autofmt_xdate(rotation=x_tick_rot)
        ax.set_yticks(np.arange(N))
        ax.set_yticklabels(channel_names)
        ax.set_ybound([-0.50, N - 0.50])
        ax.set_xbound([-0.50, N - 1.5])

    # Make the tick-marks invisible:

    for line in ax.xaxis.get_ticklines():
        line.set_markeredgewidth(0)

    for line in ax.yaxis.get_ticklines():
        line.set_markeredgewidth(0)

    ax.set_axis_off()

    if title is not None:
        ax.set_title(title)

    # The following produces the colorbar and sets the ticks

    if colorbar:

        # Set the ticks - if 0 is in the interval of values, set that, as well
        # as the maximal and minimal values:

        if min_val < 0:
            ticks = [color_min, min_val, 0, max_val, color_max]
        else:

        # Otherwise - only set the minimal and maximal value:

            ticks = [color_min, min_val, max_val, color_max]

        # This makes the colorbar:

        cb = fig.colorbar(
            im,
            cax=ax_cb,
            orientation='horizontal',
            cmap=cmap,
            norm=im.norm,
            boundaries=np.linspace(color_min, color_max, 256),
            ticks=ticks,
            format='%.2f',
            )

    # Set the current figure active axis to be the top-one, which is the one
    # most likely to be operated on by users later on

    fig.sca(ax)

    return fig


def drawmatrix_channels_nullh(
    in_m,
    channel_names=None,
    fig=None,
    x_tick_rot=0,
    size=None,
    cmap=p1.cm.RdBu_r,
    colorbar=True,
    color_anchor=None,
    title=None,
    ):
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
        fig = p1.figure()

    if size is not None:

        fig.set_figwidth(size[0])
        fig.set_figheight(size[1])

    w = fig.get_figwidth()
    h = fig.get_figheight()

    ax_im = fig.add_subplot(1, 1, 1)

    # If you want to draw the colorbar:

    if colorbar:
        divider = make_axes_locatable(ax_im)
        ax_cb = divider.new_vertical(size='10%', pad=0.1,
                pack_start=True)
        fig.add_axes(ax_cb)

    # Make a copy of the input, so that you don't make changes to the original
    # data provided

    m = in_m.copy()

    # Null the upper triangle, so that you don't get the redundant and the
    # diagonal values:
    idx_null = triu_indices(m.shape[0])
    m[idx_null] = np.nan

    # Extract the minimum and maximum values for scaling of the
    # colormap/colorbar:

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

    # The call to imshow produces the matrix plot:

    im = ax_im.imshow(
        m,
        origin='upper',
        interpolation='nearest',
        vmin=color_min,
        vmax=color_max,
        cmap=cmap,
        )

    # Formatting:

    ax = ax_im
    ax.grid(True)

    # Label each of the cells with the row and the column:

    if channel_names is not None:
        for i in xrange(0, m.shape[0]):

            if i < (m.shape[0] - 1):
                ax.text(i - 0.3, -0.3, channel_names[i], 
                        horizontalalignment='right')
               # I suspect this is the cause of the x axis ticks.

            if i > 0:
                ax.text(-1, i + 0.3, channel_names[i],
                        horizontalalignment='right')

        ax.set_axis_off()
        ax.set_xticks(np.arange(N))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(channel_formatter))
        fig.autofmt_xdate(rotation=x_tick_rot)
        ax.set_yticks(np.arange(N))
        ax.set_yticklabels(channel_names)
        ax.set_ybound([-0.50, N - 0.50])
        ax.set_xbound([-0.50, N - 1.5])

    # Make the tick-marks invisible:

    for line in ax.xaxis.get_ticklines():
        line.set_markeredgewidth(0)

    for line in ax.yaxis.get_ticklines():
        line.set_markeredgewidth(0)

    ax.set_axis_off()

    if title is not None:
        ax.set_title(title)

    # The following produces the colorbar and sets the ticks

    if colorbar:

        # Set the ticks - if 0 is in the interval of values, set that, as well
        # as the maximal and minimal values:

        if min_val < 0:
            ticks = [color_min, min_val, 0, max_val, color_max]
        else:

        # Otherwise - only set the minimal and maximal value:

            ticks = [color_min, min_val, max_val, color_max]

        # This makes the colorbar:

        cb = fig.colorbar(
            im,
            cax=ax_cb,
            orientation='horizontal',
            cmap=cmap,
            norm=im.norm,
            boundaries=np.linspace(color_min, color_max, 256),
            ticks=ticks,
            format='%.2f',
            )

    # Set the current figure active axis to be the top-one, which is the one
    # most likely to be operated on by users later on

    fig.sca(ax)

    return fig


# --------brainx utils------------------------------------------------------
# These utils were copied over from brainx - needed for viz

fig=p1.figure(fin)

fig03 = drawmatrix_channels(np.matrix(Matrix[:len(Matrix)-1]), roi_names[:len(Matrix)-1],
                            size=[10., 10.], color_anchor=0)

# ROI names and the matrix need to be the same length.
# I think I am justified using

fig03.savefig(str(h.plastic) + str(h.ff) + str(tstop) + str(ncell)
              + str(h.prunenet) + str(fin) + 'Degree_Matrix.png')

