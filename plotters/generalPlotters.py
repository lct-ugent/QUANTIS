'''
Created on 7 Oct 2016

@author: shsymoen
'''

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np


def colors_markers(no):
    ''' all unfilled markers '''
#     # Filter out filled markers and marker settings that do nothing.
#     # We use iteritems from six to make sure that we get an iterator
#     # in both python 2 and 3
#     markers = [m for m, func in iteritems(Line2D.markers)
#                         if func != 'nothing' and m not in
#                           Line2D.filled_markers]
#     # Reverse-sort for pretty. We use our own sort key which is essentially
#     # a python3 compatible reimplementation of python2 sort.
#     markers = sorted(markers,
#                               key=lambda x: (str(type(x)), str(x)))[::-1]

    markers = ['o', '^', 's', 'D', '*', '<', '>', 'v', 'p', 'd', 'H', '8', 'h']
    colors = [
        'g',
        'darkblue',
        'r',
        'deepskyblue',
        'darkmagenta',
        'coral',
        'indianred'
    ]

    return colors[:no], markers[:no]


def setAxLinesBW(ax, baw=False):
    """
    Take each Line2D in the axes, ax, and convert the line style to be
    suitable for black and white viewing.
    """
    # http://stackoverflow.com/questions/7358118/
    # matplotlib-black-white-colormap-with-dashes-dots-etc?noredirect=1&lq=1
    markersize = 3
    colormap = {
        'g': {'marker': None, 'dash': (None, None), 'color': 'g'},
        'darkblue': {'marker': None, 'dash': [5, 5], 'color': 'darkblue'},
        'r': {'marker': None, 'dash': [5, 3, 1, 3], 'color': 'r'},
        'deepskyblue': {'marker': None, 'dash': [1, 3], 'color': 'deepskyblue'},
        'coral': {'marker': None, 'dash': [5, 2, 5, 2, 5, 10], 'color': 'coral'},
        'indianred': {'marker': None, 'dash': [5, 3, 1, 2, 1, 10], 'color': 'indianred'},
        'k': {'marker': 'o', 'dash': (None, None), 'color': 'k'},
        # As of Matplotlib 2.0 the default colors have changed.
        '#1f77b4': {'marker': None, 'dash': (None, None), 'color': '#1f77b4'},
        '#ff7f0e': {'marker': None, 'dash': [5, 5], 'color': '#ff7f0e'},
        '#2ca02c': {'marker': None, 'dash': [5, 3, 1, 3], 'color': '#2ca02c'},
        '#d62728': {'marker': None, 'dash': [1, 3], 'color': '#d62728'},
        '#9467bd': {
            'marker': None,
            'dash': [5, 2, 5, 2, 5, 10],
            'color': '#9467bd'
        },
        '#8c564b': {
            'marker': None,
            'dash': [5, 3, 1, 2, 1, 10],
            'color': '#8c564b'
        },
        '#e377c2': {'marker': 'o', 'dash': (None, None), 'color': '#e377c2'}
    }

    lines_to_adjust = ax.get_lines()
    try:
        lines_to_adjust += ax.get_legend().get_lines()
    except AttributeError:
        pass

    # Adapt lines in plot
    for line in ax.get_lines():
        origColor = line.get_color()
        line_label = line.get_label()
        if line_label != "_nolegend_":
            line.set_color(colormap[origColor]['color'])
            if baw:
                line.set_color('black')
            line.set_dashes(colormap[origColor]['dash'])
            line.set_marker(colormap[origColor]['marker'])
            line.set_markersize(markersize)
        else:
            line.set_linestyle("None")

    # Adapt line in legend
    try:
        for line in ax.get_legend().get_lines():
            origColor = line.get_color()
            if baw:
                line.set_color('black')
            line.set_dashes(colormap[origColor]['dash'])
            line.set_marker(colormap[origColor]['marker'])
            line.set_markersize(markersize)
    except:
        print("No lines could be changed.")


def setFigLinesBW(fig, baw=False):
    """
    Take each axes in the figure, and for each line in the axes, make the
    line viewable in black and white.
    """
    # http://stackoverflow.com/questions/7358118/
    # matplotlib-black-white-colormap-with-dashes-dots-etc?noredirect=1&lq=1
    for ax in fig.get_axes():
        setAxLinesBW(ax, baw)


def colorbar_index(ncolors, cmap, tick_labels, cbar_label):
    """ Modifies cmap and corresponding colorbar """
    # http://stackoverflow.com/questions/18704353/
    # correcting-matplotlib-colorbar-ticks
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors+0.5)
    colorbar = plt.colorbar(mappable)
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
    colorbar.set_ticklabels(tick_labels)
    colorbar.ax.set_ylabel(cbar_label)


def create_discrete_colormap_column(column, colormapName):
    """Create a discrete colormap for a given column with repetitive values
    """
    dict_labels = {cell for cell in column}
    number_of_labels = len(dict_labels)
    cmap = create_discrete_colormap(number_of_labels, colormapName)
    cmap_labels = list(dict_labels)
    cmap_labels.sort()

    return cmap, cmap_labels


def create_discrete_colormap(number, colormapName):
    """ Creates a discrete colormap for a given number """
    import matplotlib.pyplot as plt
    cmap = plt.cm.get_cmap(colormapName, number)
    return cmap


def f_save(f, f_name, dpi=300):
    '''Saves a matplotlib figure to a given location

    Parameters
    ----------
    f : matplotlib.figure.Figure object
    f_name : string
        Name of the figure to be created.
        The string can include the folder location.
    dpi : int, optional (default: 300)
        Dots per inch

    Returns
    -------
    None
    '''
    f.savefig(
        f_name,
        format='tif',
        dpi=dpi,
    )
    return None

