#!/usr/bin/env python3


#-------------------------------------------
# All matplotlib plotting related stuff
#-------------------------------------------

import matplotlib as mpl


def setplotparams_single_plot(left=0.07, right=0.97, bottom=0.10, top=0.92, wspace=0.15, hspace=0.12):
    """
    Set rcParams for nice looking plots.
    Optimized for a figure with a single subplot.
    """
    # Plot parameters
    params = {'axes.labelsize': 10,
    'axes.titlesize': 10,
    'font.size': 12,
    'font.family': 'serif', 
    'font.serif': 'Comuter Modern',
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.subplot.left'    : left,
    'figure.subplot.right'   : right,
    'figure.subplot.bottom'  : bottom,
    'figure.subplot.top'     : top,
    'figure.subplot.wspace'  : wspace,
    'figure.subplot.hspace'  : hspace,
    'figure.figsize' : (5, 5.5), 
    'figure.dpi' : 300,
    'lines.markersize' : 6,
    'lines.linewidth' : 3.
    }
    mpl.rcParams.update(params)

    return






def setplotparams_multiple_plots(left=0.05, right=0.97, bottom=0.10, top=0.92, wspace=0.15, hspace=0.15):
    """
    Set rcParams for nice looking plots. Intended for about figsize 5 x 5.5 per subplot. Set those when
    calling plt.figure()
    """
    # Plot parameters
    params = {'axes.labelsize': 10,
    'axes.titlesize': 12,
    'font.size': 12,
    'font.family': 'serif', 
    'font.serif': 'Comuter Modern',
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.subplot.left'    : left,
    'figure.subplot.right'   : right,
    'figure.subplot.bottom'  : bottom,
    'figure.subplot.top'     : top,
    'figure.subplot.wspace'  : wspace,
    'figure.subplot.hspace'  : hspace,
    'figure.dpi' : 300,
    'lines.markersize' : 6,
    'lines.linewidth' : 3.
    }
    mpl.rcParams.update(params)

    return


