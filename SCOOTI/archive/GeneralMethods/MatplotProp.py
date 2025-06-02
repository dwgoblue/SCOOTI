##########################Packages############################

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import seaborn as sns
import scipy.stats as ss

##############################################################
def CanvasStyle(ax, square=False, remove_bottom=False, remove_left=False, lw=3, ticks_lw=2):
    ax.patch.set_facecolor('white')
    ax.grid('off')
    plt.tight_layout()
    ax.grid(False)
    ax.yaxis.set_tick_params(width=ticks_lw)
    ax.xaxis.set_tick_params(width=ticks_lw)
    if square==False:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        if remove_bottom==False:
            ax.axhline(linewidth=lw, y=ax.get_ylim()[0], color='k')
        if remove_left==False:
            ax.axvline(linewidth=lw, x=ax.get_xlim()[0], color='k')

    else:
        ax.axhline(linewidth=lw, y=ax.get_ylim()[0], color='k')
        ax.axvline(linewidth=lw, x=ax.get_xlim()[0], color='k')
        ax.axhline(linewidth=lw, y=ax.get_ylim()[1], color='k')
        ax.axvline(linewidth=lw, x=ax.get_xlim()[1], color='k')

    return ax

#############################################################

# label p value and significance.
def Significance(arrs, ax, columns=[], mode='box', test_func=ss.ttest_ind):
    
    def stars(p):
        if p < 0.0001:
            return "****"
        elif (p < 0.001):
            return "***"
        elif (p < 0.01):
            return "**"
        elif (p < 0.05):
            return "*"
        else:
            return "n.s."
    
    def ttest(arr1, arr2):
        # Calculate t-test and p value.
        # Use scipy.stats.ttest_ind.
        t, p = test_func(arr1, arr2)
        s = stars(p)
        print("ttest_ind:            t = %g  p = %g" % (t, p))
        return s
    
    trans = ax.get_xaxis_transform()
    label_min, label_max = ax.get_ylim()
    props = {'connectionstyle':"bar,fraction={0}".format(0.1),
             'arrowstyle':'-',
             'linewidth':2,
             'ec':'#000000'}
    rank = np.array([abs(v[0]-v[1]) for v in columns])
    cols = np.array([columns[ind] for ind in np.argsort(rank)])
    overlap_record = np.zeros(len(arrs))
    if mode == 'bar':
        x_position = np.arange(len(arrs))/2+0.25
        y_max = []
        for i in range(len(arrs)):
            y = np.mean(arrs[i])+np.std(arrs[i])
            y_max.append(y)
        y_max = max(y_max)
        standard = 0.025*(len(arrs)-1)
        
    else:
        x_position = np.arange(len(arrs))+0.75+0.25
        y_max = []
        DataMax = []
        for i in range(len(arrs)):
            DataMax.append(np.percentile(arrs[i], 98))
            Q1 = np.percentile(arrs[i], 75)
            Q3 = np.percentile(arrs[i], 25)
            max_bound = Q1+1.5*(Q1-Q3)
            y = np.sort(arrs[i], axis=None)[np.nonzero(np.sort(arrs[i], axis=None)<=max_bound)[0][-1]]
            y_max.append(y)
        if mode == "box":
            y_max = max(y_max)
            standard = 0.05*(len(arrs)-1)
        else:
            y_max = max(y_max)
            standard = 0.05*2
        
    for col in cols:
        # Calculate t-test and p value.
        s = ttest(arrs[col[0]], arrs[col[1]])
        # Update props.
        props['connectionstyle'] = "bar,fraction={0}".format(
            standard/(abs(x_position[col[0]]-x_position[col[1]])))
        # Process the label.
        passby = np.arange(col[0], col[1]+1)
        adj = (y_max-label_min)/(label_max-label_min)
        adj = adj+0.04*(1+np.max(overlap_record[passby]))
        overlap_record[passby] +=1
        ax.annotate("", xy=(x_position[col[0]], adj),
                    xycoords=trans,
                    xytext=(x_position[col[1]], adj),
                    textcoords=trans,
                    arrowprops=props)
        if s=="n.s.":
            text_adj = y_max+(label_max-label_min)*0.04*(
            1+np.max(overlap_record[passby]))+0.005
        else:
            text_adj = y_max+(label_max-label_min)*0.04*(
                1+np.max(overlap_record[passby]))-(label_max-label_min)*0.01
        ax.text((x_position[col[0]]+x_position[col[1]])/2,
                text_adj,
                s,
                horizontalalignment='center',
                verticalalignment='center',
                weight='bold',)
#               backgroundcolor='white')
    return ax

#########################################################

def PltProps():
    plt.rcParams['font.weight'] = 'normal'#'bold'
    plt.rcParams['font.size'] = 20
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['axes.labelsize'] = 20
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams["lines.linewidth"] = 3
    plt.rcParams["axes.titlesize"] = 42
    plt.rcParams['axes.linewidth'] = 3
    #plt.rcParams['axes.spines.left'] = True
    #plt.rcParams['axes3d.xaxis.panecolor'] = (0.95, 0.95, 0.95, 0.5)  # background pane on 3D axes
    #plt.rcParams['axes3d.yaxis.panecolor'] = (0.90, 0.90, 0.90, 0.5)  # background pane on 3D axes
    #plt.rcParams['axes3d.zaxis.panecolor'] = (0.925, 0.925, 0.925, 0.5)
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    errbarattr = {'lw':5,
                       'capthick':3,
                       'capsize':10,
                       'ecolor':'black',
                       'markersize':15,
                       'marker':'o'
                      }

def SeabornProp():
    style = {'axes.facecolor': 'white',
 'axes.edgecolor': '.8',
 'axes.grid': False,
 'axes.axisbelow': True,
 'axes.labelcolor': '.15',
 'figure.facecolor': 'white',
 'grid.color': '.8',
 'grid.linestyle': '-',
 'text.color': '.15',
 'xtick.color': '.15',
 'ytick.color': '.15',
 'xtick.direction': 'out',
 'ytick.direction': 'out',
 'lines.solid_capstyle': 'round',
 'patch.edgecolor': 'w',
 'image.cmap': 'rocket',
 'font.family': ['sans-serif'],
 'font.sans-serif': ['Arial',
  'DejaVu Sans',
  'Liberation Sans',
  'Bitstream Vera Sans',
  'sans-serif'],
 'patch.force_edgecolor': True,
 'xtick.bottom': False,
 'xtick.top': False,
 'ytick.left': False,
 'ytick.right': False,
 'axes.spines.left': False,
 'axes.spines.bottom': True,
 'axes.spines.right': False,
 'axes.spines.top': True,}
    return style
##############################################################
