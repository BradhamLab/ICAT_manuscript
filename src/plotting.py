import itertools
import json
import os
import sys

import colorcet as cc
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import ternary
from cycler import cycler
from scanpy import api as sc

from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D

from downstream.src.visualization import visualize



loc = os.path.dirname(os.path.abspath(__file__))
plt.style.use(os.path.join(loc, 'configs/figures.mplstyle'))
# plt.rc('axes', prop_cycle=cycler('color', cc.glasbey_light))

method_dictionary = {
    'icat': 'ICAT',
    'seurat311': 'Seurat 3.1',
    'scanorama': 'scanorama',
    'icat_scan': 'scanorama + ICAT',
    'seurat_icat': 'Seurat 3.1 + ICAT',
    'no-int': 'No Int.',
    'ncfs-louvain': 'NCFS + Louvain'
}

# metric_dictionary = {
#     'adjusted.mutual.info': 'AMI',
#     'adjusted.rand': 'ARI',
#     'completeness': 'Completeness',
#     'fowlkes.mallows': 'Fowlkes-Mallows',
#     'homogeneity': 'Homogeneity',
#     'kbet': 'kBET',
#     'kBET': 'kBET'
# }

metric_dictionary = {
    'adjusted.mutual.info': 'AMI',
    'adjusted.rand': 'ARI',
    'completeness': 'Comp.',
    'fowlkes.mallows': 'FM',
    'homogeneity': 'Homog.',
    'silhouette': 'Sil.',
    'calinski': 'CH',
    'davies': 'DB'
}

method_colors = ["#0a5e62", "#4bbfb8", "#fbcf5b",
                 "#ff5959", "#4f8522", "#a2d665"]

method_colors1 = ["#3c4347", "#5e6769", "#719192",
                  "#e1cec5", "#518413", "#a3c541"]

method_colors2 = ["#f77855", "#584b43", "#537d90",
                  "#a5d2c9", "#518413", "#a3c541"]

method_colors3 = ["#91919f", "#e6412b", "#f2b420",
                  "#181a1c", "#32AE65", "#1E6C7E", "#E48582"]

possible = ["#32AE65", "#11664D", "#0a5e62"]

palette = ["#e6412b", "#1E6C7E", "#32AE65", "#f2b420",
           "#E48582", "#181a1c", "#91919f"]

cyndi_palette = ['#FB4D3D', '#F3752B', '#32AE65', '#054A29',
                 '#181a1c', '#A846A0', '#72A1E5', '#345995',
                 '#242038', '#E8AEB7', '#513C2C', '#3A3335',
                 '#7E6551', '#593C8F', '#4C2A85', '#285943',
                 '#4A5043', '#596157', '#E3655B', '#0f783e',
                 '#0f5e3e']

cyndi_palette = ['#FB4D3D', '#F3752B', '#A846A0',
                 '#32AE65', '#0F5E3E',
                 '#72A1E5', '#345995']

dakota_palette = ['#FB4D3D', '#E3655B', '#596157', 
                  '#32AE65', '#0F5E3E',
                  '#72A1E5', '#345995']


            # icat   ncfs-louvain  no int
palette = ['#FB4D3D', '#F3752B', '#A846A0',
#           seurat    seurat-icat
           '#32AE65', '#0F5E3E',
#           scanorama  icat-scan
           '#72A1E5', '#345995']

    
                 # icat   ncfs-louvain  no int
palette_again = ['#d62728', '#F3752B', '#A846A0',
      #           seurat    seurat-icat
                 '#61D190', '#0f783e',
      #           scanorama  icat-scan
                 '#72A1E5', '#273f87 '
]



def radar_factory(num_vars, frame='circle'):
    """
    Taken from mpl examples <insert link>

    Create a radar chart with `num_vars` axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle', 'polygon'}
        Shape of frame surrounding axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    class RadarAxes(PolarAxes):

        name = 'radar'
        # use 1 line segment to connect specified points
        RESOLUTION = 1

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.append(x, x[0])
                y = np.append(y, y[0])
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == 'circle':
                return Circle((0.5, 0.5), 0.5)
            elif frame == 'polygon':
                return RegularPolygon((0.5, 0.5), num_vars,
                                      radius=.5, edgecolor="k")
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            elif frame == 'polygon':
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(axes=self,
                              spine_type='circle',
                              path=Path.unit_regular_polygon(num_vars))
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                    + self.transAxes)
                return {'polar': spine}
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta


def radar_plot(data, feature, color='blue', frame='polygon', ax = None):
    theta = radar_factory(data.shape[0], frame=frame)
    spoke_labels = data.index.values
    if ax is None:
        fig, ax = plt.subplots(subplot_kw=dict(projection='radar'))
        ax.set_rgrids(np.arange(0, 1.2, 0.2))
        ax.set_title(feature, weight='bold', position=(0.5, 1.1),
                     horizontalalignment='center', verticalalignment='center')
        ax.set_ylim(0, 1) 

    ax.plot(theta, data[feature], color=color, label=feature)
    ax.fill(theta, data[feature], facecolor=color, alpha=0.25, color=color)
    ax.set_varlabels(spoke_labels)
    return ax

def method_order(data):
    return data.loc[['ICAT', 'No Int.',
                     'Seurat 3.1', 'Seurat 3.1 + ICAT', 'scanorama',
                     'scanorama + ICAT'], :]


def method_radar_plot(results, metrics):
    fig, axs = plt.subplots(figsize=(21, 13), subplot_kw=dict(projection='radar'),
                            ncols=4, nrows=2)
    i = 0 
    for method in results.index.values:
        ax = axs.flat[i]
        if i == 3:
            ax.set_visible(False)
            i += 1
            ax = axs.flat[i]
        ax.set_rgrids(np.arange(0, 1.2, 0.2), fontsize='18', angle=360/10)
        ax.set_title(method, weight='bold', fontsize='large', position=(0.5, 1.15),
                     horizontalalignment='center', verticalalignment='center')
        ax.set_ylim(0, 1)
        radar_plot(results[metrics].T, method, color=colors[method], frame='polygon', ax=ax)
        i += 1
        plt.tight_layout()
    fig.subplots_adjust(wspace=0.25, hspace=0.10, top=0.85, bottom=0.05)
    plt.tight_layout()
    return fig, axs


def lollipop_plot(data, metrics, colors=palette):
    subset = data[metrics].copy()
    subset['Method'] = subset.index.values
    plotdata = subset.melt(id_vars='Method', var_name='Metric', value_name='Score')
    method_order = ['ICAT', 'NCFS + Louvain', 'No Int.',
                    'Seurat 3.1', 'Seurat 3.1 + ICAT', 'scanorama',
                    'scanorama + ICAT']
    plotdata['Method'] = pd.Categorical(plotdata['Method'], method_order)
    xmins = plotdata.groupby('Metric')['Score'].min()
    xmaxs = plotdata.groupby('Metric')['Score'].max()
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.hlines(y=range(xmins.shape[0]), xmin=xmins,
              xmax=xmaxs, color='grey', alpha=1)
    ax = sns.stripplot(data=plotdata.sort_values('Metric'),
                       y='Metric', x='Score', hue='Method',orient='h', s=16,
                       alpha=0.85, jitter=False, palette=colors, ax=ax)
    plt.grid(False, axis='y')
    plt.grid(False, axis='x')
    ax.legend(loc='upper left', bbox_to_anchor=(-0.06, -0.1), markerscale=2, fontsize='medium',
              ncol=4, labelspacing=0, columnspacing=0, handletextpad=0)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    plt.tight_layout()
    return ax

def stacked_barplot(df, label, cluster, xlabel=''):
    """
    Plot a stacked barplot.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing two columns of sample labels.
    label : str
        Name of column in `df` containing known labels.
    cluster : str
        Name of column in `df` containing predicted labels.
    xlabel : str
        Description of group plotted along the x-axis.
    """
    counts = df.groupby([cluster, label]).size().unstack().fillna(0)
    cluster_totals = counts.sum(axis=0)
    # percentages of cells belonging to each cluster for each known label 
    # (e.g. 87% of cells in known label X are in cluster Y)
    percentages = (counts / cluster_totals).T * 100
    labels = percentages.index.values
    clusters = sorted(percentages.columns.values)
    xticks = range(percentages.shape[0])
    totals = np.zeros(percentages.shape[0])
    colors = cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color'])
    # plot barplots up to calculated percentage for each cluster in known labels
    for each, color in zip(clusters, colors()):
        new_percentages = percentages.loc[:, each].values
        plt.bar(xticks, new_percentages, color=color['color'],
                width=0.85, bottom=totals, label=each)
        # update cumulate percentage for starting points
        totals += new_percentages
    legend_cols = int(np.ceil(len(clusters) / 20))
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
               ncol=legend_cols)
    plt.xticks(xticks, labels, rotation=90)
    plt.xlabel(xlabel, fontsize=24, labelpad=5)
    plt.ylabel("Percentage of Cells", fontsize=24)
    yticks, __ = plt.yticks()
    ylabels = ["{:d}%".format(int(y)) for y in yticks]
    plt.yticks(yticks, ylabels)
    plt.ylim(0, 100)
    plt.tight_layout()

def ranked_heatmap(ranks, cmap='viridis'):
    """
    Plot a heatmap of method ranks across various method.

    Low values are considered good (e.g. 1 is better than 2). Negative values
    are used during plotting for visual purposes.
    
    Parameters
    ----------
    ranks : pd.DataFrame
        A dataframe containing methods as index values and metrics as columns.
        Each cell represents the rank of the specified method in the respective
        metric compared to other methods.
    
    Returns
    -------
    mpl.Figure
        Heatmap summarizing method performance.
    """
    figsize = np.array(plt.rcParams['figure.figsize']) * 1.25
    labelsize = plt.rcParams['axes.titlesize']
    fig, ax = plt.subplots(figsize=figsize)
    plot_ranks = -1 * ranks.loc[ranks.mean(axis=1).sort_values().index, :]
    sns.heatmap(plot_ranks, cmap=cmap,
                cbar_kws={'shrink': 0.95})
    ax.set_yticklabels([x.get_text()\
                        for x in ax.get_yticklabels()],
                        fontsize=labelsize*0.75, rotation=0)
    ax.set_xticklabels([x.get_text()\
                        for x in ax.get_xticklabels()], fontsize=labelsize*0.75,
                        rotation=90)
    plt.ylim(plot_ranks.shape[0], 0)
    fig.axes[1].annotate('Worse', (-1, -0.05), xycoords='axes fraction',
                         fontsize=labelsize*0.7)
    fig.axes[1].annotate('Better', (-1, 1.01), xycoords='axes fraction',
                         fontsize=labelsize*0.7)
    fig.axes[1].set_yticklabels([''] * len(fig.axes[1].get_yticks()) )
    plt.ylabel('')
    plt.tight_layout()
    return fig


def trendplot(results, x, y, hue=None, xlabel=None):
    """
    Plot performance trends across some variable.
    
    Parameters
    ----------
    results : pd.DataFrame
        Dataframe containing performance measures and other features across
        datasets and methods.
    x : string
        Column in `results` representing the independent variable in cluster
        performance.
    y : string
        Column in `results` measuring cluster performance -- presumably affected
        by `x`.
    hue : string, optional
        Column in `results` to separate on. Default is None, and all points rows
        will be used in a single plot.
    xlabel : string, optional
        Label for x-axis, by default None, and `x` will be used.
    
    Returns
    -------
    mpl.Axes
        Plot of linear trends between `x` and `y` separated on `hue`.
    """
    lm = sns.lmplot(x=x, y=y, hue=hue, col=hue, col_wrap=3, data=results)
    xmin = min(results[x])
    xmin -= 0.1*xmin
    xmax = max(results[x])
    xmax += 0.1*xmax
    labelsize= plt.rcParams['axes.titlesize']
    plt.ylim(-0.05, 1.05)
    for ax in lm.fig.axes:
        ax.set_title(ax.get_title().replace("{} = ".format(hue), ""),
                     fontsize=int(labelsize * 0.75))
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xlim(xmin, xmax)
    if xlabel is None:
        xlabel = x
    lm.fig.text(0.5, 0.03, xlabel, ha='center', va='center',
                fontsize=labelsize)
    lm.fig.text(0.02, 0.5, y, ha='center', va='center', rotation='vertical',
                fontsize=labelsize)
    return lm


def flip(items, ncol):
    """
    Flips matplotlib legends to increment by rows before columns.

    Taken from here:
    https://stackoverflow.com/questions/10101141/matplotlib-legend-add-items-across-columns-instead-of-down
    
    Parameters
    ----------
    items : [type]
        [description]
    ncol : [type]
        [description]
    
    Returns
    -------
    [type]
        [description]
    """
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

def metric_plot(scores, errors=None, rename=True, bottom_legend=True):
    """
    Plot performance metrics across methods. 
    
    Parameters
    ----------
    scores : pd.DataFrame
        Performance data frame where each column is a different performance
        metric, and each row represents a different method.
    
    errors : pd.DataFrame, optional
        Errors associated with each method-performance pair in `scores`. Should
        be the same shape/size as `scores`.
    
    rename : bool, optional
        Whether to rename columns + rows. 
    
    bottom : bool, optional
        Whether to plot legend below x-axis.
    """
    if errors is not None:
        if not isinstance(errors, pd.DataFrame):
            raise ValueError("Expected DataFrame for `errors` parameter")
        if np.all(errors.index != scores.index):
            raise ValueError("`scores` and `errors` dataframes should have the"
                             " same index.")
        if np.all(errors.columns != scores.columns):
            raise ValueError("`scores` and `errors` dataframes should have the"
                             " same columns.")
    if rename:
        scores.rename(index=metric_dictionary, columns=method_dictionary,
                      inplace=True)
    fig, ax = plt.subplots(constrained_layout=False)
    indices = np.arange(scores.shape[0]).astype(float)
    width = np.min(np.diff(indices)) / scores.shape[1]
    indices = indices + indices * width
    colors = cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color'])
    starts = indices - width * scores.shape[1] / 2
    for method, color in zip(sorted(scores.columns), colors()):
        yerr = None
        if errors is not None:
            yerr = errors[method]
        plt.bar(starts + width, scores[method], width, color=color['color'],
                label=method,
                yerr=yerr)
        starts = starts + width
    indices = indices + width / 2
    rotation = 0
    # if not flip:
    #     rotation = 15
    plt.xticks(indices, labels=scores.index.values, rotation=rotation)
    # plt.title("Method Performance", loc='left')
    ax = plt.gca()
    if bottom_legend:
        plot_bottom_legend(ax)
    else:
        plot_right_legend(ax)
    plt.ylim(0, 1)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout()

def plot_bottom_legend(ax):
    handles, labels = ax.get_legend_handles_labels()
    legend_cols = int(np.ceil(len(labels) / 3))
    plt.legend(flip(handles, legend_cols), flip(labels, legend_cols),
               loc='upper center', bbox_to_anchor=(0.5, -0.05),
               ncol=legend_cols, frameon=False, fancybox=False)

def plot_right_legend(ax):
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1.05), frameon=False)

def close_plot():
    plt.cla()
    plt.clf() 
    plt.close()


def color_point(point, scale=9):
    """
    Transform a ternary point into an associated RGBV value.
    
    Parameters
    ----------
    point : list-like, numpy.ndarray
        A ternary point. Expected to have length equal to three.
    scale : int, optional
        Total value points sum to, by default 9.
    
    Returns
    -------
    np.ndarray
        Color associated with the provided point in RGBV space.
    """
    point = np.array(point)
    assert point.size == 3
    scaled = (point + 3) / scale
    scaled[scaled < 0] = 0
    cmaps = [plt.get_cmap('Blues'), plt.get_cmap("Reds"), plt.get_cmap('Greens')]
    rgba = np.array([cmap(p) for cmap, p in zip(cmaps, scaled)])
    out = (rgba.mean(axis=0)) * np.hstack((scaled, [1])) 
    out[-1] = 1
    return out


# [9, 0, 0] 'H1975', -> bottom right 
# [0, 9, 0] 'H2228', -> top
# [0, 0, 9] 'HCC827' -> bottom left
def ternary_plot(data, column, label, scale=9):
    """
    Summarize `benchmark` clusters by plotting cluster medians
    
    Parameters
    ----------
    data : pd.DataFrame
        Dateframe containing assigned cluster labels along with known cell
        mixtures.
    column : string
        Column in `data` containing assigned cluster labels.
    label : string
        Column in `data` containing known cell mixtures.
    scale : int, optional
        Total number of cells used in each mixture, by default 9.
    
    Returns
    -------
    mpl.Axes
        Ternary plot of cluster medians of cell mixtures.
    """
    data['c1'] = data.apply(lambda x: int(x[label].split(',')[0]), axis=1)
    data['c2'] = data.apply(lambda x: int(x[label].split(',')[1]), axis=1)
    data['c3'] = data.apply(lambda x: int(x[label].split(',')[2]), axis=1)
    by_column = data.groupby(column)
    plot_data = data.groupby(column)['c1', 'c2', 'c3'].median()
    size = by_column.size()
    plot_data['size'] = ((950 - 175) * (size - 45) /\
                         (950 - 45) + 175).astype(int)
    __, tax = ternary.figure(scale=scale)
    colors = plot_data.apply(lambda x: color_point(np.array([x.c1, x.c2, x.c3]),
                                                   scale),
                             axis=1).values
    labelsize=plt.rcParams['axes.titlesize']
    tax.gridlines(multiple=1, linewidth=1.25,
                  horizontal_kwargs={'color': color_point([0, 9, 0], scale)},
                  right_kwargs={'color': color_point([0, 0, 9], scale)},
                  left_kwargs={'color': color_point([9, 0, 0], scale)})
    tax.boundary(scale=scale, linewidth=1.25)
    tax.scatter(plot_data[['c1', 'c2', 'c3']].values, c=colors,
                s=plot_data['size'], alpha=1)
    tax.left_corner_label('HCC827', position=[-0.1, 0.075, 0],
                          fontsize=labelsize)
    tax.top_corner_label('H2228', position=[-0.02, 1.17, 0],
                         fontsize=labelsize)
    tax.right_corner_label('H1975', position=[1.02, 0.02, 0],
                           fontsize=labelsize)
    tax.ticks(axis='lbr', multiple=1, offset=0.025,
              fontsize=int(labelsize * 0.75), linewidth=1.25)
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    return tax

def xy_mean_plot(data, x, y, color, error_kwargs={}):
    """
    Plot means of `x` and `y` as a scatter plot with error bars. 

    Parameters
    ----------
    data : pd.DataFrame
        Dataframe containing data to plot.
    x : str
        Column name in `data` to plot along x-axis.
    y : str
        Column name in `data` to plot along y-axis.
    color : str
        Column name in `data` to color dots by.
    error_kwargs : dict, optional
        Key-word arguments to pass to `plt.error_bars`. Optional, deafult is an
        empty dictionary and default parameters are used.

    Returns
    -------
    matplotlib.Axes
        Scatter plot of means.
    """
    fig, ax = plt.subplots()
    for each in data[color].unique():
        subset = data[data[color] == each]
        mean_sub = subset.mean()
        sdev_sub = subset.std()
        ax.errorbar(mean_sub[x], mean_sub[y],
                    xerr=sdev_sub[x], yerr=sdev_sub[y],
                    label=each, **error_kwargs)
    ax.legend(loc='upper left', bbox_to_anchor=(1.0 , 1.05),
              frameon=False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    return ax

def jitter_values(values, val_range, bins=1):
    """Apply a uniform jitter to supplied values."""
    out = val_range / bins * 0.05 \
        * np.random.uniform(size=len(values)) \
        * np.random.choice([1, -1], size=len(values))
    return out

def error_plot(data, x, y, color, groupby=None, plot_points=False,
               jitter=False, dodge=False, error_kwargs={},
               scatter_kwargs={}, bottom_legend=False):
    """
    Plot an xy trendplot with error bars.

    Parameters
    ----------
    data : pd.DataFrame
        Dataframe containing data to plot.
    x : str
        Column in `data` to plot along the x-axis.
    y : str
        Column in `data` to plot along the y-axis.
    color : str
        Column in `data` to color trend lines by.
    groupby : str, optional
        Column in `data` to group observations by. By default None, and `x` will
        be used.
    plot_points : bool, optional
        Whether to plot points along with error bars. Default is False, and
        individual data points will not be plotted.
    jitter : bool, optional
        Whether to jitter values along the x-axis. Default is False.
    dodge : bool, optional
        Whether to dodge values by `x` or `groupby` if supplied. Default is
        False.
    error_kwargs : dict, optional
        Key-word arguments to pass to `plt.errorbar`.
    scatter_kwargs : dict, optional
        Key-word arguments to pass to `plt.scatter`. Only used if
        `plot_points==True`.

    Returns
    -------
    matplotlib.Axis
        xy trendplot with y error bars.
    """
    data = data.copy()
    xrange = data[x].max() - data[x].min()
    color_shift = {color: 0 for color in data[color].unique()}
    if groupby is None:
        means = data.groupby([color, x])[y].mean().reset_index()
        sdev = data.groupby([color, x])[y].std().reset_index()
        n_bins = len(np.unique(data[x]))
    else:
        means = data.groupby([color, groupby])[[x, y]].mean().reset_index()
        sdev = data.groupby([color, groupby])[[x, y]].std().reset_index()
        n_bins = len(np.unique(data[groupby]))
    if dodge:
        delta = xrange / (n_bins * 2) * 0.5
        shifts = np.linspace(-delta, delta,
                             num=len(data[color].unique()))
        color_shift = {color: shifts[i] for i, color in enumerate(data[color].unique())}
        
    print(sdev)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i, method in enumerate(sorted(set(means[color]))):
        subset = data[data[color] == method]
        mean_sub = means[means[color] == method].sort_values(x)
        sdev_sub = sdev[means[color] == method].loc[mean_sub.index, :]
        # print(sdev_sub)
        ax.errorbar(mean_sub[x] + color_shift[method], mean_sub[y], yerr=sdev_sub[y],
                    color=colors[i], label=method, **error_kwargs)
        if plot_points:
            jit = 0 
            if jitter:
                jit = jitter_values(subset[x], xrange, n_bins)
            ax.scatter(x=subset[x] + color_shift[method] + jit,
                       y=subset[y], color=colors[i], label=method, **scatter_kwargs)
    if bottom_legend:
        plot_bottom_legend(ax)
    else:
        plot_right_legend(ax)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(-0.05, 1)
    return ax


def lisi_by_score(means, score='ARI', groupby=None, error_kwargs={},
                  scatter_kwargs={}, bottom_legend=True):
    # fig, ax = plt.subplots(figsize=(13, 8))
    fig, ax = plt.subplots()
    if groupby is not None:
        for each in means[groupby].unique():
            method = means[means[groupby] == each]
            mean_sub = method.mean()
            sdev_sub = method.std()
            ax.errorbar(mean_sub[score], mean_sub['LISI'],
                        xerr=sdev_sub[score], yerr=sdev_sub['LISI'],
                        label=each, **error_kwargs)
    else:
        for each in sorted(means.index.values):
            ax.scatter(means.at[each, score], means.at[each, 'LISI'],
                       label=each, **scatter_kwargs)
    if bottom_legend:
        plot_bottom_legend(ax)
    else:
        plot_right_legend(ax)

    ax.set_ylabel('LISI')
    ax.set_xlabel(f'{score}')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.tight_layout()
    return ax