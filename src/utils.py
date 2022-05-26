import collections
import inspect
import os
import re
import warnings
import logging
import platform
import psutil

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from cycler import cycler
from matplotlib import pyplot as plt
from sklearn import metrics
import igraph as ig


label_dictionary = {
    'icat': 'sslouvain',
    'scanorama': 'scanorama.louvain',
    'icat_scan': 'scanorama.sslouvain',
    'seurat311': 'seurat_clusters',
    'seurat_icat': 'seurat.sslouvain',
    'no-int': 'no-int',
}

def log_system_info():
    uname = platform.uname()
    msg = "\n\nMachine Information" \
          f"System: {uname.system}\nNode Name: {uname.node}\n" \
          f"Release: {uname.release} Version: {uname.version}\n" \
          f"Machine: {uname.machine} Processor: {uname.processor}\n\n" \
          "CPU Information\n" \
          f"Physical Cores: {psutil.cpu_count(logical=False)}\n" \
          f"Total Cores: {psutil.cpu_count(logical=True)}\n\n\n\n"
    logging.info(msg)

def check_kws(reference_dict, new_dict, name):
    if not isinstance(new_dict, dict):
        raise ValueError("Expected dictionary of keyword arguments for "
                         "`{}`. Received {}.".format(name, type(new_dict)))
    for key, item in new_dict.items():
        if key not in reference_dict.keys():
            raise ValueError("Unsupported keyword argument `{}` for "
                             "{} keywords.".format(key, name))
        new_dict[key] = item
    return new_dict


def check_matching_genes(ref, new):
    return set(ref.var.index.values).difference(new.var.index.values) == 0


def get_default_kwargs(func, ignore_params=[]):
    params = inspect.signature(func).parameters
    kwargs = {x:params[x].default for x in params if x not in ignore_params}
    return kwargs


def check_np_castable(obj, name):
    """Check whether an object is castable to a numpy.ndarray."""
    if not isinstance(obj, np.ndarray):
        try:
            obj = np.array(obj)
        except:
            raise ValueError("Expected numpy.ndarray castable object for "
                            "`{}`. Got {}.".format(obj, type(obj)))
    return obj


def __evaluate_key(key, sep):
    if not isinstance(key, str):
        raise ValueError("Keys must be strings for easy concatentation.")
    if sep in key:
        raise ValueError("Cannot have `{}` in dictionary keys.".format(sep))


def flatten_dict(d, parent_key='', sep='.'):
    """
    Flatten a dictionary containing nested dictionaries.
    
    Parameters
    ----------
    d : dict
        Dictionary to flatten
    parent_key : str, optional
        Key in parent dictionary pointing to `d`. The default is '', which
        assumes `d` is the highest level nested dictionary.
    sep : str, optional
        String value to separate child and parent keys. The default is '.',
        which will place a '.' between each key. All parent and child keys
        will be assessed to ensure they do not contain a `sep` character;
        therefore, `sep` should be set to a delimiter not present in current
        keys.
    
    Returns
    -------
    dict
        Flattened dictionary with parent and child keys separted by `sep`.

    References
    ----------

    Taken shamelessly from here:
        https://stackoverflow.com/questions/6027558/flatten-nested-python-dictionaries-compressing-keys
    """

    items = []
    for k, v in d.items():
        __evaluate_key(k, sep)
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def score_population(obs, pop, label, cluster):
    """Calculate f1 index for a given population."""
    table = pd.crosstab(obs[cluster], obs[label])
    subset = obs[obs[cluster] == table[pop].idxmax()]
    all_pop = [True] * subset.shape[0]
    are_pop = subset[label] == pop
    return metrics.f1_score(are_pop, all_pop)


def performance(data, true_col, pred_col):
    """Measure the performance of a clustering partition."""
    is_adata = False
    if isinstance(data, sc.AnnData):
        obs = data.obs
        is_adata = True
    elif isinstance(data, pd.DataFrame):
        obs = data
    else:
        raise TypeError("Unsupported type: {}".format(type(data)))
    known = obs[true_col].values.astype(str)
    pred = obs[pred_col].values.astype(str)
    ar = metrics.adjusted_rand_score(known, pred)
    measures = {'ARI': float(ar),
                'ncluster': obs[pred_col].nunique()}

    # score each population by finding the cluster most labels are assigned to
    # and calculating a binary jaccard score
    for each in set(known):
        f1 = score_population(obs, each, true_col, pred_col)
        measures[f"{each}.f1"]= f1
    if is_adata:
        try:
            data.obsm['X_pca']
        except (ValueError, KeyError):
            sc.pp.pca(data)
        if len(np.unique(pred)) >= 2:
            calinski = metrics.calinski_harabasz_score(data.obsm['X_pca'], pred)
            davies = metrics.davies_bouldin_score(data.obsm['X_pca'], pred)
        else:
            calinski = np.nan
            davies = np.nan
        measures['calinski'] = calinski
        measures['DB'] = davies
    return measures

def igraph_from_adjacency(adjacency, directed=None):
    """
    Get igraph graph from adjacency matrix.
    
    Parameters
    ----------

    adjacency : numpy.ndarray, scipy.sparse.csr
        Adjacency matrix where non-zero entries represent connections between 
        samples.
    directed: boolean, optional

    Returns
    -------
    igraph.graph
        Nearest neighbor graph generated from adjacency matrix.

    References
    ----------
    Taken from: https://github.com/theislab/scanpy/blob/28498953092dc7cbecd0bd67380b1b060367d639/scanpy/_utils.py#L170
    """
    import igraph as ig
    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es['weight'] = weights
    except:
        pass
    if g.vcount() != adjacency.shape[0]:
        warnings.warn(
            f'The constructed graph has only {g.vcount()} nodes. '
            'Your adjacency matrix contained redundant nodes.'
        )
    return g

def is_none(x):
    """Check whether a value is a null value."""
    if isinstance(x, float):
        return np.isnan(x)
    return x is None

def format_labels(clusters):
    """
    Format cluster labels for sslouvain.

    Parameters
    ----------
    clusters : iterable
        List of cluster labels where None or np.nan represent previously
        unlabeled samples.
    
    Returns
    -------
    (list, list):
        labels : list
            List of labels where previously unlabeled samples are given 
            their own labels
        mutables : list
            List of samples indicating which samples were previously labelled.
    """
    if isinstance(clusters, pd.Series):
        clusters = clusters.values
    elif isinstance(clusters, list):
        clusters = np.array(clusters)
    elif not isinstance(clusters, np.ndarray):
        warnings.warn(f"Unsupport type {type(clusters)} for `clusters`")
    clusters = clusters.astype(float)
    mutables = [True] * len(clusters)
    labels = [None] * len(clusters)
    start_label = int(np.nanmax(clusters) + 1)
    for i, x in enumerate(clusters):
        if is_none(x):
            labels[i] = start_label
            start_label += 1
        else:
            labels[i] = int(x)
            mutables[i] = False
    return (labels, mutables)

def plot_umap(adata, color, shape, ax=None):
    if ax is None:
        __, ax = plt.subplots(figsize=(10, 8))
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    pallete = {}
    adata.obs[color] = adata.obs[color].astype(str)
    for i, each in enumerate(adata.obs[color].unique()):
        pallete[each] = colors[i]
    figure = sns.scatterplot(x=adata.obsm['X_umap'][:, 0],
                             y=adata.obsm['X_umap'][:, 1],
                             hue=adata.obs[color],
                             palette=pallete,
                             style=adata.obs[shape].astype(str),
                             ax=ax,
                             s=150)
    legend = ax.get_legend()
    for i, handle in enumerate(legend.legendHandles):
        if handle.get_label() == shape:
            legend.legendHandles[i].set_facecolor('white')
            legend.legendHandles[i].set_color('white')
            legend.legendHandles[i].set_edgecolor('white')

    return figure

def run_info(wildcards):
    pert_re = re.compile("Perturbation[0-9].*")
    number_regex = re.compile('[0-9]')
    out = {}
    experiment = wildcards['experiment']
    # out['Perturbation'] = pert_re.search(experiment).group(0)
    pert = pert_re.search(experiment)
    if pert is not None:
        experiment = experiment[:pert.span()[0]]
    out['Experiment'] = experiment
    out['Sim'] = int(number_regex.search(wildcards['sim']).group(0))
    out['Rep'] = int(number_regex.search(wildcards['rep']).group(0))
    return out
