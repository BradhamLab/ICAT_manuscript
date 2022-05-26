"""Find best K for phi_s and euclidean distances."""

import itertools
import sys
import os

import scanpy as sc
import matplotlib.pyplot as plt
from cycler import cycler

import json
import pandas as pd
import numpy as np
import colorcet as cc

sys.path.append(os.path.split(__file__)[0])

import utils
from ncfs import distances

try:
    loc = os.path.dirname(os.path.abspath(__file__))
    plt.style.use(os.path.join(loc, 'configs/figures.mplstyle'))
    plt.rc('axes', prop_cycle=cycler('color', cc.glasbey_light))
except:
    pass

def main(adata, label_col):
    sc.pp.pca(adata)
    performance = {'phi_s': {'n_neighbors': None,
                             'resolution': None},
                   'euclidean': {'n_neighbors': None,
                                 'resolution': None}}
    n_max = int(0.5 * adata.shape[0])
    if n_max > 50:
        n_max = 50
    n_min = 5
    scores = {'phi_s': -np.inf,
              'euclidean': -np.inf}
    resolutions = np.arange(0.2, 1.3, 0.05)
    # #
    # if label_col == 'mixture':
    #     resolutions = np.arange(0.5, 1.55, 0.05)
    # resolutions = [1]
    metrics = [distances.phi_s, 'euclidean']
    neighbors = range(n_min, n_max + 5, 5)
    weights = np.ones(adata.X.shape[1])
    sc.pp.pca(adata)
    for metric, n, r in itertools.product(metrics, neighbors, resolutions):
        metric_kwds = {}
        if metric == distances.phi_s:
            metric_kwds = {'w': weights}
        sc.pp.neighbors(adata, n_neighbors=n, metric=metric,
                        metric_kwds=metric_kwds)
        sc.tl.louvain(adata, resolution=r, key_added='louvain')
        # add sil score using eucldiean distance on pcas
        measures = utils.performance(adata, label_col, 'louvain')
        if metric == distances.phi_s:
            metric = 'phi_s'
        # if scores[metric] < measures['adjusted.rand']:
        #     performance[metric].update(measures)
        #     performance[metric]['n_neighbors'] = n
        #     performance[metric]['resolution'] = r
        #     scores[metric] = measures['adjusted.rand']
        if np.isnan(measures['calinski']):
            measures['calinski'] = -1
        if scores[metric] < measures['calinski']:
            performance[metric].update(measures)
            performance[metric]['n_neighbors'] = n
            performance[metric]['resolution'] = r
            scores[metric] = measures['calinski']
    return performance


if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        X = np.loadtxt(snakemake.input['X'], delimiter=',')
        obs = pd.read_csv(snakemake.input['obs'], index_col=0)
        var = pd.read_csv(snakemake.input['var'], index_col=0)
        adata = sc.AnnData(X=X, obs=obs, var=var)
        
        treat_col = snakemake.params['treatment']
        ctrl_val = snakemake.params['control']
        label = snakemake.params['label']

        if not snakemake.params['simulated']:
            treat_col = snakemake.params['treatment'][snakemake.wildcards['dataset']]
            ctrl_val = snakemake.params['control'][snakemake.wildcards['dataset']]
            label = snakemake.params['label'][snakemake.wildcards['dataset']]
            
        controls = adata[adata.obs[treat_col] == ctrl_val, :]
        fit_data = main(controls, label)
        with open(snakemake.output['json'], 'w') as f:
            json.dump(fit_data, f, indent=4)
