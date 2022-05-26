import os
import re
import json
import sys

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from cycler import cycler
from scanpy import api as sc
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns

from ncfs import distances
from downstream.src.visualization import visualize

sys.path.append(os.path.dirname(__file__))
import plotting, utils
sns.color_palette()

def color_point(point, scale):
    point = np.array(point)
    assert point.size == 3
    scaled = (point) / scale
    cmaps = [plt.get_cmap('Blues'), plt.get_cmap("Reds"), plt.get_cmap('Greens')]
    rgba = np.array([cmap(p) for cmap, p in zip(cmaps, scaled)])
    out = (np.median(rgba, axis=0)) * np.hstack((scaled, [1])) 
    return out

def get_benchmark_palette(obs):
    lines = obs.groupby('mixture')[['H1975', 'H2228', 'HCC827']].mean()
    palette = []
    for each in lines.index.values:
        palette.append(color_point(lines.loc[each, :].values, 9))
    return palette

def plot_umaps(adata, plotdir, label, treatment, n_neighbor, metric,
               dataset):
    # create AnnData object
    adata = sc.AnnData(X=X, obs=obs)
    # check to see if UMAP projection is already saved in obs data
    # calculate otherwise
    umap_cols = set(['UMAP1', 'UMAP2'])
    if umap_cols.intersection(obs.columns) != umap_cols:
        metric_kwds = {}
        if metric != 'euclidean':
            metric_kwds = {'w': np.ones(adata.shape[1])}

        sc.pp.pca(adata)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors,
                        metric=metric, metric_kwds=metric_kwds)
        sc.tl.umap(adata, min_dist=0.0)
        adata.obs['PC1'] = adata.obsm['X_pca'][:, 0]
        adata.obs['PC2'] = adata.obsm['X_pca'][:, 1]
        adata.obs['UMAP1'] = adata.obsm['X_umap'][:, 0]
        adata.obs['UMAP2'] = adata.obsm['X_umap'][:, 1]
    sc.settings.figdir = plotdir
    adata.obs[label] = adata.obs[label].astype('str').astype('category')
    print(f'labels: {adata.obs[label].unique()}')
    # plot umap colored by known cell type
    if dataset == 'benchmark':
        sc.pl.umap(adata, color=label, palette=get_benchmark_palette(adata.obs))
        plotting.close_plot()
    else:
        sc.pl.umap(adata, color=label, save='_known_cells.png')
        plotting.close_plot()
    # plot umap colored by cluster
    sc.pl.umap(adata, color='Cluster', save='_clusters.png')
    plotting.close_plot()
    # plot umap by 'treatment'
    # visualize.plot_umap(adata, color_col=treatment)
    # if dataset != 'benchmark':
    sc.pl.umap(adata, color=treatment, save='_treatment.png')
    plotting.close_plot()
    return adata.obs

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        dataset = snakemake.wildcards['dataset']
        identity = snakemake.params['identity'][dataset]
        treatment = snakemake.params['treatment'][dataset]
        control_id = snakemake.params['controls'][dataset]
        plotdir = snakemake.params['plotdir']
        sc.settings.figdir = plotdir
        with open(snakemake.input['fit'], 'r') as f:
            fit = json.load(f)
        method_dist = {'icat': ('phi_s', distances.phi_s),
                       'seurat311': ('euclidean', 'euclidean'),
                       'scanorama': ('euclidean', 'euclidean'),
                       'icat_scan': ('phi_s', distances.phi_s),
                       'seurat_icat': ('phi_s', distances.phi_s),
                       'no-int': ('phi_s', distances.phi_s)}
        
        performances = {}
        prefix = os.path.commonpath(snakemake.input['obs']).replace('/', r'\/')
        method_regex = re.compile(r'(?<={}\/)(.*?)(?=\/)'.format(prefix))
        obs_files = sorted(snakemake.input['obs'])
        X_files = sorted(snakemake.input['X'])
        for i, each in enumerate(obs_files):
            method = method_regex.search(each).group()
            X = np.loadtxt(X_files[i], delimiter=',')
            obs = pd.read_csv(each, index_col=0)
            adata = sc.AnnData(X=X, obs=obs)
            # standardize cluster column to Cluster
            # rename cluster column to be consistent between methods
            dist_name = method_dist[method][0]
            dist_func = method_dist[method][1]
            n_neighbors = fit[dist_name]['n_neighbors']
            col_name = utils.label_dictionary[method]
            obs.rename(columns={col_name: 'Cluster'}, inplace=True)
            obs['Cluster'] = obs['Cluster'].astype(str)
            # plot umaps
            umap_dir = os.path.join(plotdir, method)
            if not os.path.exists(umap_dir):
                os.makedirs(umap_dir)
            obs = plot_umaps(adata, umap_dir, identity, treatment, n_neighbors,
                             dist_func, dataset)
            # obs.to_csv(f"data/results/{dataset}/{method}_cells.csv")
            obs.to_csv(os.path.join('data', 'results', f'{dataset}', 'obs',
                                    f"{method}_cells.csv"))
                
            
            performances[method] = utils.performance(adata, identity, 'Cluster')
            # if dataset == 'Kang':
            #     for each in ['davies', 'silhouette']:
            #         performances[each] = labelless[each]
        plt.rcParams['axes.prop_cycle'] = cycler(color=plotting.palette)
        scores = pd.DataFrame(performances).T
        # scale calinksi between 0 and 1 for same axis as other metrics
        scores['calinski'] = MinMaxScaler().fit_transform(
                                scores['calinski'].values.reshape(-1, 1))
        # smaller values are better for davies, so take the negative before
        # min-max scaling 
        scores['DB'] = MinMaxScaler().fit_transform(
                                -1 * scores['DB'].values.reshape(-1, 1))
        scores.rename(columns=plotting.metric_dictionary,
                      index=plotting.method_dictionary,
                      inplace=True)
        lisi = pd.read_csv(snakemake.input['lisi'])\
                 .set_index('method', drop=True)\
                 .rename(index=plotting.method_dictionary)
        print(f"lisi: {lisi.index}")
        scores['LISI'] = lisi['LISI']
        scores = plotting.method_order(scores)
        metrics = ['ARI', 'DB', 'LISI']
        # plot lollipop plot
        plotting.lollipop_plot(scores, metrics)
        plt.savefig(snakemake.output['metrics'])
        plotting.close_plot()
        

        plotting.lisi_by_score(scores, 'ARI',
                               scatter_kwargs={'s': 200},
                               bottom_legend=False)
        plt.savefig(snakemake.output['lisi'])
        scores.to_csv(snakemake.output['csv'])
        

        
        
