"""Run louvain clustering."""
import json

import scanpy as sc
import pandas as pd
import numpy as np

from ncfs import distances


if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        X = np.loadtxt(snakemake.input['X'], delimiter=',')
        obs = pd.read_csv(snakemake.input['obs'], index_col=0)
        var = pd.read_csv(snakemake.input['var'], index_col=0)
        colname = snakemake.params['cluster']
        adata = sc.AnnData(X=X, obs=obs, var=var)

        with open(snakemake.input['json'], 'r') as f:
            fit_data = json.load(f)
        
        # cluster cells 
        sc.pp.pca(adata)
        sc.pp.neighbors(adata, n_neighbors=fit_data['phi_s']['n_neighbors'],
                        metric=distances.phi_s,
                        metric_kwds={'w': np.ones(adata.shape[1])})
        sc.tl.umap(adata)
        sc.tl.louvain(adata, key_added=colname)
        adata.write_csvs(dirname=snakemake.params['outdir'],
                         skip_data=False)
