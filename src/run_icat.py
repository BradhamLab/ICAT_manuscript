import pandas as pd
import numpy as np
import scanpy as sc
from icat import models
import json

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        var = None
        if snakemake.input['var'] is not None:
            var=pd.read_csv(snakemake.input['var'], index_col=0)
        adata = sc.AnnData(X=np.loadtxt(snakemake.input['X'],
                                        delimiter=','),
                           obs=pd.read_csv(snakemake.input['obs'],
                                           index_col=0),
                           var=var)
        treatment = snakemake.params['treatment']
        ctrl_val = snakemake.params['controls']
        colname = snakemake.params['cluster']
        with open(snakemake.input['icat'], 'r') as f:
            icat_kws = json.load(f)
        if not snakemake.params['simulated']:
            treatment = snakemake.params['treatment'][snakemake.wildcards['dataset']]
            ctrl_val = snakemake.params['controls'][snakemake.wildcards['dataset']]
            if snakemake.wildcards['dataset'] in ['Kang', 'Kagohara']:
                icat_kws['train_size'] = 1000 / (adata.obs[treatment] == ctrl_val).sum()
        with open(snakemake.input['json'], 'r') as f:
            fit_data = json.load(f)
        icat_kws['neighbor_kws']['n_neighbors'] = fit_data['phi_s']['n_neighbors']
        icat_kws['cluster_kws']['resolution'] = fit_data['phi_s']['resolution']
        if 'seurat' in snakemake.params['outdir'] or 'scanorama' in snakemake.params['outdir']:
            icat_kws['neighbor_kws']['n_neighbors'] = fit_data['euclidean']['n_neighbors']
            icat_kws['neighbor_kws']['metric'] = 'euclidean'
            icat_kws['cluster_kws']['resolution'] = fit_data['euclidean']['resolution']

        icat_model = models.icat(ctrl_val, **icat_kws)
        out = icat_model.cluster(adata, adata.obs[treatment], verbose=True)
        # sc.tl.uamp(out)
        out.X = out.obsm['X_icat'] # write to data for ease of snakemake output
        out.obs.rename(columns={'sslouvain': colname},
                       inplace=True)
        out.write_csvs(dirname=snakemake.params['outdir'],
                       skip_data=False)
        
        
