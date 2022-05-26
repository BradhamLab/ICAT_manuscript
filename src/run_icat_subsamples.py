import pandas as pd
import numpy as np
import scanpy as sc
from icat import models
import json
import os
import re

sys.path.append(os.path.basename(__file__))
import utils


def measure_performance(data, labels, clusters, cells):
    performance = utils.performance(data, labels, clusters)
    performance['cells'] = cells
    return pd.DataFrame(performance, index=[cells])

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        dataset = snakemake.wildcards['large_ds']
        treatment = snakemake.params['treatment'][dataset]
        ctrl_val = snakemake.params['controls'][dataset]
        population = snakemake.params['label'][dataset]
        with open(snakemake.input['icat'], 'r') as f:
            icat_kws = json.load(f)
        adata = sc.AnnData(X=np.loadtxt(snakemake.input['X'],
                                        delimiter=','),
                           obs=pd.read_csv(snakemake.input['obs'],
                                            index_col=0))

        with open(snakemake.input['fit'], 'r') as f:
            fit_data = json.load(f)
        icat_kws['neighbor_kws']['n_neighbors'] = fit_data['phi_s']['n_neighbors']
        icat_kws['cluster_kws']['resolution'] = fit_data['phi_s']['resolution']
        ctrl_size = (adata.obs[treatment] == ctrl_val).sum()
        cells = int(snakemake.wildcards['cells'])
        # for i, cells in enumerate(snakemake.params['n_cells']):
        icat_kws['train_size'] = cells / ctrl_size
        # print(icat_kws['train_size'])
        icat_model = models.icat(ctrl_val, **icat_kws)
        clustered = icat_model.cluster(adata.copy(), adata.obs[treatment],
                                        verbose=True)
        performance = measure_performance(clustered,
                                          population,
                                          'sslouvain', cells)
        performance['dataset'] = dataset
        performance.to_csv(snakemake.output['csv'])
        clustered.write_csvs(snakemake.params['outdir'], skip_data=False)      
        
