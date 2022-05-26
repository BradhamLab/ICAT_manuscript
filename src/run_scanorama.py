import pandas as pd
import scanpy as sc
import json
import numpy as np
import scanorama


def run_scanorama(adatas, k, r):
    integrated_X = scanorama.integrate_scanpy(adatas)
    integrated = []
    for i, each in enumerate(integrated_X):
        adata = sc.AnnData(X=each, obs=adatas[i].obs)
        integrated.append(adata)
    data = integrated[0].concatenate(integrated[1:])
    sc.pp.pca(data)
    sc.pp.neighbors(data, n_neighbors=k)
    sc.tl.umap(data)
    sc.tl.louvain(data, resolution=r, key_added='scanorama.louvain')
    return data

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        adata = sc.AnnData(X=np.loadtxt(snakemake.input['X'], delimiter=','),
                           obs=pd.read_csv(snakemake.input['obs'], index_col=0),
                           var=pd.read_csv(snakemake.input['var'], index_col=0))
        treatment = snakemake.params['treatment']
        ctrl_val = snakemake.params['controls']
        if not snakemake.params['simulated']:
            treatment = snakemake.params['treatment'][snakemake.wildcards['dataset']]
            ctrl_val = snakemake.params['controls'][snakemake.wildcards['dataset']]
        adatas = []
        for each in adata.obs[treatment].unique():
            subset = adata[adata.obs[treatment] == each, :].copy()
            print(subset.shape)
            if ctrl_val == each:
                adatas = [subset] + adatas
            else:
                adatas.append(subset)
        with open(snakemake.input['json'], 'r') as f:
            fit_data = json.load(f)
        k = fit_data['euclidean']['n_neighbors']
        r = fit_data['euclidean']['resolution']
        out = run_scanorama(adatas, k, r)
        out.var.index = ['scan-{}'.format(i + 1) for i in range(out.shape[1])]
        out.write_csvs(dirname=snakemake.params['outdir'], skip_data=False)