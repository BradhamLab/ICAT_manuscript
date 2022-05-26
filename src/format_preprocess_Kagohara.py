import pandas as pd
import numpy as np
import scanpy as sc
import os
import glob

from scipy import sparse



def adata_from_line(datadir, line):
    X_df = pd.read_csv(glob.glob(os.path.join(datadir, f"*{line}Matrix.csv"))[0], index_col=0).T
    obs = pd.read_csv(glob.glob(os.path.join(datadir, f"*phenoData{line}.csv"))[0], index_col=0)
    var = pd.read_csv(os.path.join(datadir, "GSE137524_featureData.csv"), index_col=0)
    adata = sc.AnnData(X=sparse.csr_matrix(X_df.values),
                       obs=obs, var=var)
    return adata

def combine_matrices(datadir):
    lines = ['SCC1', 'SCC25', 'SCC6']
    adata = None
    var = None
    for line in lines:
        new = adata_from_line(datadir, line)
        if adata is None:
            adata = new.copy()
            var = adata.var
        else:
            adata = adata.concatenate(new)
        del new
    adata.var = var
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, batch_key='treatment')
    adata = adata[:, adata.var['highly_variable_nbatches'] > 0].copy()
    adata.highly_variable = True
    sc.pp.scale(adata)
    sc.pp.combat(adata, 'replicate', covariates=['treatment'])
    return adata

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        adata = combine_matrices(snakemake.params['datadir'])
        adata.write_csvs(snakemake.params['outdir'], skip_data=False)
    