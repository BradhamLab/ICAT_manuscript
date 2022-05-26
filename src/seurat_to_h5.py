import pandas as pd
import numpy as np
import scanpy as sc

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
                adata = sc.AnnData(X=np.loadtxt(snakemake.input['X'], delimiter=','),
                                   obs=pd.read_csv(snakemake.input['obs'], index_col=0),
                                   var=pd.read_csv(snakemake.input['var'], index_col=0))
                sc.pp.pca(adata, n_comps=20)
                sc.pp.neighbors(adata, n_neighbors=15)
                sc.tl.umap(adata)
                adata.write(snakemake.output['adata'])