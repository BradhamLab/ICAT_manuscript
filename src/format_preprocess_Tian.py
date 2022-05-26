import pandas as pd
import os
import numpy as np
from scanpy import api as sc
import re

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        adatas = []
        genes = set()
        # lines in CelSeq benchmark data
        lines = ['H1975', 'H2228', 'HCC827']
        # cell mixtures are mixtures of 9 cells
        line_to_mixture = {lines[0]:[9,0,0],
                           lines[1]:[0,9,0],
                           lines[2]:[0,0,9]}
        bench_regex = re.compile('^(.*?)\.')
        for x, meta in zip(snakemake.input['counts'], snakemake.input['meta']):
            count = pd.read_csv(x, index_col=0)
            obs = pd.read_csv(meta, index_col=0)
            # flag for dataset
            bench = os.path.basename(x[:bench_regex.search(x).end() - 1])
            obs['benchmark'] = bench
            obs['index'] = obs.apply(lambda x: '-'.join([x.benchmark, x.name]),
                                     axis=1)
            obs.set_index('index', inplace=True)

            if len(genes) == 0:
                genes = set(count.index.values)
            else:
                genes = genes.intersection(count.index.values)
            
            adatas.append(sc.AnnData(X=count.T.values, obs=obs,
                                     var=pd.DataFrame(count.index.values,
                                                      index=count.index.values,
                                                      columns=['gene'])))
        genes = list(genes)
        # subset genes down to shared genes
        for i in range(len(adatas)):
            adata = adatas[i][:, genes]
            # convert number of cells of each cell type in a mixture to a
            # mixture id. Remove mixtures with less than 9 cells to avoid
            # unclear identity.
            if all([x in adata.obs.columns for x in lines]):
                adata.obs['n_cells'] = adata.obs[lines].apply(lambda x: sum(x),
                                                              axis=1)
                adata.obs['mixture'] = adata.obs[lines].apply(lambda x:
                                                  ','.join([str(y) for y in x]),
                                                  axis=1)
                adata = adata[adata.obs['n_cells'] == 9, :].copy()
            # convert cell lines to mixture id of pure cell line
            else:
                adata.obs['mixture'] = adata.obs['cell_line'].apply(lambda x:
                                    ','.join([str(y) for y in line_to_mixture[x]]))
                adata.obs[lines[0]] = None
                adata.obs[lines[1]] = None
                adata.obs[lines[2]] = None 
                # if you know how to apply for multipel column assignment, lmk
                for x in adata.obs.index.values:
                    cell_type = adata.obs.loc[x, 'cell_line']
                    adata.obs.loc[x, lines] = line_to_mixture[cell_type] 
            adatas[i] = adata
        combined = adatas[0].concatenate(adatas[1:])
        # select mixtures by 3s, and 4,5 + 2,7 for more discrete populations
        if True:
            import itertools
            allowable_mixtures = ['9,0,0', '0,0,9', '0,9,0', '3,3,3']
            def keep_mixture(x):
                # return sum(x[lines] % 3) == 0 or x['mixture'] in allowable_mixtures
                return x['mixture'] in allowable_mixtures
            keep = combined.obs.apply(lambda x: keep_mixture(x), axis=1)
            combined = combined[keep, :]
            
        sc.pp.filter_genes(combined, min_cells=3)
        # sc.pp.filter_cells(combined, min_genes=5)
        # normalize data
        sc.pp.normalize_total(combined)
        # log transform counts to detected highly variable genes
        sc.pp.log1p(combined)
        # select highly variable
        sc.pp.highly_variable_genes(combined, flavor='seurat',
                                    batch_key='benchmark')
        combined.var['highly_variable'] = combined.var['highly_variable_nbatches'] > 0
        # subset to variable genes
        hvgs = combined[:, combined.var['highly_variable']].copy()
        hvgs.obs['mixture'] = hvgs.obs['mixture'].astype('category')
        # write data
        hvgs.write_csvs(dirname=snakemake.params['outdir'], skip_data=False)
        
                                        
        

