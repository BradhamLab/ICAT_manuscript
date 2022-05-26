"""
Functions to create single-cell datasets.

Combines data from Kang et al. 2018.

author: Dakota Y. Hawkins
contact: dyh0110@bu.edu
"""
import pandas as pd
from scanpy import api as sc
from scipy.io import mmread
from glob import glob
import os
import pickle as pkl
from scipy import sparse

from downstream.src.analysis import utils as dutils

def create_count_matrix(matrix_file, barcode_file, genes):
    """
    Create a count matrix from necessary files.
    
    Parameters
    ----------
    matrix_file : str
        A sparse matrix file with 3 columns separated by spaces. In order, the
        expected columns are gene number, cell number, and umi count.
    barcode_file : str
        A single-column file containing all cell barcodes. Expected to be in
        order, such that "cell 1" listed in the matrix file will be the first
        entry.
    genes : pd.DataFrame
        A pandas dataframe containing gene annotations. Order is expected, such
        that 'gene 1' in the matrix file will be the first indexed gene.
    
    Returns
    -------
    pd.DataFrame
        An (n x p) data frame where n is the number of cells and p is the
        number of genes with at least one non-zero read.
    """
    # read in cell barcodes
    barcodes = pd.read_csv(barcode_file, names=['cell.barcode'])
    # match index numbers to names, indices started with 1
    idx_to_gene = {i + 1:x for i, x in enumerate(genes.index)}
    idx_to_bc = {i + 1:barcodes.loc[i, 'cell.barcode'] for i in barcodes.index}
    data = pd.read_csv(matrix_file, delimiter=' ', skiprows=3,
                       names=['gene', 'cell', 'count'])
    matrix = pd.pivot_table(data=data, index='cell', columns='gene',
                            values='count', fill_value=0)
    print(matrix.shape)
    matrix = matrix.rename(index=idx_to_bc, columns=idx_to_gene)
    return matrix

# def assign_cell_type(x):
#     perturbed_states = ['FCGR3A+ Monocytes', 'CD4 T cells']
#     if x['stim'] == 'stim' and x['cell_type'] in perturbed_states:
#         return 'Stimulated ' + x['cell_type']
#     return x['cell_type']

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        matrices = []
        barcodes = []
        genes = pd.read_csv(snakemake.input['genes'],
                            names=['ensmbl.id', 'name'], delimiter='\t',
                            index_col=0, header=None)
        cells = pd.read_csv(snakemake.input['cells'],
                            index_col=0, delimiter='\t')
        for mtx, bc in zip(snakemake.input['mtx'], snakemake.input['barcodes']):
            matrices.append(mmread(mtx).T)
            barcodes.append(pd.read_csv(bc, header=None, index_col=0))
        shared = list(set(barcodes[0].index).intersection(set(barcodes[1].index)))
        barcodes[0].rename(index={x:x+'1' for x in shared}, inplace=True)
        barcodes = pd.concat(barcodes)
        barcodes = barcodes.join(cells)
        X = sparse.vstack(matrices)
        adata = sc.AnnData(X=X.toarray(),
                           obs=barcodes,
                           var=genes)
        adata = adata[adata.obs['multiplets'] == 'singlet', :].copy()
        adata.obs.index.name = 'cell-barcode'
        adata.obs.rename(columns={'cell': 'cell_type'}, inplace=True)
        adata = adata[adata.obs['cell_type'].notnull(), :].copy()
        with open(snakemake.input['filter_cells'], 'r') as f:
            remove_cells = [line.strip() for line in f]
        adata = adata[adata.obs.index[~adata.obs.index.isin(remove_cells)], :]
        # adata.obs['perturbed_cell_type'] = adata.obs.apply(lambda x: assign_cell_type(x),
        #                                                    axis=1)
        sc.pp.filter_genes(adata, min_cells=50)
        sc.pp.filter_cells(adata, min_genes=50)
        # normalize data
        sc.pp.normalize_total(adata)
        # log transform counts to detected highly variable genes
        sc.pp.log1p(adata)
        # select highly variable
        sc.pp.highly_variable_genes(adata, flavor='seurat', batch_key='stim')
        adata.var.highly_variable = adata.var.highly_variable_nbatches > 0
        adata = adata[:, adata.var.highly_variable]
        adata.write_csvs(dirname=snakemake.params['outdir'], skip_data=False)

    
