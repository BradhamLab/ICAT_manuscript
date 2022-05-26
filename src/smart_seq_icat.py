import scanpy as sc
import pandas as pd
import numpy as np
from icat import models
from ncfs import distances
from scipy import stats
from statsmodels.stats.multitest import multipletests



bad_cells = ["ASW_C07_2016_07_29", "ASW_C01_2016_07_29",
             "ASW_C09_2016_07_29", "ASW_A02_2016_07_29",
             "ASW_G06_2016_07_29", "ASW_E05_2016_07_29",
             "ASW_A12_2016_07_29", "ASW_G03_2016_07_29",
             "ASW_E04_2016_07_29", "ASW_E11_2016_07_29",
             "ASW_C02_2016_07_29", "Chlorate_PMCs_1_G05_2018_07_01",
             "Chlorate_PMCs_1_E06_2018_07_01", "MK886_PMCs_2_A05_2018_07_01",
             "MK886_PMCs_2_F10_2018_07_01", "MK886_PMCs_3_E05_2018_07_01",
             "MK886_PMCs_3_D01_2018_07_01", "MK886_PMCs_3_B11_2018_07_01",
             "MK886_PMCs_3_G10_2018_07_01", "MK886_PMCs_3_G05_2018_07_01"]

region_markers = [
    "evm.model.scaffold118_len3932288_cov177.83",
    "evm.model.scaffold1591_len1589934_cov177.43",
    "evm.model.scaffold2456_len1195570_cov177.18",
    "evm.model.scaffold2778_len965975_cov166.38",
    "evm.model.scaffold2799_len3283164_cov183.71",
    "evm.model.scaffold316954_len180040_cov163.4",
    "evm.model.scaffold463_len517687_cov169.9",
    "evm.model.scaffold486_len2917596_cov191.59",
    "evm.model.scaffold5303_len678482_cov176.3",
    "evm.model.scaffold679_len1722255_cov171.40",
    "evm.model.scaffold928_len929390_cov196.19",
    "evm.model.scaffold964_len748964_cov165.5",
    "evm.model.scaffold1102_len820463_cov161.12",
    "evm.model.scaffold1168_len1380116_cov183.17",
    "evm.model.scaffold1168_len1380116_cov183.18",
    "evm.model.scaffold150_len1782623_cov165.27",
    "evm.model.scaffold2184_len1117722_cov172.9",
    "evm.model.scaffold2799_len3283164_cov183.31",
    "evm.model.scaffold3917_len267786_cov197.2",
    "evm.model.scaffold1591_len1589934_cov177.39",
    "evm.model.scaffold316991_len500887_cov178.12",
    "evm.model.scaffold4274_len1671731_cov191.2",
    "evm.model.scaffold2456_len1195570_cov177.19",
    "evm.model.scaffold316377_len465659_cov199.9",
    "evm.model.scaffold7953_len991057_cov189.23",
    "evm.model.scaffold3329_len1799558_cov178.23",
    "evm.model.scaffold1678_len1363396_cov182.24",
]


def compare_clusters(anno_df, comparisons=None, cluster_col='louvain',
                     compare_col='treatment', merge_cols=None):
    """
    Compare cell make-up between identified clusters.
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe containing cluster membership for each observation
            and a second categorical variable to compare frequency against.
    comparisons : list, optional
        Post-hoc pairwise comparisons to perform. Elements of the list should be
            length 2 containers containing categories to compare (e.g.
            ('Control', 'Expermintanl')). Frequencies between conditions
            will be compared in each cluster using a pairwise G-Test (the
            default is None, which results in no post-hoc analysis). 
    cluster_col : str, optional
        Name of the column in `anno_df.obs` containing observation cluster/group
            assignment (the default is 'louvain'.)
    compare_col : str, optional
        Name of the column in `anno_df.obs` partitioning observations into
            groups of interest (the default is 'treatment', which compares
            distribution of treatments within clusters.)
    merge_cols : dict, optional
        Name of columns to merge together. Keys should point to list of column
            names to merge. The merged column name will be set the key (e.g.
            {'Control': ['ASW', 'DMSO']} will merge 'ASW' and 'DMSO' counts to 
            a single 'Control' column). Default is None, and no merging is
            performed.
    
    Returns
    -------
    dictionary
        Dictionary containing comparison results between clusters including
        G-tests of independence and, if a `comparisons` argument was provided,
        pairwise Fisher exact test results.
        key-value pairs:
            'gtest': dictionary of Gtest results
                key-value pairs:
                    'observed': matrix of observed counts.
                    'expected': matrix of expected counts.
                    'pvals': pvalues for a GTest of independence following
                        index order of 'observed' and 'expected' tables.
                    'pvals.adj': adjusted p-values using the bonferonni
                        correction.
            'fisher': dictionary of pairwise comparison results. Results are
                contained in a dictionary keyed by comparison groups separated
                by a hyphen (e.g. if ['X', 'Y'] was provided as a comparison,
                results of the comparison would be keyed 'X-Y').
                
                key-value pairs:
                    'odds': numpy.array of odds ratios between comparisons
                        ordered by cluster.
                    'pvals': numpy.array of calculated p-values orderd by
                        cluster.
                    'pvals.adj': numpy.array of adjusted p-values by benjamin-hochberg fdr.
                    'cluster': id denoting which cluster comparisons. 
    """
    # check for cluster column
    if cluster_col not in anno_df.obs.columns:
        raise ValueError('No {} column in observation data.'.format(
                         cluster_col))

    # check for comparison column
    if compare_col not in anno_df.obs.columns:
        raise ValueError('No {} column in observations data.'.format(
                         compare_col))
    
    #TODO check for comparisons in comparison column

    count_table = pd.crosstab(anno_df.obs[cluster_col],
                              anno_df.obs[compare_col])
    if merge_cols is not None and isinstance(merge_cols, dict):
        for key, columns in merge_cols.items():
            try:
                count_table[key] = count_table[columns].sum(axis=1)
            except KeyError:
                raise('Unknown columns: {}'.format(columns))
            count_table.drop(columns, axis=1, inplace=True)
        count_table = count_table[sorted(count_table.columns.values)]
        
    
    # calculate probabilities for comparison values
    probabilities = np.sum(count_table, axis=0) / np.sum(count_table.values)

    # calculate total number of cells per cluster in comparison
    group_counts = np.sum(count_table, axis=1)
    
    # matrix multipy cluster counts and condition probabilities to get
    # expected counts
    expectation = group_counts.values.reshape((count_table.shape[0], 1))\
                  @ probabilities.values.reshape((1, count_table.shape[1]))
    expectation = pd.DataFrame(data=expectation, index=count_table.index,
                               columns=count_table.columns)

    # perform tests of independence between treatments and 
    results = stats.power_divergence(count_table, expectation,
                                     axis=1, lambda_='log-likelihood')
    gtest_results = {'cluster': count_table.index.values,
                     'pvals': results.pvalue,
                     'pvals.adj': multipletests(results.pvalue, method='fdr_bh')[1]}

    out = {'gtest': pd.DataFrame.from_dict(gtest_results)}
    if comparisons is not None: # perform odds-ratio/fisher exact tests
        fisher_dfs = []
        for each in comparisons:
            if len(each) != 2:
                msg = ('Comparisons must be pairwise. Received'
                       ' {} groups: {}'.format(len(each), each))
                raise ValueError(msg)
            
            pairwise = count_table.loc[:, list(each)].copy()
            pairwise.index = pairwise.index.astype(str)
            out_key = '-'.join(each)
            results  = {'odds': np.ones(count_table.shape[0]),
                        'pvals': np.ones(count_table.shape[0]),
                        'pvals.adj': np.ones(count_table.shape[0]),
                        'cluster': count_table.index.values,
                        'comparison': [out_key] * count_table.shape[0]}
            for i, cluster in enumerate(pairwise.index):

                # create a 2 x 2 contigency table between treatment comparisons
                # and cluster membership. Counts are pairwise between treatments
                # and cluster X membership vs. not X
                test_cluster = pairwise.loc[cluster, :]
                other_clusters = [x for x in pairwise.index if x != cluster]
                not_cluster = pairwise.loc[other_clusters, :].sum()
                contingency = pd.concat((test_cluster, not_cluster), axis=1)

                # perform fisher exact's test
                odds, pval = stats.fisher_exact(contingency.values)
                results['odds'][i] = odds
                results['pvals'][i] = pval
            results['pvals.adj'] = multipletests(results['pvals'], method='fdr_bh')[1]
            fisher_dfs.append(pd.DataFrame.from_dict(results))
        out['fisher'] = pd.concat(fisher_dfs)
    return out


def get_gene_identifier(scaffold, gene_df):
    """
    Get gene name associated with scaffold id.

    Parameters
    ----------
    scaffold : str
        Scaffold transcript id in `gene_df` as an index value.
    gene_df : pd.DataFrame
        A (gene x feamodelre) dataframe with 'UniProt.Name', 'UniProt.ID', 'SPU',
        and 'NCBI.ID' columns containing different gene id/name values.

    Remodelrns
    -------
    str
        Gene name associated with `scaffold` id. Name priority follows
        'UniProt.Name', 'Uniprot.ID', 'SPU', 'NCBI.ID', `scaffold` in order. 
    """

    if not pd.isnull(gene_df.at[scaffold, 'UniProt.Name']) and\
    gene_df.at[scaffold, 'UniProt.Name'] != 'nan':
        return gene_df.at[scaffold, 'UniProt.Name'].split('_')[0]
    if not pd.isnull(gene_df.loc[scaffold, 'UniProt.ID']) and\
    gene_df.at[scaffold, 'UniProt.ID'] != 'nan':
        return gene_df.at[scaffold, 'UniProt.ID']
    if not pd.isnull(gene_df.loc[scaffold, 'SPU']) and\
    gene_df.at[scaffold, 'SPU'] != 'nan':
        return gene_df.at[scaffold, 'SPU']
    if not pd.isnull(gene_df.loc[scaffold, 'NCBI.ID']) and\
    gene_df.at[scaffold, 'NCBI.ID'] != 'nan':
        return gene_df.at[scaffold, 'NCBI.ID']
    return scaffold


if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        var = pd.read_csv(snakemake.input['var'],
                          index_col=0)
        var['Name'] = var.apply(lambda x: get_gene_identifier(x.name, var), axis=1)

        obs = pd.read_csv(snakemake.input['obs'],
                          index_col=0)
        obs.set_index(obs.apply(lambda x: x.name.replace('-', '_'),
                                axis=1),
                      inplace=True)
        obs = obs.loc[~obs.index.isin(bad_cells), :].copy()
        obs.replace({'treatment': {'ASW': 'Control', 'DMSO': 'Control'}},
                    inplace=True)
        norm_counts = pd.read_csv(snakemake.input['counts'],
                                  index_col=0)
        norm_counts.set_index(norm_counts.apply(lambda x: x.name.replace('TU', 'model'),
                                                axis=1),
                              inplace=True)
        norm_counts.columns = [x.replace('-', '_') for x in norm_counts.columns]
        matched = list(set(var.index.values).intersection(norm_counts.index.values))
        X = norm_counts.loc[matched, obs.index].T.values
        adata = sc.AnnData(X=X, obs=obs, var=var.loc[matched, :].copy())
        # adata.obs = obs
        print(f"obs columns: {adata.obs.columns}")
        print(f"adata dimensions: {adata.shape}")
        ctrls = adata[adata.obs['treatment'] == 'Control', :].copy()
        # forcibly consider known regional markers
        sc.pp.highly_variable_genes(ctrls,
                                    n_top_genes=2000,
                                    flavor='seurat',)
        ctrls.var.loc[region_markers, 'highly_variable'] = True
        model = models.icat('Control',
                            reference='all',
                            ncfs_kws={'reg': 0.5, 'sigma': 3},
                            pca_kws={'n_comps': 20},
                            cluster_kws={'resolution': 1},
                            neighbor_kws={'n_neighbors': 15})
        out = model.cluster(adata[:, ctrls.var.highly_variable].copy(),
                            adata.obs.treatment)
        out.obs['sslouvain'].astype(int, inplace=True)
        
        sc.tl.umap(out)
        import matplotlib.pyplot as plt
        ax = sc.pl.umap(out, color='sslouvain', show=False)
        plt.savefig('/projectnb/bradham/analysis/evaluate-icat/smart_seq_umap.png')
        dif_abundance = compare_clusters(out,
                                         comparisons=[['Chlorate', 'Control'],['MK886', 'Control']],
                                         cluster_col='sslouvain',
                                         compare_col='treatment',
                                         merge_cols=None)
        dif_abundance['gtest'].to_csv(snakemake.output['gtest'])
        dif_abundance['fisher'].to_csv(snakemake.output['fisher'])
        out.write(snakemake.output['adata'])
        np.savetxt(snakemake.output['X'], out.X, delimiter=',')
        out.obs.to_csv(snakemake.output['obs'])
        out.var.to_csv(snakemake.output['var'])
