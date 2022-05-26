import re
import pandas as pd
import numpy as np

import plotting

def experiment_group(x):
    return re.sub('[0-9]', '', x)
    
# method_dictionary = {
#     'icat': 'ICAT',
#     'seurat233': 'Seurat 2.3',
#     'seurat.aligned': 'Seurat 2.3 - Aligned',
#     'seurat311': 'Seurat 3.1',
#     'scanorama': 'scanorama',
#     'icat_scan': 'scanorama + ICAT',
#     'seurat_icat': 'Seurat 3.1 + ICAT',
#     'no-int': 'No Integration',
#     'ncfs-louvain': 'NCFS + Louvain'
# }

# # metric_dictionary = {
# #     'adjusted.mutual.info': 'AMI',
# #     'adjusted.rand': 'ARI',
# #     'completeness': 'Completeness',
# #     'fowlkes.mallows': 'Fowlkes-Mallows',
# #     'homogeneity': 'Homogeneity',
# # }

def format_simulated_output(output, method_dict=plotting.method_dictionary,
                            column_dict=plotting.metric_dictionary):
    output = output.copy()
    output['method'] = output['method'].replace(method_dict)
    output['Percent Perturbed'] *= 100
    output.rename(columns = {'Percent Perturbed': 'Perturbed Genes (%)'},
                  inplace=True)
    output['by_sim'] = output.apply(lambda x: x['Experiment'] + str(x['Sim']),
                                    axis=1)
    by_sim = output.groupby(['by_sim', 'method'])
    str_data = by_sim.agg(lambda col: ';'.join(np.unique(col)))
    perf = by_sim.mean()
    
    merged = perf.merge(str_data, on=['by_sim', 'method'])
    merged.rename(columns=column_dict, inplace=True)
    merged['ExpType'] = merged.apply(lambda x: experiment_group(x['Experiment']),
                                     axis=1)
#     merged['Perturbed'] = merged.apply(lambda x: "Perturbation" in x['Experiment'],
#                                        axis=1)
#     merged['Run'] = merged.apply(lambda x: x['Experiment'] + x['Perturbation'],
#                                  axis=1)
    drop_cols = ['by_sim', 'Rep', 'Run']
    out = merged.reset_index().drop(columns=drop_cols)
    out.rename(columns={'method': "Method"}, inplace=True)
    if any(out['Method'].isin(['NCFS-SSLouvain'])):
        out = out[out['Method'] != 'NCFS-SSLouvain']
    return out

if __name__ == '__main__':
    results = pd.read_csv(snakemake.input['csv'], index_col=0)
    out = format_simulated_output(results)
    out.to_csv(snakemake.output['csv'])