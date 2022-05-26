import re
import sys
import os

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import MinMaxScaler

print(os.path.basename(__file__))
sys.path.append(os.path.basename(__file__))
import utils

def parse_input_path(filepath):
    directory = os.path.split(filepath)[0]
    path_split = directory.split(os.path.sep)
    simrep = path_split[5]
    sim_re = re.compile('Sim*.[0-9]')
    return({'method': path_split[3],
            'experiment': path_split[4],
            'sim': sim_re.search(simrep).group(0),
            'rep': simrep[sim_re.search(simrep).span()[1]:]})

def parse_performance_path(filepath):
    filename = os.path.basename(filepath)
    sim_re = re.compile('Sim*.[0-9]')
    rep_re = re.compile('Rep*.[0-9]')
    path_split = filepath.split(os.path.sep)
    return({'method': path_split[3],
            'experiment': path_split[4],
            'sim': sim_re.search(filename).group(0),
            'rep': rep_re.search(filename).group(0)})

def gini(pops):
    value = 0
    for x, y in itertools.product(pops, pops):
        value += abs(x - y)
    return value / (2 * len(pops) ** 2 * np.mean(pops))

def measure_performance(data, labels, wildcards):
    method = wildcards['method']
    performance = utils.performance(data, labels,
                                    utils.label_dictionary[method])
    performance['method'] = method
    run_info = utils.run_info(wildcards)
    for key, value in run_info.items():
        performance[key] = value
    out = pd.DataFrame(performance, index=["{}{}{}".format(wildcards['experiment'],
                                                           wildcards['sim'],
                                                           wildcards['rep'])])
    out['Run'] = out.index.values
    return out

# gine vs variance or something for genes to identify spatially isoloationed genes

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        labels = snakemake.params['identity']
        dfs = []
        # measure performance for each method + run
        for csv in snakemake.input['obss']:
            obs = pd.read_csv(csv, index_col=0)
            X = np.loadtxt(os.path.join(os.path.split(csv)[0], 'X.csv'),
                           delimiter=',')
            lisi = pd.read_csv(csv.replace('obs.csv', 'lisi.csv'), index_col=0)
            adata = sc.AnnData(X=X, obs=obs)
            wildcards = parse_input_path(csv)
            key = ''.join([wildcards[x]\
                           for x in ['method', 'experiment', 'sim', 'rep']])
            perf = measure_performance(adata, labels, wildcards)
            perf['LISI'] = lisi['LISI'].values[0]
            perf['Activated.LISI'] = lisi['LISI'].values[0]
            dfs.append(perf)
        meta = pd.read_csv(snakemake.input['meta'])
        combined = pd.concat(dfs, axis=0, sort=False).reset_index(drop=True)
        merged = combined.merge(meta, on='Run')
        merged.to_csv(snakemake.output['csv'])