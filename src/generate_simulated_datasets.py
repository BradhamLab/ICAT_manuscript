import json
import itertools
import os

import numpy as np
import pandas as pd

from icat import simulate
from icat import utils

def parse_params(sim_params):
    if 'dispersion' in sim_params and sim_params['dispersion'] == 'random':
        sim_params['dispersion'] = simulate.dispersions(sim_params['genes'],
                                                        a=1, b=4)
    return sim_params


def main(config, sims=1, reps=1, outdir='.', expname=''):
    csv_dict = dict()
    c_params = parse_params(config['control_kwargs'])
    p_params = parse_params(config['perturb_kwargs'])

    # create experiment object for control-perturbation pairing
    experiment = simulate.Experiment(control_kwargs=c_params,
                                     perturb_kwargs=p_params)
    # run simulations over all sims + reps 
    simmed = experiment.run(sims, reps)
    flattened =  dict(utils.flatten_dict(c_params),  # move flatten dict out 
                      **utils.flatten_dict({'perturbation': p_params}))
    for k, v in flattened.items():
        if isinstance(v, list):
            v = ';'.join([str(x) for x in v])
        flattened[k] = v
    # would be quicker to simulate data in loop
    for i, j in itertools.product(range(sims), range(reps)):
        sim_rep_data = flattened.copy()
        data = simmed[i][j]
        assert np.all(data.X.astype(int) == data.X)
        markers = np.hstack(list(simulate.population_markers(data).values()))
        sim_rep_data['n_markers'] = len(markers)
        sim_rep_data['dropout'] = np.sum(data.X == 0) / data.X.size
        sim_rep_data['Experiment'] = expname
        sim_rep_data['Sim'] = i + 1
        sim_rep_data['Rep'] = j + 1 
        sim_rep = "Sim{}Rep{}".format(i + 1, j + 1)
        # filename = os.path.join(outdir, f"{exp_rep_sim}.h5ad")
        data.obs['Population'].replace({'1': 'C1',
                                        '2': 'C2',
                                        '3': 'C3',
                                        'Perturbed-1': 'C1+',
                                        'Perturbed-added-1': 'P4',
                                        'Perturbed-added-2': 'P5'},
                                        inplace=True)
        data.write_csvs(os.path.join(outdir, sim_rep), skip_data=False)
        # print(f"FILENAME {filename}")
        # data.write(filename=filename)
        csv_dict["{}{}".format(expname, sim_rep)] = sim_rep_data

    return pd.DataFrame(csv_dict).T
    

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        with open(snakemake.input['json']) as f:
            configs = json.load(f)
        print(snakemake.wildcards['experiment'])
        try:
            expname = snakemake.wildcards['experiment']
        except KeyError:
            try:
                expname = snakemake.params['expname']
            except:
                raise ValueError("No parameter for expname.")
        metadata = main(configs,
                        snakemake.params['sims'],
                        snakemake.params['reps'],
                        outdir=snakemake.params['outdir'],
                        expname=snakemake.wildcards['experiment'])
        metadata.to_csv(snakemake.output['csv'])
                