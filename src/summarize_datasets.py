import json 
import re
import itertools
import os

import numpy as np
import pandas as pd

def get_run_data(filepath):
    directory = os.path.split(filepath)[0]
    path_split = directory.split(os.path.sep)
    return({'run': "{}{}".format(path_split[-2], path_split[-1]),
            'experiment': path_split[-2]})

def gini(pops):
    value = 0
    for x, y in itertools.product(pops, pops):
        value += abs(x - y)
    return value / (2 * len(pops) ** 2 * np.mean(pops))

def summarize_var(var):
    marker_cols = [x for x in var.columns if 'Marker' in x]
    base_col = 'Base.Mu'
    perturb_col = 'Perturbation.Shift'
    n_markers = 0
    marker_shift = 0
    perturb_shift = (var[perturb_col] - 1).sum()
    for col in marker_cols:
        exp_col = col.replace('Marker', 'Mu')
        n_markers += var[col].sum()
        # marker shifts aren't included for stimulated populations since they
        # are the same -- avoid double counting
        try:
            marker_shift += abs(var[exp_col][var[col]]\
                                - var[base_col][var[col]]).sum()
        except KeyError:
            pass
    return {"Markers (Total)": n_markers,
            "Markers (%)": n_markers / var.shape[0] * 100,
            "Total Average Marker Shift": marker_shift,
            "Perturbation Shift": perturb_shift}

def population_data(obs):
    ctrl = obs[obs['Treatment'] == 'Control'].groupby('Population').size()
    prtb = obs[obs['Treatment'] != 'Control'].groupby('Population').size()
    combined = ctrl.copy()
    for each in prtb.index.values:
        if each in combined.index.values:
            combined[each] += prtb[each]
        else:
            combined[each] = prtb[each]
    added = pd.concat([ctrl, prtb])
    return {"Gini (Control)": gini(ctrl.values),
            "Gini (Perturbed)": gini(prtb.values),
            "Gini (Merged)": gini(combined.values),
            "Gini (Combined)": gini(added.values)}

if __name__ == "__main__":
    dfs = []
    for filepath in snakemake.input['var']:
        obs = pd.read_csv(filepath.replace('var.csv', 'obs.csv'), index_col=0)
        pop_data = population_data(obs)
        data = pd.read_csv(filepath, index_col=0)
        var_data = summarize_var(data)
        run_data = get_run_data(filepath)
        for each in snakemake.input['json']:
            if run_data['experiment'] in each:
                with open(each) as f:
                    exp_data = json.load(f)
                    var_data['Percent Perturbed'] = exp_data['perturb_kwargs']["percent_perturb"]
                    var_data['Markers (Mean)'] = data.shape[0] * exp_data['control_kwargs']['p_marker']
        var_data['Run'] = run_data['run']
        var_data.update(pop_data)
        dfs.append(pd.DataFrame(var_data, index=[run_data['run']]))
    out_df = pd.concat(dfs, axis=0, sort=False).reset_index(drop=True)
    out_df.to_csv(snakemake.output['csv'])
        
