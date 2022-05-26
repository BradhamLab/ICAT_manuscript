import pandas as pd

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        dfs = []
        for each in snakemake.input['perf']:
            df = pd.read_csv(each, index_col=0)
            df.rename(columns={'dataset': 'Dataset'},
                      inplace=True)
            dfs.append(df)
        out = pd.concat(dfs)
        out.to_csv(snakemake.output['csv'])