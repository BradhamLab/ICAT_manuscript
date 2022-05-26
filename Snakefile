import glob
import os
import itertools


configfile: "./snakemake-config.yaml"


FILES = ["X.csv", "obs.csv", "var.csv"]

METHODS = ["icat", "seurat311", "scanorama", "icat_scan", "seurat_icat", "no-int"]

(EXPERIMENTS,) = glob_wildcards("data/external/simulations/{experiment}_params.json")

SIMS = [f"Sim{i}" for i in range(1, int(config["simulations"]["sims"]) + 1)]
REPS = [f"Rep{i}" for i in range(1, int(config["simulations"]["reps"]) + 1)]
RUNS = [
    f"{exp}/{sim}{rep}" for exp, sim, rep in itertools.product(EXPERIMENTS, SIMS, REPS)
]
SIMULATED = [
    "data/processed/simulated/{run}/{out}".format(run=run, out=out)
    for run, out in itertools.product(RUNS, FILES)
]

BENCHMARK = ["cellmix1", "cellmix2", "cellmix3", "cellmix4", "sc_celseq2"]
MIXES = BENCHMARK[:-1]

DATASETS = ["Tian", "Kang", "Kagohara"]


rule all:
    input:
        "data/results/simulated/final/results.csv",
        "data/results/subsample.csv",
        expand(
            os.path.join("data", "results", "{dataset}", "results.csv"),
            dataset=DATASETS,
        ),


# ---------------------------- Generate Simulated Data -------------------------
rule simulate_data:
    input:
        json=os.path.join("data", "external", "simulations", "{experiment}_params.json"),
    params:
        sims=config["simulations"]["sims"],
        reps=config["simulations"]["reps"],
        outdir=os.path.join("data", "processed", "simulated", "{experiment}"),
    output:
        data=expand(
            os.path.join(
                "data",
                "processed",
                "simulated",
                "{{experiment}}",
                "{sim}{rep}",
                "{type}.csv",
            ),
            sim=SIMS,
            rep=REPS,
            type=["X", "obs", "var"],
        ),
        csv=os.path.join(
            "data", "processed", "simulated", "{experiment}", "simulations.csv"
        ),
    script:
        "src/generate_simulated_datasets.py"


# ---------------------- Fit and Analyze Simulated Data ------------------------
rule fit_simulated:
    input:
        X=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "X.csv"
        ),
        obs=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "obs.csv"
        ),
        var=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "var.csv"
        ),
    params:
        treatment="Treatment",
        control="Control",
        label="Population",
        simulated=True,
    output:
        json="data/interim/fits/simulated/{experiment}/{sim}{rep}_fit.json",
    script:
        "src/fit_louvain.py"


rule simulated_louvain:
    input:
        X=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "X.csv"
        ),
        obs=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "obs.csv"
        ),
        var=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "var.csv"
        ),
        json="data/interim/fits/simulated/{experiment}/{sim}{rep}_fit.json",
    params:
        outdir="data/results/simulated/no-int/{experiment}/{sim}{rep}/",
        cluster="no-int",
    output:
        X="data/results/simulated/no-int/{experiment}/{sim}{rep}/X.csv",
        obs="data/results/simulated/no-int/{experiment}/{sim}{rep}/obs.csv",
        var="data/results/simulated/no-int/{experiment}/{sim}{rep}/var.csv",
    script:
        "src/run_louvain.py"


rule simulated_icat:
    input:
        X=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "X.csv"
        ),
        obs=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "obs.csv"
        ),
        var=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "var.csv"
        ),
        json="data/interim/fits/simulated/{experiment}/{sim}{rep}_fit.json",
        icat="data/external/icat_params.json",
    params:
        treatment="Treatment",
        controls="Control",
        plotdir="reports/figures/simulated/{experiment}/icat/{sim}{rep}",
        outdir="data/results/simulated/icat/{experiment}/{sim}{rep}/",
        simulated=True,
        cluster="sslouvain",
    output:
        X="data/results/simulated/icat/{experiment}/{sim}{rep}/X.csv",
        obs="data/results/simulated/icat/{experiment}/{sim}{rep}/obs.csv",
        var="data/results/simulated/icat/{experiment}/{sim}{rep}/var.csv",
    script:
        "src/run_icat.py"


rule simulated_seurat311:
    input:
        X=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "X.csv"
        ),
        obs=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "obs.csv"
        ),
        json="data/interim/fits/simulated/{experiment}/{sim}{rep}_fit.json",
    params:
        treatment="Treatment",
        seurat=config["libraries"]["seurat3.1.1"],
        python=config["libraries"]["conda_env"],
        simulated=True,
    output:
        X="data/results/simulated/seurat311/{experiment}/{sim}{rep}/X.csv",
        obs="data/results/simulated/seurat311/{experiment}/{sim}{rep}/obs.csv",
    script:
        "src/run_seurat3-1-1.R"


rule simulated_seurat_icat:
    input:
        X="data/results/simulated/seurat311/{experiment}/{sim}{rep}/X.csv",
        obs="data/results/simulated/seurat311/{experiment}/{sim}{rep}/obs.csv",
        var="data/processed/simulated/{experiment}/{sim}{rep}/var.csv",
        json="data/interim/fits/simulated/{experiment}/{sim}{rep}_fit.json",
        icat="data/external/icat_params.json",
    params:
        treatment="Treatment",
        controls="Control",
        outdir="data/results/simulated/seurat_icat/{experiment}/{sim}{rep}",
        cluster="seurat.sslouvain",
        simulated=True,
    output:
        X="data/results/simulated/seurat_icat/{experiment}/{sim}{rep}/X.csv",
        obs="data/results/simulated/seurat_icat/{experiment}/{sim}{rep}/obs.csv",
    script:
        "src/run_icat.py"


rule simulated_scanorama:
    input:
        X=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "X.csv"
        ),
        obs=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "obs.csv"
        ),
        var=os.path.join(
            "data", "processed", "simulated", "{experiment}", "{sim}{rep}", "var.csv"
        ),
        json=os.path.join(
            "data",
            "interim",
            "fits",
            "simulated",
            "{experiment}",
            "{sim}{rep}_fit.json",
        ),
    output:
        X=os.path.join(
            "data",
            "results",
            "simulated",
            "scanorama",
            "{experiment}",
            "{sim}{rep}",
            "X.csv",
        ),
        obs=os.path.join(
            "data",
            "results",
            "simulated",
            "scanorama",
            "{experiment}",
            "{sim}{rep}",
            "obs.csv",
        ),
        var=os.path.join(
            "data",
            "results",
            "simulated",
            "scanorama",
            "{experiment}",
            "{sim}{rep}",
            "var.csv",
        ),
    params:
        simulated=True,
        treatment="Treatment",
        controls="Control",
        outdir=os.path.join(
            "data", "results", "simulated", "scanorama", "{experiment}", "{sim}{rep}"
        ),
    script:
        "src/run_scanorama.py"


rule simulated_scanorama_icat:
    input:
        X=os.path.join(
            "data",
            "results",
            "simulated",
            "scanorama",
            "{experiment}",
            "{sim}{rep}",
            "X.csv",
        ),
        obs=os.path.join(
            "data",
            "results",
            "simulated",
            "scanorama",
            "{experiment}",
            "{sim}{rep}",
            "obs.csv",
        ),
        var=os.path.join(
            "data",
            "results",
            "simulated",
            "scanorama",
            "{experiment}",
            "{sim}{rep}",
            "var.csv",
        ),
        json=os.path.join(
            "data",
            "interim",
            "fits",
            "simulated",
            "{experiment}",
            "{sim}{rep}_fit.json",
        ),
        icat="data/external/icat_params.json",
    output:
        X=os.path.join(
            "data",
            "results",
            "simulated",
            "icat_scan",
            "{experiment}",
            "{sim}{rep}",
            "X.csv",
        ),
        obs=os.path.join(
            "data",
            "results",
            "simulated",
            "icat_scan",
            "{experiment}",
            "{sim}{rep}",
            "obs.csv",
        ),
    params:
        outdir=os.path.join(
            "data", "results", "simulated", "icat_scan", "{experiment}", "{sim}{rep}"
        ),
        treatment="Treatment",
        controls="Control",
        cluster="scanorama.sslouvain",
        simulated=True,
    script:
        "src/run_icat.py"


rule summarize_simulated_data:
    input:
        json=expand(
            os.path.join(
                "data", "external", "simulations", "{experiment}_params.json"
            ),
            experiment=EXPERIMENTS,
        ),
        var=expand(
            os.path.join(
                "data",
                "processed",
                "simulated",
                "{experiment}",
                "{sim}{rep}",
                "var.csv",
            ),
            experiment=EXPERIMENTS,
            sim=SIMS,
            rep=REPS,
        ),
    output:
        csv=os.path.join("data", "results", "simulated", "final", "dataset_meta.csv"),
    script:
        "src/summarize_datasets.py"


rule evaluate_treatment_integration:
    input:
        X="data/results/simulated/{method}/{experiment}/{sim}{rep}/X.csv",
        obs="data/results/simulated/{method}/{experiment}/{sim}{rep}/obs.csv",
    params:
        treatment="Treatment",
        simulated=True,
        label="Population",
        rpath=config["libraries"]["lisi"],
    output:
        csv="data/results/simulated/{method}/{experiment}/{sim}{rep}/lisi.csv",
    script:
        "src/calculate_lisi.R"


rule evaluate_methods_simulated:
    input:
        xs=expand(
            "data/results/simulated/{method}/{experiment}/{sim}{rep}/X.csv",
            method=METHODS,
            experiment=EXPERIMENTS,
            sim=SIMS,
            rep=REPS,
        ),
        obss=expand(
            "data/results/simulated/{method}/{experiment}/{sim}{rep}/obs.csv",
            method=METHODS,
            experiment=EXPERIMENTS,
            sim=SIMS,
            rep=REPS,
        ),
        lisi=expand(
            "data/results/simulated/{method}/{experiment}/{sim}{rep}/lisi.csv",
            method=METHODS,
            experiment=EXPERIMENTS,
            sim=SIMS,
            rep=REPS,
        ),
        meta=os.path.join("data", "results", "simulated", "final", "dataset_meta.csv"),
    output:
        csv="data/results/simulated/final/results.csv",
    params:
        identity="Population",
    script:
        "src/evaluate_clusters.py"


rule summarize_simulated_performance:
    input:
        csv="data/results/simulated/final/results.csv",
    output:
        csv="data/results/simulated/final/formatted_results.csv",
    script:
        "src/format_results.py"


# --------------------------- Process Tian Data ---------------------------


# normalize across cells, subset to highly variable genes, and ln(X +1)
# transform cells
rule format_tian_data:
    input:
        counts=[
            "data/raw/Tian/{bench}.count.csv".format(bench=bench)
            for bench in BENCHMARK
        ],
        meta=[
            "data/raw/Tian/{bench}.metadata.csv".format(bench=bench)
            for bench in BENCHMARK
        ],
    params:
        outdir="data/processed/Tian/hvgs",
        plotdir="figures/Tian/",
    output:
        X="data/processed/Tian/hvgs/X.csv",
        obs="data/processed/Tian/hvgs/obs.csv",
        var="data/processed/Tian/hvgs/var.csv",
    script:
        "src/format_preprocess_Tian.py"


# ---------------------------- Format Kang Data -------------------------------
rule format_kang_data:
    input:
        mtx=["data/raw/Kang/GSM2560248_2.1.mtx", "data/raw/Kang/GSM2560249_2.2.mtx"],
        barcodes=[
            "data/raw/Kang/GSM2560248_barcodes.tsv",
            "data/raw/Kang/GSM2560249_barcodes.tsv",
        ],
        genes="data/raw/Kang/GSE96583_batch2.genes.tsv",
        cells="data/raw/Kang/GSE96583_batch2.total.tsne.df.tsv",
        filter_cells="data/external/remove_cells.txt",
    params:
        datadir="data/raw/Kang/",
        outdir="data/processed/Kang/hvgs/",
    output:
        "data/processed/Kang/hvgs/X.csv",
        "data/processed/Kang/hvgs/obs.csv",
        "data/processed/Kang/hvgs/var.csv",
    script:
        "src/format_preprocess_Kang.py"


kago_input = [
    "GSE137524_exprsSCC1Matrix.csv",
    "GSE137524_exprsSCC6Matrix.csv",
    "GSE137524_phenoDataSCC1.csv",
    "GSE137524_phenoDataSCC6.csv",
    "GSE137524_exprsSCC25Matrix.csv",
    "GSE137524_featureData.csv",
    "GSE137524_phenoDataSCC25.csv",
]


rule format_kagohara_data:
    input:
        [os.path.join("data/raw/Kagohara", x) for x in kago_input],
    output:
        "data/processed/Kagohara/hvgs/X.csv",
        "data/processed/Kagohara/hvgs/obs.csv",
        "data/processed/Kagohara/hvgs/var.csv",
    params:
        datadir="data/raw/Kagohara/",
        outdir="data/processed/Kagohara/hvgs/",
    script:
        "src/format_preprocess_Kagohara.py"


rule fit_real_data:
    input:
        X="data/processed/{dataset}/hvgs/X.csv",
        obs="data/processed/{dataset}/hvgs/obs.csv",
        var="data/processed/{dataset}/hvgs/var.csv",
    params:
        treatment=config["params"]["treatment"],
        control=config["params"]["control"],
        label=config["params"]["label"],
        simulated=False,
    output:
        json="data/interim/fits/{dataset}/isolated_fits.json",
    script:
        "src/fit_louvain.py"


rule no_integration:
    input:
        X="data/processed/{dataset}/hvgs/X.csv",
        obs="data/processed/{dataset}/hvgs/obs.csv",
        var="data/processed/{dataset}/hvgs/var.csv",
        json="data/interim/fits/{dataset}/isolated_fits.json",
    params:
        outdir="data/results/{dataset}/no-int/",
        cluster="no-int",
    output:
        X="data/results/{dataset}/no-int/X.csv",
        obs="data/results/{dataset}/no-int/obs.csv",
        var="data/results/{dataset}/no-int/var.csv",
    script:
        "src/run_louvain.py"


rule run_icat:
    input:
        X="data/processed/{dataset}/hvgs/X.csv",
        obs="data/processed/{dataset}/hvgs/obs.csv",
        var="data/processed/{dataset}/hvgs/var.csv",
        json="data/interim/fits/{dataset}/isolated_fits.json",
        icat="data/external/{dataset}_icat_params.json",
    params:
        treatment=config["params"]["treatment"],
        controls=config["params"]["control"],
        outdir="data/results/{dataset}/icat/",
        cluster="sslouvain",
        simulated=False,
    output:
        X="data/results/{dataset}/icat/X.csv",
        obs="data/results/{dataset}/icat/obs.csv",
        var="data/results/{dataset}/icat/var.csv",
    script:
        "src/run_icat.py"


rule run_seurat311:
    input:
        X="data/processed/{dataset}/hvgs/X.csv",
        obs="data/processed/{dataset}/hvgs/obs.csv",
        var="data/processed/{dataset}/hvgs/var.csv",
        json="data/interim/fits/{dataset}/isolated_fits.json",
    params:
        treatment=config["params"]["treatment"],
        seurat=config["libraries"]["seurat3.1.1"],
        python=config["libraries"]["conda_env"],
        simulated=False,
    output:
        X="data/results/{dataset}/seurat311/X.csv",
        obs="data/results/{dataset}/seurat311/obs.csv",
        var="data/results/{dataset}/seurat311/var.csv",
    script:
        "src/run_seurat3-1-1.R"


rule run_seurat_icat:
    input:
        X="data/results/{dataset}/seurat311/X.csv",
        obs="data/results/{dataset}/seurat311/obs.csv",
        var="data/results/{dataset}/seurat311/var.csv",
        json="data/interim/fits/{dataset}/isolated_fits.json",
        icat="data/external/{dataset}_icat_params.json",
    params:
        outdir="data/results/{dataset}/seurat_icat/",
        treatment=config["params"]["treatment"],
        controls=config["params"]["control"],
        cluster="seurat.sslouvain",
        simulated=False,
    output:
        X="data/results/{dataset}/seurat_icat/X.csv",
        obs="data/results/{dataset}/seurat_icat/obs.csv",
    script:
        "src/run_icat.py"


rule run_scanorama:
    input:
        X="data/processed/{dataset}/hvgs/X.csv",
        obs="data/processed/{dataset}/hvgs/obs.csv",
        var="data/processed/{dataset}/hvgs/var.csv",
        json="data/interim/fits/{dataset}/isolated_fits.json",
    params:
        treatment=config["params"]["treatment"],
        controls=config["params"]["control"],
        outdir="data/results/{dataset}/scanorama/",
        simulated=False,
    output:
        X="data/results/{dataset}/scanorama/X.csv",
        obs="data/results/{dataset}/scanorama/obs.csv",
        var="data/results/{dataset}/scanorama/var.csv",
    script:
        "src/run_scanorama.py"


rule run_scanorama_icat:
    input:
        X="data/results/{dataset}/scanorama/X.csv",
        obs="data/results/{dataset}/scanorama/obs.csv",
        var="data/results/{dataset}/scanorama/var.csv",
        json="data/interim/fits/{dataset}/isolated_fits.json",
        icat="data/external/{dataset}_icat_params.json",
    output:
        X="data/results/{dataset}/icat_scan/X.csv",
        obs="data/results/{dataset}/icat_scan/obs.csv",
        var="data/results/{dataset}/icat_scan/var.csv",
    params:
        outdir="data/results/{dataset}/icat_scan/",
        treatment=config["params"]["treatment"],
        controls=config["params"]["control"],
        cluster="scanorama.sslouvain",
        simulated=False,
    script:
        "src/run_icat.py"


rule test_subsample:
    input:
        X="data/processed/{large_ds}/hvgs/X.csv",
        obs="data/processed/{large_ds}/hvgs/obs.csv",
        fit="data/interim/fits/{large_ds}/isolated_fits.json",
        icat="data/external/{large_ds}_icat_params.json",
    output:
        csv=os.path.join(
            "data", "results", "{large_ds}", "subsample", "{cells}_performance.csv"
        ),
        X=os.path.join("data", "results", "{large_ds}", "subsample", "{cells}", "X.csv"),
        obs=os.path.join(
            "data", "results", "{large_ds}", "subsample", "{cells}", "obs.csv"
        ),
        var=os.path.join(
            "data", "results", "{large_ds}", "subsample", "{cells}", "var.csv"
        ),
    params:
        treatment=config["params"]["treatment"],
        controls=config["params"]["control"],
        label=config["params"]["label"],
        outdir=os.path.join("data", "results", "{large_ds}", "subsample", "{cells}"),
    script:
        "src/run_icat_subsamples.py"


rule combine_subsample_results:
    input:
        perf=expand(
            os.path.join(
                "data", "results", "{large_ds}", "subsample", "{cells}_performance.csv"
            ),
            large_ds=["Kagohara", "Kang"],
            cells=config["params"]["subsample"]["cells"],
        ),
    output:
        csv=os.path.join("data", "results", "subsample.csv"),
    script:
        "src/combine_subsample_results.py"


rule dataset_integration:
    input:
        X=expand(
            os.path.join("data", "results", "{{dataset}}", "{method}", "X.csv"),
            method=METHODS,
        ),
        obs=expand(
            os.path.join("data", "results", "{{dataset}}", "{method}", "obs.csv"),
            method=METHODS,
        ),
    params:
        treatment=config["params"]["treatment"],
        simulated=False,
        label=config["params"]["label"],
        rpath=config["libraries"]["lisi"],
    output:
        csv=os.path.join("data", "results", "{dataset}", "lisi.csv"),
    script:
        "src/calculate_lisi.R"


rule summarize_dataset:
    input:
        X=expand(
            os.path.join("data", "results", "{{dataset}}", "{method}", "X.csv"),
            method=METHODS,
        ),
        obs=expand(
            os.path.join("data", "results", "{{dataset}}", "{method}", "obs.csv"),
            method=METHODS,
        ),
        lisi=os.path.join("data", "results", "{dataset}", "lisi.csv"),
        fit=os.path.join("data", "interim", "fits", "{dataset}", "isolated_fits.json"),
    params:
        identity=config["params"]["label"],
        treatment=config["params"]["treatment"],
        controls=config["params"]["control"],
        plotdir=os.path.join("reports", "figures", "{dataset}"),
    output:
        csv="data/results/{dataset}/results.csv",
        metrics="reports/figures/{dataset}/metrics.svg",
        lisi="reports/figures/{dataset}/lisi.svg",
        radar="reports/figures/{dataset}/spiderplot.svg",
        obs=[
            os.path.join("data", "results", "{dataset}", "obs", f"{x}_cells.csv")
            for x in METHODS
        ],
        svg=[
            os.path.join("reports", "figures", "{dataset}", x, "umap_clusters.png")
            for x in METHODS
        ],
    script:
        "src/summarize_dataset.py"


rule icat_smart_seq:
    input:
        counts="data/raw/smart-seq/normalized_log_matrix.csv",
        var="data/raw/smart-seq/evms_w_length.csv",
        obs="data/raw/smart-seq/filtered_metadata.csv",
    output:
        adata="data/processed/smart-seq/hvgs.h5ad",
        X="data/processed/smart-seq/X.csv",
        obs="data/processed/smart-seq/obs.csv",
        var="data/processed/smart-seq/var.csv",
        gtest="data/results/smart-seq/icat_gtest.csv",
        fisher="data/results/smart-seq/icat_fisher.csv",
    script:
        "src/smart_seq_icat.py"


rule scanorama_smart_seq:
    input:
        X="data/processed/smart-seq/X.csv",
        obs="data/processed/smart-seq/obs.csv",
        var="data/processed/smart-seq/var.csv",
    params:
        treatment="treatment",
        controls="Control",
        outdir="data/processed/smart-seq/scanorama",
    output:
        adata="data/processed/smart-seq/scanorama/adata.h5ad",
    script:
        "src/smart_seq_scan.py"


rule seurat_smart_seq:
    input:
        X="data/processed/smart-seq/X.csv",
        obs="data/processed/smart-seq/obs.csv",
        var="data/processed/smart-seq/var.csv",
    params:
        treatment="treatment",
    output:
        X="data/processed/smart-seq/seurat/X.csv",
        obs="data/processed/smart-seq/seurat/obs.csv",
        var="data/processed/smart-seq/seurat/var.csv",
    script:
        "src/smart_seq_seurat.R"


rule seurat_to_adata:
    input:
        X="data/processed/smart-seq/seurat/X.csv",
        obs="data/processed/smart-seq/seurat/obs.csv",
        var="data/processed/smart-seq/seurat/var.csv",
    output:
        adata="data/processed/smart-seq/seurat/adata.h5ad",
    script:
        "src/seurat_to_h5.py"
