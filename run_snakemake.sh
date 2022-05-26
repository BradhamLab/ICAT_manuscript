#! /bin/bash -l
#$ -N evaluate_icat
#$ -M dyh0110@bu.edu
#$ -m eas
#$ -l h_rt=24:00:00

export PYTHONPATH="{$PYTHONPATH}:/projectnb/bradham/PythonModules"
module load R/3.6.2;
export R_LIBS="{$R_LIBS}:/projectnb/bradham/RPackages/3.6.2"
source activate icat;
snakemake data/processed/smart-seq/hvgs.h5ad data/processed/smart-seq/seurat/adata.h5ad data/processed/smart-seq/scanorama/adata.h5ad --force \
--cluster "qsub -v PYTHONPATH={cluster.python} -v LD_LIBRARY_PATH={cluster.gcc} -P {cluster.project} -pe omp {cluster.processors} -e {cluster.error} -o {cluster.out} -l cpu_arch={cluster.cpus} -l mem_per_core={cluster.memory} -v NUMBA_NUM_THREADS={cluster.processors} -l h_rt={cluster.runtime}" --jobs 100 --latency-wait 60 --cluster-config cluster.json
