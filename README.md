# ICAT 2022 Manuscript
Repository for generating results in the manuscript "ICAT: A Novel Algorithm to
Robustly Identify Cell States Following Perturbations in Single Cell Transcriptomes"

This analysis was run on a high performance cluster with a CentOS operating
system.

It is not guaranteed to run on your system out-of-the-box.

To reproduce analysis, simply unzip the archived data files in the `archived_data.tar.gz`
file, install the `conda` environment using either the `environment.yaml` or `spec-file.txt`
files, activate the environemnt with the command, `conda activate icat`, and finally
run the pipeline using snakemake (i.e. simply execute `snakemake` in the terminal)

