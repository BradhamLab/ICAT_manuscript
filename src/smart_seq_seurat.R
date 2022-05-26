if (exists("snakemake")) {
  .libPaths(c(snakemake@params$seurat,
              .libPaths()))
}

# check conda installation
if (exists("snakemake")) {
  # Sys.setenv(RETICULATE_PYTHON="/projectnb/bradham/dyh0110/.conda/envs/icat/bin/python")
  reticulate::use_condaenv("icat", required=TRUE)
  reticulate::py_run_string('import umap')
}

suppressPackageStartupMessages({
  library(Seurat)
  library(rjson)
  library(reticulate)
})



create_seurat <- function(X, obs, var=NULL) {
  X_data = as.data.frame(t(read.csv(X, header=FALSE)))
  obs_data = read.csv(obs, row.names=1, check.names=FALSE,
                      stringsAsFactors=FALSE)
  # drop columns with NA b/c tibble binding in Seurat breaks
  keep_cols <- which(colSums(is.na(obs_data)) == 0)
  obs_data <- data.frame(obs_data[ , keep_cols],
                         row.names=row.names(obs_data))
  names(obs_data) <- names(keep_cols)
  # force to population
  # obs_data$Population <- as.factor(obs_data$Population)
  row.names(obs_data) <- colnames(X_data)
  if (!is.null(var)) {
    var_data <- read.csv(var, row.names=1, check.names=FALSE,
                         stringsAsFactors=FALSE)
    row.names(X_data) <- row.names(var_data)
  }
  data <- Seurat::CreateSeuratObject(counts=X_data, meta.data=obs_data)
  return(data)
}

integrate_cells <- function(data, treatment, k, r) {
  by_treatment <- Seurat::SplitObject(data, split.by=treatment)
  k.filter <- 100
  print('COLUMNS: ')
  print(sapply(by_treatment, ncol))
  if (any(sapply(by_treatment, function(x) ncol(x) < 75))) {
    k.filter <- 35
  }
  anchors <- Seurat::FindIntegrationAnchors(object.list=by_treatment, dims=1:20,
                                            k.filter=k.filter)
  combined <- Seurat::IntegrateData(anchorset=anchors, dims=1:20)
  Seurat::DefaultAssay(combined) <- "integrated"
  combined <- Seurat::ScaleData(combined, verbose=FALSE)
  combined <- Seurat::RunPCA(combined, npcs=20, verbose=FALSE)
  combined <- Seurat::RunUMAP(combined, reduction='pca', dims=1:20)
  combined <- Seurat::FindNeighbors(combined, reduction='pca', dims=1:20,
                                    k.param=k)
  combined <- Seurat::FindClusters(combined, resolution=r)
  return(combined)
}

if (exists('snakemake')) {
  seurat <- create_seurat(snakemake@input[['X']],
                          snakemake@input[['obs']],
                          snakemake@input[['var']])
  r <- 1
  k <- 15
  treatment <- 'treatment'
  integrated <- integrate_cells(seurat, treatment, k, r)
  write.table(t(as.matrix(integrated@assays$integrated@data)),
              snakemake@output[['X']], sep=',',
              row.names=FALSE, col.names=FALSE)
  write.csv(integrated@meta.data, snakemake@output[['obs']])
  write.csv(integrated[["integrated"]][[]], snakemake@output[['var']])
}
