library(SingleCellExperiment)
library(scran)

create_sce <- function(X, obs, var) {
    X <- read.csv(X, row.names=NULL, header=FALSE)
    obs <- read.csv(obs, row.names=1,
                    stringsAsFactors=FALSE)
    var <- read.csv(var, row.names=1,
                    stringsAsFactors=FALSE)
    row.names(X) <- row.names(obs)
    colnames(X) <- row.names(var)
    sce <- SingleCellExperiment(list("counts"=t(X)),
                                colData=obs,
                                rowData=var)
    return(sce)
}

preprocess_by_treatment <- function(sce, treatment) {
    normalized <- lapply(unique(colData(sce)[ , treatment]), function(x) {
        by_treatment <- sce[, colData(sce)[, treatment]==x]
        clusters <- scran::quickCluster(by_treatment)
        by_treatment <- scran::computeSumFactors(by_treatment,
                                                 clusters=clusters)
        return(normalize(by_treatment))
    })
    return(do.call(cbind, normalized))
}


if (exists('snakemake')) {
    sce <- create_sce(snakemake@input[['X']],
                      snakemake@input[['obs']],
                      snakemake@input[['var']])
    treatment <- snakemake@params[['treatment']][[snakemake@wildcards[['dataset']]]]
    print(paste0("Normalzing counts separated by ", treatment, "."))
    sce <- preprocess_by_treatment(sce, treatment)
    print(paste0("Normalized single cell experiment table:", "\n"))
    print(sce)
    write.table(t(SingleCellExperiment::logcounts(sce)),
                snakemake@output[['X']], sep=',',
                row.names=FALSE, col.names=FALSE)
    write.csv(colData(sce), snakemake@output[['obs']])
    write.csv(rowData(sce), snakemake@output[['var']])
}