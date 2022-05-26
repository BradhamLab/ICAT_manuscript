if (exists("snakemake")) {
  .libPaths(snakemake@params[['rpath']]) 
}
library(lisi)

METHOD_CLUSTERS = list(
    'icat'='sslouvain',
    'scanorama'='scanorama.louvain',
    'icat_scan'='scanorama.sslouvain',
    'seurat311'='seurat_clusters',
    'seurat_icat'='seurat.sslouvain',
    'no-int'='no.int',
    'ncfs-louvain'='ncfs.louvain'
)


assess_mixing <- function(X, obs, treatment, cluster, label, matched=TRUE) {
  # assess mixing in PCA space
  pca <- prcomp(X)$x[, c(1:50)]
  shared <- unique(obs[, label])
  # assess only shared cell-types by set intersection cell labels between
  # treatments
  if (matched) {
    for (each in unique(obs[, treatment])) {
      shared <- intersect(shared, obs[obs[, treatment] == each, label])
    }
  }
  rows <- obs[, label] %in% shared
  perplexity <- 30
  if (sum(obs[, label] %in% shared) < 30 && ~matched) {
    perplexity <- 5
  }
  lisi <- compute_lisi(pca[rows, ], obs[rows, ], treatment, perplexity=perplexity)[, 1]
  n_groups <- length(unique(obs[rows, treatment]))
  # calculate mean lisi, scale between 0 and 1
  score <- (mean(lisi) - 1) / (n_groups - 1)
  return(score)
}

get_number <- function(string) {
  return(as.numeric(gsub('[A-Z]|[a-z]', '', string)))
}

if (exists("snakemake")) {
  X <- snakemake@input[['X']]
  obs <- snakemake@input[['obs']]
  treatment <- snakemake@params[["treatment"]]
  label <- snakemake@params[["label"]]
#   cluster <- snakemake@params[['cluster']]
  if (! snakemake@params[['simulated']]) {
    X_files <- sort(X)
    obs_files <- sort(obs)
    treatment <- treatment[[snakemake@wildcards[['dataset']]]]
    label <- label[[snakemake@wildcards[['dataset']]]]
    df <- do.call(rbind, lapply(seq(length(X_files)), function(i) {
      X <- as.matrix(read.table(X_files[[i]], sep=','), stringsAsFactors=FALSE)
      obs <- read.csv(obs_files[[i]], row.names=1, stringsAsFactors=FALSE)
      method <- unlist(strsplit(X_files[[i]], '/'))[[4]]
      out <- data.frame(list(dataset=snakemake@wildcards[['dataset']],
                             method=method))
          
      print(paste0(method, ' assessing mixture...'))
      out['LISI'] <- assess_mixing(X, obs, treatment,
                                   METHOD_CLUSTERS[[method]],
                                   label)
      return(out)
    }))
  } else {
      filesplit <- unlist(strsplit(X, '/'))
      method <- filesplit[4]
      X <- as.matrix(read.table(X, sep=','), stringsAsFactors=FALSE)
      obs <- read.csv(obs, row.names=1, stringsAsFactors=FALSE)
      df <- data.frame(list(method=method,
                             experiment=filesplit[5],
                             sim=get_number(substr(filesplit[6], 1, 4)),
                             rep=get_number(substr(filesplit[6], 5, 8))))

      print(paste0(method, ' assessing mixture...'))
      df['LISI'] <- assess_mixing(X, obs, treatment,
                                  METHOD_CLUSTERS[[method]], label)
      df['Activated.LISI'] = NA
      # assess separation between activated and non-activate cell types
      if ('C1+' %in% obs[, label]) {
          rows = obs[, label] %in% c('C1', 'C1+')
          # pass label as treatment to assess separation between C1 and C1+
          # labels
          df['Activated.LISI'] = assess_mixing(X[rows, ], obs[rows, ],
                                               label, METHOD_CLUSTERS[[method]],
                                               label, matched=FALSE)
      }
  }
  write.csv(df, snakemake@output[['csv']])
}