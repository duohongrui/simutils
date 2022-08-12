#' Calculate Batch Quantifying Metrics
#'
#' @param data A matrix of gene expression profile.
#' @param batch_info Batch assignment of every cells in columns of matrix.
#' @param cluster_info Cluster(or group) assignment of every cells in columns of matrix.
#' @param verbose Whether messages are returned during the process.
#' @importFrom SingleCellExperiment SingleCellExperiment logcounts reducedDim colData reducedDim<-
#' @importFrom scater logNormCounts
#' @importFrom BiocGenerics var
#' @importFrom stats prcomp
#' @export
calculate_batch_properties <- function(
  data,
  batch_info,
  cluster_info = NULL,
  verbose = FALSE
){
  ## 1) cms, 2) iLISI, 3) Mixing metric, 4) Shannon entropy
  if(!requireNamespace("CellMixS", quietly = TRUE)){
    stop("Package \"CellMixS\" must be installed by \"BiocManager::install('CellMixS')\" command.")
  }
  colData <- data.frame("batch" = batch_info)
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = data),
                                                    colData = colData)
  sce <- scater::logNormCounts(sce)
  var_gene <- apply(data, 1, BiocGenerics::var)
  sort_index <- sort(var_gene, decreasing = TRUE)[1:2000]
  pca <- stats::prcomp(t(SingleCellExperiment::logcounts(sce)[names(sort_index), ]),
                       center = TRUE,
                       scale. = FALSE,
                       rank. = 50)
  SingleCellExperiment::reducedDim(sce, "PCA") <- pca$x
  if(verbose){
    message("Calculate cms, ILSI, Mixing Metric and Shannon entropy...")
  }
  result <- CellMixS::evalIntegration(metrics = c("cms",
                                                  "isi",
                                                  "mixingMetric",
                                                  "entropy"),
                                      sce = sce,
                                      group = "batch",
                                      k = 5)
  col_result <- as.data.frame(SingleCellExperiment::colData(result))
  cms <- mean(col_result$cms)
  LISI <- mean(col_result$isi)
  mm <- mean(col_result$mm)
  shannon_entropy <- mean(col_result$entropy)
  ## 5) kBET
  if(!requireNamespace("kBET", quietly = TRUE)){
    stop("Package \"kBET\" must be installed by \"devtools::install_github('theislab/kBET')\" command.")
  }
  if(verbose){
    message("Calculate kBET...")
  }
  batch_estimate <- kBET::kBET(df = t(data),
                               batch_info,
                               plot = FALSE)
  rejection <- batch_estimate[["results"]]
  kBET <- mean(rejection[, 1])
  ## 6) Average Silouette Width (ASW for batch)
  if(verbose){
    message("Calculate Average Silouette Width for batch..")
  }
  AWS_batch <- kBET::batch_sil(pca.data = pca, batch = as.factor(batch_info))
  ## 7) Principal Component Regression (PCR)
  if(verbose){
    message("Calculate principal component regression...")
  }
  batch_pca <- kBET::pcRegression(pca.data = pca,
                                  batch = batch_info,
                                  n_top = 50)
  pcr <- batch_pca$R2Var
  ## 8) Graph connectivity
  if(!requireNamespace("bluster", quietly = TRUE)){
    stop("Package \"bluster\" must be installed by \"BiocManager::install('bluster')\" command.")
  }
  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("Package \"igraph\" must be installed by \"install.packages('igraph')\" command.")
  }
  if(is.null(cluster_info)){
    gc <- NULL
  }else{
    if(verbose){
      message("Calculate graph connectivity...")
    }
    cluster <- unique(cluster_info)
    gc <- c()
    for(i in 1:length(cluster)){
      cluster_name <- cluster[i]
      index <- cluster_info == cluster_name
      ## Filter data
      sub_data <- SingleCellExperiment::logcounts(sce)[names(sort_index), index]
      ## Perform PCA
      sub_pca <- stats::prcomp(t(sub_data),
                               center = TRUE,
                               scale. = FALSE,
                               rank. = 50)
      ## Make KNN graph
      knn_graph <- bluster::makeKNNGraph(sub_pca$x, k = 5)
      ## Generate connected component graph
      com_graph <- igraph::components(knn_graph)
      new_gc_cluster <- max(com_graph$csize)/sum(com_graph$csize)
      gc <- append(gc, new_gc_cluster)
    }
    gc <- sum(gc)/length(cluster)
  }
  dplyr::lst(cms,
             LISI,
             mm,
             shannon_entropy,
             kBET,
             AWS_batch,
             pcr,
             gc)
}

