#' Calculate Batch Quantifying Metrics
#'
#' @param data A matrix of gene expression profile.
#' @param batch_info Batch assignment of every cells in columns of matrix.
#' @param k k nearest neighborhoods of every cell.
#' @param verbose Whether messages are returned during the process.
#' @importFrom SingleCellExperiment SingleCellExperiment logcounts reducedDim colData reducedDim<-
#' @importFrom scater logNormCounts
#' @importFrom BiocGenerics var
#' @importFrom stats prcomp
#' @importFrom BiocManager install
#' @export
calculate_batch_properties <- function(
  data,
  batch_info,
  k,
  verbose = TRUE
){
  ## 1) cms, 2) ISI, 3) Mixing metric, 4) Shannon entropy
  if(!requireNamespace("CellMixS", quietly = TRUE)){
    message("Install CellMixS...")
    BiocManager::install('CellMixS')
  }
  colData <- data.frame("batch" = batch_info)
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = data),
                                                    colData = colData)
  sce <- scater::logNormCounts(sce)
  var_gene <- apply(data, 1, BiocGenerics::var)
  if(length(var_gene) >= 2000){
    variable_num <- 2000
  }else{
    variable_num <- length(var_gene)/2
  }
  sort_index <- sort(var_gene, decreasing = TRUE)[1:variable_num]
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
                                      k = k)
  col_result <- as.data.frame(SingleCellExperiment::colData(result))
  cms <- col_result$cms
  ISI <- col_result$isi
  mm <- col_result$mm
  shannon_entropy <- col_result$entropy
  ## 5) kBET
  if(!requireNamespace("kBET", quietly = TRUE)){
    message("Install kBET...")
    devtools::install_github('theislab/kBET')
  }
  if(verbose){
    message("Calculate kBET...")
  }
  batch_estimate <- kBET::kBET(df = t(data),
                               batch_info,
                               plot = FALSE,
                               k0 = k)
  rejection <- batch_estimate[["results"]]
  kBET <- rejection[, 1]
  ## 6) Average Silouette Width (ASW for batch)
  if(verbose){
    message("Calculate Average Silouette Width for batch..")
  }
  AWS_batch <- kBET::batch_sil(pca.data = pca,
                               batch = as.factor(batch_info),
                               nPCs = 50)
  ## 7) Principal Component Regression (PCR)
  if(verbose){
    message("Calculate principal component regression... \n")
  }
  batch_pca <- kBET::pcRegression(pca.data = pca,
                                  batch = batch_info,
                                  n_top = 50)
  pcr <- batch_pca$R2Var
  ## 8) LISI
  if(!requireNamespace("lisi", quietly = TRUE)){
    message("Install lisi...")
    devtools::install_github("immunogenomics/lisi")
  }
  LISI <- lisi::compute_lisi(pca$x,
                             meta_data = colData,
                             label_colnames = "batch",
                             perplexity = k)
  LISI <- mean(LISI$batch)
  dplyr::lst(cms,
             ISI,
             LISI,
             mm,
             shannon_entropy,
             kBET,
             AWS_batch,
             pcr)
}

