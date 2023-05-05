#' Calculate Cell Properties
#'
#' @param data A count matrix.
#' @param verbose Whether messages are returned during the process.
#' @importFrom edgeR cpm DGEList calcNormFactors
#' @importFrom dplyr lst
#' @importFrom stats IQR
#' @return A list
#' @export
#' @examples
#' data <- matrix(rpois(100*100, 2),
#'                nrow = 100,
#'                dimnames = list(paste0('Gene', 1:100),
#'                                paste0('Cell', 1:100)))
#' result <- cell_properties(data, verbose = TRUE)
cell_properties <- function(data,
                            verbose = FALSE){
  ## 1) library size
  if(verbose){
    message("Calculating library size of cells...")
  }
  library_size <- colSums(data)
  ## 2) fraction of zero in cells
  if(verbose){
    message("Calculating fraction of zero in cells...")
  }
  zero_fraction_cell <- apply(data, 2, function(x){
    sum(x == 0)/length(x)
  })
  ## 3) cell correlation
  if(verbose){
    message("Calculating cell correlation...")
  }
  if(!requireNamespace("WGCNA", quietly = TRUE)){
    message("WGCNA is not installed on your device...")
    message("Installing WGCNA...")
    utils::install.packages("WGCNA")
    cell_cor <- WGCNA::cor(data, method = "pearson", nThreads = 1)
  }else{
    cell_cor <- WGCNA::cor(data, method = "pearson", nThreads = 1)
  }
  ## 4) TMM normalization factor
  if(verbose){
    message("Calculating TMM normalization factor...")
  }
  dge <- edgeR::DGEList(counts = data)
  error_detect <- try(dge <- edgeR::calcNormFactors(dge, method = "TMM"),
                      silent = TRUE)
  if("try-error" %in% class(error_detect)){
    TMM_factor <- NA
  }else{
    TMM_factor <- dge[["samples"]][["norm.factors"]]
  }
  ## 5) Effective library size
  if(verbose){
    message("Calculating effective library size...")
  }
  effective_library_size <- library_size * TMM_factor
  ## 6) Proportion of cell outliers
  if(verbose){
    message("Calculating proportion of cell outliers...")
  }
  q <- quantile(library_size)
  iqr <- stats::IQR(library_size)
  outlier_value_min <- q[2] - 1.5*iqr
  outlier_value_max <- q[4] + 1.5*iqr
  prop_outliers_cell <- sum(library_size < outlier_value_min | library_size > outlier_value_max)/ncol(data)
  ### list
  cell_properties <- dplyr::lst(library_size,
                                zero_fraction_cell,
                                cell_cor,
                                TMM_factor,
                                effective_library_size,
                                prop_outliers_cell)
  if(verbose){
    message("Done...")
  }
  return(cell_properties)
}



#' Calculate Gene Properties
#'
#' @param data A count matrix.
#' @param cpm_norm Whether CPM normalization will be used. Default is TRUE.
#' @param verbose Whether messages are returned during the process.
#' @importFrom Seurat CreateSeuratObject FindVariableFeatures HVFInfo
#' @return A list
#' @export
#' @examples
#' data <- matrix(rpois(100*100, 2),
#'                nrow = 100,
#'                dimnames = list(paste0('Gene', 1:100),
#'                                paste0('Cell', 1:100)))
#' result <- gene_properties(data, verbose = TRUE)
gene_properties <- function(data,
                            cpm_norm = TRUE,
                            verbose = FALSE){
  if(cpm_norm){
    if(verbose){
      message("Performing log2 CPM nomalization...")
    }
    norm_data <- log2(edgeR::cpm(data)+1)
  }
  ## 1) mean expression of log2 CPM of genes
  if(verbose){
    message("Calculating mean expression for genes...")
  }
  mean_expression <- apply(norm_data, 1, mean)
  ## 2) standard variance of log2 CPM of genes
  if(verbose){
    message("Calculating standard variance of genes...")
  }
  sd <- apply(norm_data, 1, sd)
  ## 3) coefficient of variance
  if(verbose){
    message("Calculating coefficient of variance...")
  }
  cv <- apply(data, 1, sd)/apply(data, 1, mean) * 100
  ## 4) fraction of zero in genes
  if(verbose){
    message("Calculating fraction of zero in genes...")
  }
  zero_fraction_gene <- apply(data, 1, function(x){
    sum(x == 0)/length(x)
  })
  ## 5) dispersion
  if(verbose){
    message("Calculating gene dispersions using Seurat...")
  }
  data_seurat <- Seurat::CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
  data_seurat <- Seurat::NormalizeData(data_seurat,
                                       normalization.method = "LogNormalize",
                                       scale.factor = 1e6)
  data_seurat <- Seurat::FindVariableFeatures(data_seurat, selection.method = "disp")
  dispersion <- Seurat::HVFInfo(data_seurat)$dispersion
  ## 6) Proportion of gene outliers
  if(verbose){
    message("Calculating proportion of gene outliers...")
  }
  q <- quantile(rowSums(data))
  iqr <- stats::IQR(rowSums(data))
  outlier_value_min <- q[2] - 1.5*iqr
  outlier_value_max <- q[4] + 1.5*iqr
  prop_outliers_gene <- sum(rowSums(data) < outlier_value_min | rowSums(data) > outlier_value_max)/nrow(data)
  ### list
  gene_properties <- dplyr::lst(mean_expression,
                                sd,
                                cv,
                                zero_fraction_gene,
                                dispersion,
                                prop_outliers_gene)
  if(verbose){
    message("Done...")
  }
  return(gene_properties)
}
