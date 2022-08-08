#' Calculate Cell Properties
#'
#' @param data A count matrix.
#' @param verbose Whether messages are returned during the process.
#' @importFrom edgeR cpm DGEList calcNormFactors
#' @importFrom WGCNA cor
#' @importFrom dplyr lst
#' @return A list
#' @export
#' @examples
#' data <- matrix(rpois(100*100, 2),
#'                nrow = 100,
#'                dimnames = list(paste0('Gene', 1:100),
#'                                paste0('Cell', 1:100)))
#' result <- cell_properties(data, verbose = TRUE)
#' str(result)
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
  cell_cor <- WGCNA::cor(data, method = "pearson", nThreads = 1)
  ## 4) TMM normalization factor
  if(verbose){
    message("Calculating TMM normalization factor...")
  }
  dge <- edgeR::DGEList(counts = data)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  TMM_factor <- dge[["samples"]][["norm.factors"]]
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
  iqr <- IQR(library_size)
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
#' str(result)
gene_properties <- function(data,
                            cpm_norm = TRUE,
                            verbose = FALSE){
  if(cpm_norm){
    if(verbose){
      message("Performing log2 CPM nomalization...")
    }
    data <- log2(edgeR::cpm(data)+1)
  }
  ## 1) mean expression of log2 CPM of genes
  if(verbose){
    message("Calculating mean expression for genes...")
  }
  mean_expression <- apply(data, 1, mean)
  ## 2) standard variance of genes
  if(verbose){
    message("Calculating standard variance of genes...")
  }
  sd <- apply(data, 1, sd)
  ## 3) coefficient of variance
  if(verbose){
    message("Calculating coefficient of variance...")
  }
  cv <- sd/mean_expression * 100
  ## 4) gene correlation
  if(verbose){
    message("Calculating gene correlation...")
  }
  gene_cor <- WGCNA::cor(t(data), method = "spearman", nThreads = 1)
  ## 5) fraction of zero in genes
  if(verbose){
    message("Calculating fraction of zero in genes...")
  }
  zero_fraction_gene <- apply(data, 1, function(x){
    sum(x == 0)/length(x)
  })
  ## 6) dispersion
  if(verbose){
    message("Calculating gene dispersions using Seurat...")
  }
  data_seurat <- Seurat::CreateSeuratObject(counts = data)
  data_seurat <- Seurat::FindVariableFeatures(data_seurat, selection.method = "disp")
  dispersion <- Seurat::HVFInfo(data_seurat)$dispersion
  ## 7) Proportion of gene outliers
  if(verbose){
    message("Calculating proportion of gene outliers...")
  }
  q <- quantile(rowSums(data))
  iqr <- IQR(rowSums(data))
  outlier_value_min <- q[2] - 1.5*iqr
  outlier_value_max <- q[4] + 1.5*iqr
  prop_outliers_gene <- sum(rowSums(data) < outlier_value_min | rowSums(data) > outlier_value_max)/nrow(data)
  ### list
  gene_properties <- dplyr::lst(mean_expression,
                                sd,
                                cv,
                                gene_cor,
                                zero_fraction_gene,
                                dispersion,
                                prop_outliers_gene)
  if(verbose){
    message("Done...")
  }
  return(gene_properties)
}