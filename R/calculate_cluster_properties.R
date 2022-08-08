#' Calculate Clustering Deviation Index (CDI)
#'
#' @param data A matrix of gene expression profile.
#' @param cluster_info Cluster(or group) assignment of every cells in columns of matrix.
#' @param batch_info Default is NULL. Batch information of every cells.
#' @importFrom CDI feature_gene_selection size_factor calculate_CDI
#' @return A dataframe of CDI values
#' @export
#' @references
#' Fang J, Chan C, Owzar K, et al. Clustering Deviation Index (CDI): A robust and accurate unsupervised measure for evaluating scRNA-seq data clustering. bioRxiv, 2022. <https://doi.org/10.1101/2022.01.03.474840>
#' Github URL: <https://github.com/jichunxie/CDI>
calculate_CDI <- function(
  data,
  cluster_info,
  batch_info = NULL
){
  ## cluster
  cell_label <- data.frame("clsuter" = as.numeric(as.factor(cluster_info)))
  ## feature selection
  feature_gene_index <- CDI::feature_gene_selection(
    gcmat = data,
    batch_label = batch_info,
    method = "wds",
    nfeature = 1000
  )
  ## subset
  sub_data <- data[feature_gene_index, ]
  ## size factor
  size_factor <- CDI::size_factor(data)
  ## calculate CDI
  t1 <- Sys.time()
  CDI_result <- CDI::calculate_CDI(
    sub_gcmat = sub_data,
    cand_lab_df = cell_label,
    batch_label = batch_info,
    cell_size_factor = size_factor
  )
  t2 <- Sys.time()
  print(difftime(t2, t1))
  return(CDI_result)
}



#' Calculate ROUGE Value
#' Calculate entropy-based metric for assessing the purity of single cell populations.
#'
#' @param data A matrix of gene expression profile.
#' @param cluster_info Cluster(or group) assignment of every cells in columns of matrix.
#' @importFrom ROGUE matr.filter SE_fun CalculateRogue
#' @return A numeric ranged from 0 to 1.
#' @export
#' @references
#' Liu B, Li C, Li Z, et al. An entropy-based metric for assessing the purity of single cell populations. Nature communications, 2020, 11(1): 1-13. <https://doi.org/10.1038/s41467-020-16904-3>
#' Github URL: <https://github.com/PaulingLiu/ROGUE>
calculate_ROGUE <- function(
  data,
  cluster_info
){
  ## cluster name
  cluster_name <- unique(cluster_info)
  rogue <- c()
  ## cluster proportions
  proportion_cluster <- table(cluster_info)/ncol(data)
  ## calculate ROGUE value for every cluster
  for(i in seq_len(length(cluster_name))){
    name <- cluster_name[i]
    index <- which(cluster_info == name)
    sub_data <- data[, index]
    expr <- ROGUE::matr.filter(sub_data, min.cells = 10, min.genes = 10)
    ent_res <- ROGUE::SE_fun(expr)
    rogue_value <- ROGUE::CalculateRogue(ent_res, platform = "UMI")
    rogue <- append(rogue, rogue_value)
  }
  ## summary
  rogue_result <- data.frame("cluster/group" = cluster_name,
                             "prop_cluster" = proportion_cluster,
                             "rogue_value" = rogue)
  final_rogue <- sum(rogue_result$prop_cluster.Freq * rogue_result$rogue_value)
  return(final_rogue)
}

