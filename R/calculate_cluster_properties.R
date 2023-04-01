#' Calculate Clustering Deviation Index (CDI)
#'
#' @param data A matrix of gene expression profile.
#' @param cluster_info Cluster(or group) assignment of every cells in columns of matrix.
#' @param batch_info Default is NULL. Batch information of every cells.
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
  if(!requireNamespace("CDI", quietly = TRUE)){
    message("Installing CDI package...")
    remotes::install_github('jichunxie/CDI')
  }
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
  CDI_result <- CDI::calculate_CDI(
    sub_gcmat = sub_data,
    cand_lab_df = cell_label,
    batch_label = batch_info,
    cell_size_factor = size_factor
  )
  return(CDI_result)
}



#' Calculate ROUGE Value
#' Calculate entropy-based metric for assessing the purity of single cell populations.
#'
#' @param data A matrix of gene expression profile.
#' @param cluster_info Cluster(or group) assignment of every cells in columns of matrix.
#' @return A numeric ranged from 0 to 1.
#' @export
#' @references
#' Liu B, Li C, Li Z, et al. An entropy-based metric for assessing the purity of single cell populations. Nature communications, 2020, 11(1): 1-13. <https://doi.org/10.1038/s41467-020-16904-3>
#' Github URL: <https://github.com/PaulingLiu/ROGUE>
calculate_ROGUE <- function(
  data,
  cluster_info
){
  if(!requireNamespace("ROGUE", quietly = TRUE)){
    message("Installing ROGUE package...")
    devtools::install_github('PaulingLiu/ROGUE')
  }
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



#' Calculate Average Silouette Width
#'
#' @param data A matrix of gene expression profile or the distance matrix.
#' @param cluster_info Cluster(or group) assignment of every cells in columns of matrix.
#' @importFrom parallelDist parDist
#' @importFrom utils install.packages
#' @export
calculate_silhouette <- function(
  data,
  cluster_info
){
  if(!requireNamespace("cluster", quietly = TRUE)){
    message("Installing cluster package...")
    utils::install.packages("cluster")
  }
  if(!requireNamespace("parallelDist", quietly = TRUE)){
    message("Installing parallelDist package...")
    utils::install.packages("parallelDist")
  }
  if(class(data) == "dist"){
    dist <- data
  }else{
    dist <- parallelDist::parDist(t(data))
  }
  silhouette_width <- cluster::silhouette(x = as.numeric(as.factor(cluster_info)),
                                          dist)
  average_silhouette <- mean(silhouette_width[, 3])
  return(average_silhouette)
}



#' Calculate Dunn Index
#'
#' @param data A matrix of gene expression profile or the distance matrix.
#' @param cluster_info Cluster(or group) assignment of every cells in columns of matrix.
#'
#' @export
calculate_dunn <- function(
  data,
  cluster_info
){
  if(!requireNamespace("clValid", quietly = TRUE)){
    message("Installing clValid package...")
    utils::install.packages("clValid")
  }
  if(!requireNamespace("parallelDist", quietly = TRUE)){
    message("Installing parallelDist package...")
    utils::install.packages("parallelDist")
  }
  if(class(data) == "dist"){
    dist <- data
  }else{
    dist <- parallelDist::parDist(t(data))
  }
  cluster_info <- as.numeric(as.factor(cluster_info))
  dunn <- clValid::dunn(distance = dist, clusters = cluster_info)
  return(dunn)
}



#' Calculate Connectivity
#'
#' @param data A matrix of gene expression profile or the distance matrix.
#' @param cluster_info Cluster(or group) assignment of every cells in columns of matrix.
#'
#' @export
calculate_connectivity <- function(
  data,
  cluster_info
){
  if(!requireNamespace("clValid", quietly = TRUE)){
    message("Installing clValid package...")
    utils::install.packages("clValid")
  }
  if(!requireNamespace("parallelDist", quietly = TRUE)){
    message("Installing parallelDist package...")
    utils::install.packages("parallelDist")
  }
  if(class(data) == "dist"){
    dist <- data
  }else{
    dist <- parallelDist::parDist(t(data))
  }
  cluster_info <- as.numeric(as.factor(cluster_info))
  con <- clValid::connectivity(distance = dist, clusters = cluster_info)
  return(con)
}



#' Calculate Davies-Bouldin Index
#'
#' @param data A matrix of gene expression profile.
#' @param cluster_info Cluster(or group) assignment of every cells in columns of matrix.
#'
#' @export
#'
calculate_DB_index <- function(
  data,
  cluster_info
){
  if(!requireNamespace("clusterSim", quietly = TRUE)){
    message("Installing clusterSim package...")
    utils::install.packages("clusterSim")
  }
  cluster_info <- as.numeric(as.factor(cluster_info))
  DB <- clusterSim::index.DB(t(data), cluster_info)$DB
  return(DB)
}



#' Calculate Calinski-Harabasz Index
#'
#' @param data A matrix of gene expression profile.
#' @param cluster_info Cluster(or group) assignment of every cells in columns of matrix.
#'
#' @export
#'
calculate_CH_index <- function(
    data,
    cluster_info
){
  if(!requireNamespace("fpc", quietly = TRUE)){
    message("Installing fpc package...")
    utils::install.packages("fpc")
  }
  cluster_info <- as.numeric(as.factor(cluster_info))
  CH <- fpc::calinhara(t(data), cluster_info)
  return(CH)
}



#' Summarize the Ability of Simulating Clusters
#'
#' @param data A matrix or dataframe of gene expression
#' @param dist Optionally, if NULL, the distance matrix is computed
#' @param cluster_info The vector of characters which each cell belongs to
#' @param threads How many cores used for parallel computation
#' @param verbose Whether the messages of execution process are returned
#'
#' @return A list of metric results
#' @export
#'
calculate_cluster_properties <- function(
  data,
  dist = NULL,
  cluster_info,
  threads = 1,
  verbose = TRUE
){
  if(is.data.frame(data)){
    data <- as.matrix(data)
  }
  if(is.null(dist)){
    if(!requireNamespace("parallelDist", quietly = TRUE)){
      install.packages("parallelDist")
    }
    dist <- parallelDist::parDist(t(data), threads = threads)
  }

  if(verbose){
    message("1-Calculating CDI...")
  }
  error <- try(
    CDI <- simutils::calculate_CDI(data, cluster_info = cluster_info),
    silent = TRUE
  )
  if("try-error" %in% class(error)){
    warning("The CDI calculation failed")
    CDI <- NA
  }else{
    CDI <- min(CDI[1, 1], CDI[1, 2])
  }

  if(verbose){
    message("2-Calculating ROUGE...")
  }
  error <- try(
    ROUGE <- simutils::calculate_ROGUE(data, cluster_info = cluster_info),
    silent = TRUE
  )
  if("try-error" %in% class(error)){
    warning("The ROUGE calculation failed")
    ROUGE <- NA
  }

  if(verbose){
    message("3-Calculating silhouette...")
  }
  silhouette <- simutils::calculate_silhouette(dist, cluster_info = cluster_info)

  if(verbose){
    message("4-Calculating dunn...")
  }
  dunn <- simutils::calculate_dunn(dist, cluster_info = cluster_info)

  if(verbose){
    message("5-Calculating connectivity...")
  }
  connectivity <- simutils::calculate_connectivity(dist, cluster_info = cluster_info)

  if(verbose){
    message("6-Calculating DB index...")
  }
  DB_index <- simutils::calculate_DB_index(data, cluster_info = cluster_info)

  group_metrics <- dplyr::lst(CDI,
                              ROUGE,
                              silhouette,
                              dunn,
                              connectivity,
                              DB_index)
  return(group_metrics)
}



