#' Trace the Depths of All Leaf Nodes
#'
#' @param edge The edge element in Phylo tree object
#'
#' @return A numeric vector of every leaf node
#' @export
#'
leaf_node_depth <- function(edge){
  ## root node
  root_node <- edge[1, 1]
  ## leaf nodes
  leaf_node <- 1:(root_node - 1)
  ## function for node trace
  node_trace <- function(node, depth){
    which_row <- which(edge$V2 == node)
    parent_node <- edge[which_row, "V1"]
    if(parent_node == root_node){
      return(depth)
    }else{
      depth <- depth + 1
      depth <- node_trace(parent_node, depth)
    }
  }
  ## blanck container
  depth_count <- c()
  ## processing bar
  if(!requireNamespace("progress")){
    install.packages("progress")
  }
  pb <- progress::progress_bar$new(
    format = "  Tracing :what [:bar] :percent eta: :eta",
    clear = FALSE, total = length(leaf_node), width = 100)
  ## trace each of the leaf nodes and search the depths of them
  for(i in leaf_node){
    ## process bar
    pb$tick(tokens = list(what = "every leaf nodes' depth "))
    ## default depth
    depth <- 1
    ## trace
    depth <- node_trace(i, depth)
    ## append
    depth_count <- append(depth_count, depth)
  }
  return(depth_count)
}



#' Calculate the Discrepency of Nodes' Height
#'
#' @param ref_tree Phylo tree derived from reference data
#' @param sim_tree Phylo tree derived from simulated data
#' @importFrom stats as.dendrogram
#'
#' @return An RMSE value
#' @export
#'
node_height_dis <- function(ref_tree, sim_tree){
  if(!requireNamespace("dendextend")){
    install.packages("dendextend")
  }
  if(!requireNamespace("MLmetrics")){
    install.packages("MLmetrics")
  }
  ref_dendro <- stats::as.dendrogram(ref_tree)
  ref_height <- dendextend::get_branches_heights(ref_dendro)
  sim_dendro <- stats::as.dendrogram(sim_tree)
  sim_height <- dendextend::get_branches_heights(sim_dendro)
  discrepency <- MLmetrics::RMSE(ref_height, sim_height)
  return(discrepency)
}

