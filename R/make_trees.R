#' Convert An Expression Matrix into Tree Format
#'
#' This function is used to convert the gene expression matrix into the tree
#' format.
#'
#'
#' @param ref_data A Matrix. The gene expression matrix where every row represents
#' a cell and every column represents a gene.
#' @param group A Vector. The group information of the cells in ref_data
#' @param is_Newick If TRUE, return the Newick format.
#' @param is_parenthetic If TRUE, return the parenthetic format.
#' @return A file of parenthetic format or Newick format
#' @export
#' @importFrom  Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#' ScaleData AverageExpression
#' @importFrom stringr str_remove_all str_replace_all
#' @importFrom stats hclust dist
#'
#' @examples
#' ref_data <- matrix(rpois(n = 2500, lambda = 2), nrow = 50)
#' rownames(ref_data) <- paste0("cell_", 1:ncol(ref_data))
#' colnames(ref_data) <- paste0("gene_", 1:nrow(ref_data))
#' group <- c(rep("A", 15), rep("B", 5), rep("C", 20), rep("D", 10))
#' tree_format <- make_trees(ref_data, group = group)
make_trees <- function(ref_data, group, is_Newick=TRUE, is_parenthetic=FALSE){

  data <- CreateSeuratObject(counts = t(ref_data), verbose=FALSE)

  data <- NormalizeData(data,
                        normalization.method = "LogNormalize",
                        scale.factor = 10000,
                        verbose=FALSE)

  data <- FindVariableFeatures(data,
                               selection.method = "vst",
                               nfeatures = 2000,
                               verbose=FALSE)
  all.genes <- rownames(data)

  data <- ScaleData(data, features = all.genes, verbose=FALSE)

  data@meta.data$'group' <- str_remove_all(group, pattern = "[(].+[)]")

  #Get state tree by hierachical clustering on the state means
  exp_data <- AverageExpression(data, slot = 'scale.data', group.by = 'group')

  clu <- hclust(dist(t(exp_data$RNA)), method = 'ward.D')

  for(i in 1:length(clu[["labels"]])){

    clu[["labels"]][i] <- str_replace_all(clu[["labels"]][i], ",", "_")

  }

  phyla <- ctc::hc2Newick(clu, flat=TRUE)

  if(is_Newick){

    return(phyla)

  }

  phyla <- ape::read.tree(text = phyla)

  if(is_parenthetic){

    return(list(phyla))

  }

  phyla$edge.length <- ceiling(phyla$edge.length)

  phyla$edge.length[phyla$edge.length == 0] <- 1

  return(phyla)

}