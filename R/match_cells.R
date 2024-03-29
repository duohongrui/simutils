#' Match Cells in Reference Data with Those in Simulation Data
#'
#' This function is used to match cells in reference data with those in simulation
#' data. Each cell in reference data corresponds to only one cell in
#' simulation data and the correlation between them is also calculated and returned.
#'
#' @param ref_data,sim_data A matrix or dynwrap object created by \code{\link[dynwrap]{wrap_expression}}.
#' Note that every row represents a cell and every column represents a gene.
#' Simulation and reference data must have the same size.
#' @param t Whether the data matrix should be transposed when performing the PCA.
#' @param algorithm Optional. Which algorithm used for matching cells in simulated and real data. Improved_Hungarian (default), Hungarian.
#' @return A list contains the correlation matrix and a dataframe of the paired
#' cells in reference and simulation data.
#' @export
#' @importFrom stats prcomp
#' @importFrom utils head
#' @examples
#' # Generate a reference data
#' set.seed(1)
#' ref_data <- matrix(rpois(n = 2500, lambda = 2), nrow = 50)
#' rownames(ref_data) <- paste0("cell_", 1:ncol(ref_data))
#' colnames(ref_data) <- paste0("gene_", 1:nrow(ref_data))
#' # Generate a simulation data
#' set.seed(1)
#' sim_data <- matrix(rpois(n = 2500, lambda = 1.5), nrow = 50)
#' rownames(sim_data) <- paste0("fcell_", 1:ncol(sim_data))
#' colnames(sim_data) <- paste0("fgene_", 1:nrow(sim_data))
#' # Match cells
#' match_result <- match_cells(ref_data = ref_data, sim_data)
match_cells <- function(ref_data,
                        sim_data,
                        t = FALSE,
                        algorithm = "Improved_Hungarian")
{
  if(!requireNamespace("harmony", quietly = TRUE)){
    message("harmony is not installed on your device...")
    message("Installing harmony...")
    install.packages("harmony")
  }
  ### Remove batch effect using Harmony
  cat("Performing PCA...\n")
  if(dynwrap::is_data_wrapper(ref_data)){
    if(is.null(ref_data$expression)){
      ref_data <- ref_data$counts
      pca_ref <- prcomp(ref_data, rank. = 50)
    }else{
      ref_data <- ref_data$expression
      pca_ref <- prcomp(ref_data, rank. = 50)
    }
  }
  if(is.matrix(ref_data)){
    if(t){
      pca_ref <- prcomp(t(ref_data), rank. = 50)
    }else{
      pca_ref <- prcomp(ref_data, rank. = 50)
    }
  }
  if(dynwrap::is_data_wrapper(sim_data)){
    if(is.null(sim_data$expression)){
      sim_data <- sim_data$counts
      pca_sim <- prcomp(sim_data, rank. = 50)
    }else{
      sim_data <- sim_data$expression
      pca_sim <- prcomp(sim_data, rank. = 50)
    }
  }
  if(is.matrix(sim_data)){
    if(t){
      pca_sim <- prcomp(t(sim_data), rank. = 50)
    }else{
      pca_sim <- prcomp(sim_data, rank. = 50)
    }
  }
  pca_input <- rbind(pca_ref$x, pca_sim$x)
  meta_data <- data.frame('Classification' = c(rep('reference', nrow(ref_data)),
                                               rep('simulation', nrow(sim_data))))
  cat("Performing Harmony...\n")
  set.seed(1)
  harmony_embeddings <- harmony::HarmonyMatrix(pca_input,
                                               meta_data,
                                               'Classification',
                                               do_pca=FALSE,
                                               verbose=FALSE)
  ### Calculate correlation matrix
  cat("Calculate correlation matrix...\n")
  if(!requireNamespace("WGCNA", quietly = TRUE)){
    message("WGCNA is not installed on your device...")
    message("Installing WGCNA...")
    utils::install.packages("WGCNA")
    cor_result <- WGCNA::cor(t(harmony_embeddings), method = 'spearman')
  }else{
    cor_result <- WGCNA::cor(t(harmony_embeddings), method = 'spearman')
  }
  index <- dim(cor_result)[1]/2
  cor_result <- cor_result[(index+1):(index*2), 1:index]

  if(algorithm == "Improved_Hungarian"){
    cat("Match simulated and real cells using improved Hungarian...\n")
    a <- apply(cor_result, 2, function(x) return(unname(which(x==max(x)))))
    if(is.list(a)){
      a<- lapply(a, dplyr::first)
      a <- unlist(a)
    }
    record_ref_cell_rank <- data.frame('ref_cell' = colnames(cor_result),
                                       'rank' = rep(1, length(colnames(cor_result))))
    for(iteration in 1:length(a)^3){
      rep_cell <- as.numeric(names(which((table(a)>1) == T)))
      if(length(rep_cell) == 0){
        break
      }
      for(i in rep_cell){
        rep_ref_cell <- cor_result[, which(a == i)]
        max_ref_cell <- which(rep_ref_cell[i, ] == max(rep_ref_cell[i, ]))
        if(length(max_ref_cell) >1){
          max_ref_cell <- max_ref_cell[1]
        }
        second_rank_ref_cell <- which(a == i)[-max_ref_cell]
        rank_change_index <- which(record_ref_cell_rank[, 1] %in% names(second_rank_ref_cell))
        record_ref_cell_rank[rank_change_index, 2] <- record_ref_cell_rank[rank_change_index, 2]+1
        rep_ref_cell <- as.data.frame(rep_ref_cell[, -max_ref_cell])
        if(dim(rep_ref_cell)[2] == 1){
          colnames(rep_ref_cell) <- names(second_rank_ref_cell)
        }
        ncell <- length(colnames(rep_ref_cell))
        for(w in 1:ncell){
          name <- colnames(rep_ref_cell)[w]
          b <- sort(rep_ref_cell[, w], decreasing = T)
          names(b)<- c(1:length(b))
          allo_index <- record_ref_cell_rank[which(record_ref_cell_rank$ref_cell == name),2]
          allo_num <- b[which(names(b) == allo_index)]
          if(length(b[which(b==allo_num)]) >1){
            same_index <- which(names(b[which(b == allo_num)]) == allo_index)
            allo_index <- which(rep_ref_cell[, w] == allo_num)
            allo_index <- allo_index[same_index]
          }else{allo_index <- which(rep_ref_cell[, w] == allo_num)}
          names(allo_index) <- name
          a[which(names(a) == name)] <- allo_index
        }
      }
    }
    cell_pair <- data.frame('reference'=names(a),
                            'simulation'=rownames(cor_result)[a])
    value <- function(x){
      index <- c(which(colnames(cor_result) == x[1]),
                 which(rownames(cor_result) == x[2]))
      value_tmp <- cor_result[index[2], index[1]]
    }
    match_value <- apply(cell_pair, 1, value)
    cell_pair <- cbind(cell_pair, match_value)
    print(utils::head(cell_pair))
    cost <- sum(-match_value)
  }
  if(algorithm == "Hungarian"){
    if(!requireNamespace("RcppHungarian", quietly = TRUE)){
      message("RcppHungarian is not installed on your device...")
      message("Installing RcppHungarian...")
      install.packages("RcppHungarian")
    }
    cat("Match simulated and real cells using Hungarian...\n")
    Hungarian_result <- RcppHungarian::HungarianSolver(-cor_result)
    match_value <- purrr::map(1:ncol(cor_result), .f = function(x){
      cor_result[Hungarian_result[["pairs"]][x, 1], Hungarian_result[["pairs"]][x, 2]]
    }) %>% unlist()
    cell_pair <- data.frame("reference" = colnames(cor_result)[Hungarian_result[["pairs"]][, 2]],
                            "simulation" = rownames(cor_result)[Hungarian_result[["pairs"]][, 1]],
                            "match_value" = match_value)
    print(utils::head(cell_pair))
    cost <- Hungarian_result[["cost"]]
  }
  return(list("PCA_raw" = pca_input,
              "cor_result" = cor_result,
              "harmony_embeddings" = harmony_embeddings,
              "cell_pair" = cell_pair,
              "cost" = cost))
}
