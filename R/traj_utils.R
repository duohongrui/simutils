#' Match Cells in Reference Data with Those in Simulation Data
#'
#' This function is used to match cells in reference data with those in simulation
#' data. Each cell in reference data corresponds to only one cell in
#' simulation data and the correlation between them is also calculated and returned.
#'
#' @param ref_data,sim_data A matrix or dynwrap object created by `dynwrap::wrap_expression()`.
#' Note that every row represents a cell and every column represents a gene.
#' Simulation and reference data must have the same size.
#' @return A list contains the correlation matrix and a dataframe of the paired
#' cells in reference and simulation data.
#' @export
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
match_cells <- function(ref_data, sim_data){

  ### Remove batch effect using Harmony

  if("dynwrap::data_wrapper" %in% class(ref_data)){


    if(is.null(ref_data$expression)){

      ref_data <- ref_data$counts

      pca_ref <- stats::prcomp(ref_data, rank. = 50)

    }else{

      ref_data <- ref_data$expression

      pca_ref <- stats::prcomp(ref_data, rank. = 50)

    }
  }

  if("matrix" %in% class(ref_data)){

    pca_ref <- stats::prcomp(ref_data, rank. = 50)

  }

  if("dynwrap::data_wrapper" %in% class(sim_data)){


    if(is.null(sim_data$expression)){

      sim_data <- sim_data$counts

      pca_ref <- stats::prcomp(sim_data, rank. = 50)

    }else{

      sim_data <- sim_data$expression

      pca_ref <- stats::prcomp(sim_data, rank. = 50)

    }
  }

  if("SingleCellExperiment" %in% class(sim_data)){

    sim_data <- t(SingleCellExperiment::logcounts(sim_data))

    pca_sim <- stats::prcomp(sim_data, rank. = 50)

  }

  if("list" %in% class(sim_data)){

    pca_sim <- stats::prcomp(sim_data$expression, rank. = 50)

  }

  if("matrix" %in% class(sim_data)){

    pca_sim <- stats::prcomp(sim_data, rank. = 50)

  }

  pca_input <- rbind(pca_ref$x, pca_sim$x)

  meta_data <- data.frame('Classification' = c(rep('reference', nrow(ref_data)),
                                               rep('simulation', nrow(ref_data))))
  set.seed(1)
  harmony_embeddings <- harmony::HarmonyMatrix(pca_input,
                                               meta_data,
                                               'Classification',
                                               do_pca=FALSE,
                                               verbose=FALSE)


  ### Calculate correlation matrix
  cor_result <- stats::cor(t(harmony_embeddings), method = 'spearman')

  index <- dim(cor_result)[1]/2

  cor_result <- cor_result[(index+1):(index*2), 1:index]

  print(cor_result[1:5,1:5])

  a <- base::apply(cor_result, 2, function(x) return(unname(which(x==max(x)))))

  if(class(a)=='list'){

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


  match_value <- base::apply(as.matrix(cell_pair), 1, value)


  cell_pair <- cbind(cell_pair, match_value)

  return(list('cor_result'=cor_result,
              'cell_pair'=cell_pair))

}









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
#'
#' @examples
#' ref_data <- matrix(rpois(n = 2500, lambda = 2), nrow = 50)
#' rownames(ref_data) <- paste0("cell_", 1:ncol(ref_data))
#' colnames(ref_data) <- paste0("gene_", 1:nrow(ref_data))
#' group <- c(rep("A", 15), rep("B", 5), rep("C", 20), rep("D", 10))
#' tree_format <- make_tree(ref_data, group = group)
make_tree <- function(ref_data, group, is_Newick=TRUE, is_parenthetic=FALSE){

  data <- Seurat::CreateSeuratObject(counts = t(ref_data), verbose=FALSE)

  data <- Seurat::NormalizeData(data,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000,
                                verbose=FALSE)

  data <- Seurat::FindVariableFeatures(data,
                                       selection.method = "vst",
                                       nfeatures = 2000,
                                       verbose=FALSE)

  all.genes <- rownames(data)

  data <- Seurat::ScaleData(data, features = all.genes, verbose=FALSE)

  data@meta.data$'group' <- stringr::str_remove_all(group, pattern = "[(].+[)]")

  #Get state tree by hierachical clustering on the state means
  exp_data <- Seurat::AverageExpression(data, slot = 'scale.data', group.by = 'group')

  clu <- stats::hclust(stats::dist(t(exp_data$RNA)), method = 'ward.D')


  for(i in 1:length(clu[["labels"]])){

    clu[["labels"]][i] <- stringr::str_replace_all(clu[["labels"]][i], ",", "_")

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





.change_info <- function(model_ref, match_result){

  model_ref$cell_ids <- match_result[["cell_pair"]][["simulation"]]

  model_ref[["progressions"]][["cell_id"]] <- match_result[["cell_pair"]][["simulation"]]

  times <- dim(model_ref[["milestone_percentages"]])[1]/length(model_ref[["cell_ids"]])

  model_ref[["milestone_percentages"]][["cell_id"]] <- rep(match_result[["cell_pair"]][["simulation"]], times)

  return(model_ref)

}



#' Title Calculate the correlation between geodesic distances
#'
#' This function calculate the correlation between geodesic distances which refer
#' to the relative distance of one cell to all other cell in the trajectory. The
#' result is obtained from the mean values of choosing 0.05, 0.1, 0.15, 0.2,
#' 0.3, 0.35, 0.4, 0.45, 0.5 percent cells as waypoints.
#'
#'
#' @param model_ref,model_sim A matrix.
#' @param match_result The result generated by `match_cells` function.
#' @return A value ranged from 0 to 1
#' @export
#'
#' @examples
#' # Check the docker status
#' dynwrap::test_docker_installation(detailed = TRUE)
#'
#' # Generate a reference data
#' set.seed(1)
#' a <- matrix(rpois(n = 2500, lambda = 2), nrow = 50)
#' rownames(a) <- paste0("cell_", 1:ncol(a))
#' colnames(a) <- paste0("gene_", 1:nrow(a))

#' dataset_ref <- dynwrap::wrap_expression(
#'   counts = a,
#'   expression = log2(a+1)
#' )

#' model_ref <- dynwrap::infer_trajectory(dataset_ref, 'slingshot')

# Generate a simulation data
#' set.seed(1)
#' b <- matrix(rpois(n = 2500, lambda = 1.5), nrow = 50)
#' rownames(b) <- paste0("fcell_", 1:ncol(b))
#' colnames(b) <- paste0("fgene_", 1:nrow(b))

#' dataset_sim <- dynwrap::wrap_expression(
#'   counts = b,
#'   expression = log2(b+1)
#' )

#' model_sim <- dynwrap::infer_trajectory(dataset_sim, 'slingshot')

# Match cells
#' match_result <- simutils::match_cells(ref_data = dataset_ref, sim_data = dataset_sim)

#' cor_dist <- cal_corr(model_ref = model_ref, model_sim = model_sim, match_result = match_result)
cal_corr <- function(model_ref, model_sim, match_result){

  m <- dplyr::arrange(match_result$cell_pair,
                      dplyr::desc(match_result$cell_pair$match_value))
  print(utils::head(m))

  model_ref <- .change_info(model_ref = model_ref,
                            match_result = match_result)

  model_sim <- dynwrap::add_cell_waypoints(trajectory = model_sim,
                                           num_cells_selected = 50)

  model_ref <- dynwrap::add_cell_waypoints(trajectory = model_ref,
                                           num_cells_selected = 50)

  waypoints_num <- round(stats::quantile(1:length(model_ref[["cell_ids"]]), seq(0.05, 0.5, 0.05)))
  print(waypoints_num)

  record <- data.frame()

  for(cell_num in waypoints_num){

    waypoints_cell <- m$simulation[1:cell_num]

    model_sim$waypoint_cells <- waypoints_cell

    model_ref$waypoint_cells <- waypoints_cell

    Correlation <- dyneval::calculate_metrics(dataset = model_ref,
                                              model = model_sim,
                                              metrics = "correlation")

    print(Correlation)

    record <- rbind(record, Correlation)

  }

  return(mean(record$correlation))

}










