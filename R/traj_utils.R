#' Match Cells in Reference Data with Those in Simulation Data
#'
#' This function is used to match cells in reference data with those in simulation
#' data. Each cell in reference data corresponds to only one cell in
#' simulation data and the correlation between them is also calculated and returned.
#'
#' @param ref_data,sim_data A matrix or dynwrap object created by \code{\link[dynwrap]{wrap_expression}}.
#' Note that every row represents a cell and every column represents a gene.
#' Simulation and reference data must have the same size.
#' @return A list contains the correlation matrix and a dataframe of the paired
#' cells in reference and simulation data.
#' @export
#' @importFrom stats prcomp cor
#' @importFrom harmony HarmonyMatrix
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

  if(dynwrap::is_wrapper_with_expression(ref_data)){

    if(is.null(ref_data$expression)){

      ref_data <- ref_data$counts

      pca_ref <- prcomp(ref_data, rank. = 50)

    }else{

      ref_data <- ref_data$expression

      pca_ref <- prcomp(ref_data, rank. = 50)

    }
  }

  if(is.matrix(ref_data)){

    pca_ref <- prcomp(ref_data, rank. = 50)

  }

  if(dynwrap::is_wrapper_with_expression(sim_data)){

    if(is.null(sim_data$expression)){

      sim_data <- sim_data$counts

      pca_sim <- prcomp(sim_data, rank. = 50)

    }else{

      sim_data <- sim_data$expression

      pca_sim <- prcomp(sim_data, rank. = 50)

    }
  }

  if("SingleCellExperiment" %in% class(sim_data)){

    sim_data <- t(SingleCellExperiment::logcounts(sim_data))

    pca_sim <- prcomp(sim_data, rank. = 50)

  }

  if(is.list(sim_data)){

    pca_sim <- prcomp(sim_data$expression, rank. = 50)

  }

  if(is.matrix(sim_data)){

    pca_sim <- prcomp(sim_data, rank. = 50)

  }

  pca_input <- rbind(pca_ref$x, pca_sim$x)
  meta_data <- data.frame('Classification' = c(rep('reference', nrow(ref_data)),
                                               rep('simulation', nrow(ref_data))))
  set.seed(1)
  harmony_embeddings <- HarmonyMatrix(pca_input,
                                      meta_data,
                                      'Classification',
                                      do_pca=FALSE,
                                      verbose=FALSE)


  ### Calculate correlation matrix
  cor_result <- cor(t(harmony_embeddings), method = 'spearman')

  index <- dim(cor_result)[1]/2

  cor_result <- cor_result[(index+1):(index*2), 1:index]

  print(cor_result[1:5,1:5])

  a <- apply(cor_result, 2, function(x) return(unname(which(x==max(x)))))

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

  match_value <- apply(as.matrix(cell_pair), 1, value)

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
#' tree_format <- make_tree(ref_data, group = group)
make_tree <- function(ref_data, group, is_Newick=TRUE, is_parenthetic=FALSE){

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


.change_info <- function(model_ref, match_result){

  model_ref$cell_ids <- match_result[["cell_pair"]][["simulation"]]
  model_ref[["progressions"]][["cell_id"]] <- match_result[["cell_pair"]][["simulation"]]
  times <- dim(model_ref[["milestone_percentages"]])[1]/length(model_ref[["cell_ids"]])
  model_ref[["milestone_percentages"]][["cell_id"]] <- rep(match_result[["cell_pair"]][["simulation"]], times)

  return(model_ref)

}


#' Calculate the Correlation Between Geodesic Distances
#'
#' This function calculate the correlation between geodesic distances which refer
#' to the relative distance of one cell to all other cell in the trajectory. The
#' result is obtained from the mean values of choosing 0.05, 0.1, 0.15, 0.2,
#' 0.3, 0.35, 0.4, 0.45, 0.5 percent cells as waypoints.
#'
#'
#' @param model_ref,model_sim A matrix.
#' @param match_result The result generated by \code{\link[simutils]{match_cells}} function.
#' @return A value ranged from 0 to 1
#' @export
#' @importFrom dplyr arrange desc
#' @importFrom stats quantile
#'
#' @examples
#' # Check the docker status
#' # dynwrap::test_docker_installation(detailed = TRUE)
#'
#' # Open Terminal and execute the command
#' # docker pull dynverse/ti_slingshot:v1.0.3
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
#' # Trajectory inference

#'
#' #Generate a simulation data
#' set.seed(1)
#' b <- matrix(rpois(n = 2500, lambda = 1.5), nrow = 50)
#' rownames(b) <- paste0("fcell_", 1:ncol(b))
#' colnames(b) <- paste0("fgene_", 1:nrow(b))
#' dataset_sim <- dynwrap::wrap_expression(
#'   counts = b,
#'   expression = log2(b+1)
#' )
#' # Trajectory inference

#'
#' # Match cells
#' match_result <- match_cells(ref_data = dataset_ref,
#'                             sim_data = dataset_sim)

cal_cor_dist <- function(model_ref, model_sim, match_result){

  m <- arrange(match_result$cell_pair,
               desc(match_result$cell_pair$match_value))

  model_ref <- .change_info(model_ref = model_ref,
                            match_result = match_result)

  model_sim <- dynwrap::add_cell_waypoints(trajectory = model_sim,
                                           num_cells_selected = 50)

  model_ref <- dynwrap::add_cell_waypoints(trajectory = model_ref,
                                           num_cells_selected = 50)

  waypoints_num <- round(quantile(1:length(model_ref[["cell_ids"]]), seq(0.05, 0.5, 0.05)))

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


#' Calculate Four Metrics to Compare Two Trajectories
#'
#' @param model_ref,model_sim The object generated by \code{\link[dynwrap]{infer_trajectory}}.
#'
#' @return A list containing the results of four metrics.
#' @export
#'
#' @examples
#' # Check the docker status
#' # dynwrap::test_docker_installation(detailed = TRUE)
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

#' #Generate a simulation data
#' set.seed(1)
#' b <- matrix(rpois(n = 2500, lambda = 1.5), nrow = 50)
#' rownames(b) <- paste0("fcell_", 1:ncol(b))
#' colnames(b) <- paste0("fgene_", 1:nrow(b))
#' dataset_sim <- dynwrap::wrap_expression(
#'   counts = b,
#'   expression = log2(b+1)
#' )


calculate_traj_metrics <- function(model_ref,
                                   model_sim){
  ## Match cells
  match_result <- match_cells(ref_data = model_ref, sim_data = model_sim)

  ## Calculate correlation
  cor_dist <- cal_cor_dist(model_ref = model_ref,
                           model_sim = model_sim,
                           match_result = match_result)

  ## Change information in reference data
  model_ref <- .change_info(model_ref = model_ref,
                            match_result = match_result)

  ## Calculate metrics
  him <- dyneval::calculate_metrics(dataset = model_ref,
                                    model = model_sim,
                                    metrics = "him")

  f1_branches <- dyneval::calculate_metrics(dataset = model_ref,
                                            model = model_sim,
                                            metrics = "F1_branches")

  f1_milestones <- dyneval::calculate_metrics(dataset = model_ref,
                                              model = model_sim,
                                              metrics = "F1_milestones")

  return(list(him, f1_branches, f1_milestones, cor_dist))

}


#' Add Gene Expression Data to Model
#'
#' This function is used to add the gene expression data into the result of
#' trajectory inference generated by \code{\link[dynwrap]{infer_trajectory}} function.
#'
#' @param model A dynwrap::with_trajectory object generated by \code{\link[dynwrap]{infer_trajectory}}.
#' function.
#' @param dataset A dynwrap::data_wrapper object generated by \code{\link[dynwrap]{infer_trajectory}}.
#' function.
#'
#' @return A new dynwrap::with_trajectory object.
#' @export
#'
#' @examples
#' # Generate a reference data
#' set.seed(1)
#' a <- matrix(rpois(n = 2500, lambda = 2), nrow = 50)
#' rownames(a) <- paste0("cell", 1:ncol(a))
#' colnames(a) <- paste0("gene", 1:nrow(a))
#' dataset_ref <- dynwrap::wrap_expression(
#'   counts = a,
#'   expression = log2(a+1)
#' )
#' # Trajectory inference

add_data_to_model <- function(model, dataset){

  checkmate::assertClass(model, "dynwrap::with_trajectory")
  checkmate::assertClass(dataset, "dynwrap::data_wrapper")

  if(is.null(model[["counts"]])){

    model <- dynwrap::add_expression(dataset = model,
                                     counts = dataset[['counts']])

  }

  if(is.null(model[["expression"]])){

    model <- dynwrap::add_expression(dataset = model,
                                     expression = dataset[['expression']])

  }

  return(model)

}



#' Synthesize Fake Cells
#'
#' This function is used to synthesize fake cells via linear combination when the
#' number of cells is not the power of 2. First, we need cell group information
#' where real cells are sampled from. If no group information input, we perform
#' k-means algorithm on the data and use \code{\link[NbClust]{NbClust}} funciton
#' to determin the best cluster number. Finally, we merge the synthesized and
#' real data as the output result.
#'
#' @param dataset A matrix or the result generated by \code{\link[dynwrap]{wrap_expression}}
#' @param group A vector. Default is NULL.
#' @param seed Integer. A random seed.
#'
#' @return A list generated by \code{\link[dynwrap]{wrap_expression}}
#' @importFrom stats runif
#' @export
#'
#' @examples
#' set.seed(1)
#' a <- matrix(rpois(n = 2500, lambda = 2), nrow = 50)
#' rownames(a) <- paste0("cell_", 1:ncol(a))
#' colnames(a) <- paste0("gene_", 1:nrow(a))
#' dataset_ref <- dynwrap::wrap_expression(
#'   counts = a,
#'   expression = log2(a+1)
#' )
#' result <- syn_cell(dataset = a, seed = 2)
syn_cell <- function(dataset, group = NULL, seed){

  if(is.matrix(dataset)){

    dataset <- dynwrap::wrap_expression(counts = dataset,
                                        expression = log2(dataset+1))

  }

  if(!dynwrap::is_wrapper_with_grouping(dataset)){

    message("Performing k-means and determin the best number of clusters...")

    if(is.null(dataset[['expression']])){

      stop("No expression data is detected. Please use dynwrap::wrap_expression function\nto specify expression parameter.")

    }

    clust <- NbClust::NbClust(data = dataset[['expression']],
                              distance = 'euclidean',
                              min.nc = 2,
                              max.nc = sqrt(nrow(dataset[['expression']])),
                              method = "kmeans",
                              index = "silhouette")

    dataset <- dynwrap::add_grouping(dataset = dataset,
                                     grouping = clust[["Best.partition"]])

  }

  if(is.null(dataset[['counts']])){

    stop("No counts data is detected. Please use dynwrap::wrap_expression function\nto specify counts parameter.")

  }

  ncells <- length(dataset$cell_ids)

  if(log2(ncells) != as.integer(log2(ncells))){

    ## The number of additional cells which should be synthesized.
    diff_num <- 2^ceiling(log2(ncells))-ncells

    group <- dataset$group_ids

    group_num <- length(group)

    if(diff_num < group_num){

      group <- group[(group_num-diff_num+1):group_num]

    }


    ## Allocate cell number for every group
    num_allo <- c(rep(round(diff_num/group_num), group_num-1),
                  diff_num-sum(rep(round(diff_num/group_num), group_num-1)))

    add_syn_result <- matrix(ncol = dim(dataset[["counts"]])[2])


    ### Synthesize fake cells
    message("Synthesize fake cells...")
    for(i in 1:group_num){

      index <- which(dataset$grouping == group[i])

      set.seed(seed)

      add_syn <- matrix(data = runif(10*num_allo[i], min = 0, max = 0.2),
                        nrow = num_allo[i])


      set.seed(seed)

      tmp <- dataset[["counts"]][sample(index, 10, replace = TRUE), ]

      add_syn_matrix <- add_syn %*% tmp

      add_syn_matrix <- round(add_syn_matrix)

      add_syn_result <- rbind(add_syn_result, add_syn_matrix)

    }

    ## Add synthesized information in reference data
    add_syn_result <- add_syn_result[-1, ]

    message("Add the synthesized data to the real data...")
    syncell_id <- paste0(rep('syncell', diff_num), seq(1, diff_num))

    rownames(add_syn_result) <- syncell_id

    dataset[["cell_ids"]] <- c(dataset[["cell_ids"]], syncell_id)

    group_tmp <- rep(group, num_allo)

    names(group_tmp) <- syncell_id

    dataset[["grouping"]] <- c(dataset[["grouping"]], group_tmp)

    dataset[["counts"]] <- rbind(dataset[["counts"]], add_syn_result)

    dataset[["expression"]] <- rbind(dataset[["expression"]], log2(add_syn_result+1))

    message("Done")

  }

  return(dataset)

}


