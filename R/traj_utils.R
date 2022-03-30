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

  ### Match cells from simulation datasets with thoes in reference datasets
  ### according to correlation values.

  # cor_tmp <- cor_result
  #
  # match_result <- vector()
  #
  # for(i in 1:(nrow(cor_result)-1)){
  #
  #   print(max(cor_tmp))
  #
  #   print(i)
  #
  #   position <- which(cor_tmp==max(cor_tmp), arr.ind = TRUE)
  #
  #   if(dim(position)[1]!=1){
  #
  #     position <- position[1, ]
  #
  #   }
  #
  #   match_index <- rownames(cor_tmp)[position[1]]
  #
  #   names(match_index) <- colnames(cor_tmp)[position[2]]
  #
  #   match_result <- append(match_result, match_index)
  #
  #   if(dim(cor_tmp)[1]==2){
  #
  #     match_index_final <- base::setdiff(rownames(cor_tmp), match_index)
  #
  #     names(match_index_final) <- base::setdiff(colnames(cor_tmp), names(match_index))
  #
  #     match_result <- append(match_result, match_index_final)
  #
  #
  #   }else{
  #
  #     cor_tmp <- cor_tmp[-position[1], -position[2]]
  #
  #   }
  #
  # }




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
