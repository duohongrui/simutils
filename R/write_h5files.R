#' Write R object to .h5 Files
#'
#' This function facilitates the progress of converting any R object into .h5 file
#' to store information and read it in Python easily.
#' @param data An R object. Depend on our test, you would like to convert large
#' dataframes filled with numeric values into sparse matrix or matrix to save time.
#' @param file_path The file path of the h5file. Default is NULL and h5 file will
#' be generated automatically
#' @importFrom dynutils is_sparse
#' @importFrom hdf5r h5attr H5File
#' @importFrom methods as
#' @return An h5 file
#' @export
#'
#' @examples
#' # Generate a matrix
#' set.seed(1)
#' ref_data <- matrix(rpois(n = 10^4, lambda = 0.5), nrow = 1000)
#' colnames(ref_data) <- paste0("cell_", 1:ncol(ref_data))
#' rownames(ref_data) <- paste0("gene_", 1:nrow(ref_data))
#'
#' # Generate a sparse matrix
#' sparse_matrix <- Matrix::Matrix(data = ref_data, sparse = TRUE)
#'
#' # Generate a dataframe
#' df <- as.data.frame(sparse_matrix)
#'
#' # Generate a list
#' x <- list("length" = 2,
#'           "integer" = 2L,
#'           "logical" = TRUE,
#'           "sparse_matrix" = sparse_matrix,
#'           "df" = df,
#'           "matrix" = ref_data)
#'
#' # Set a path
#' file_path <- tempfile(fileext = ".h5")
#' # Write an h5 file
#' write_h5files(x, file_path = file_path)
#' # Load the h5 file
#' file.h5 <- hdf5r::H5File$new(file_path, mode = "r")
write_h5files <- function(data, file_path = NULL){
  if(is.null(file_path)){
    file_name <- tempfile(fileext = ".h5")
  }else{
    file_name <- file_path
  }
  cat(paste0("Your .h5 file is in: ", file_path, ".\n"))
  file.h5 <- hdf5r::H5File$new(file_name, mode = "w")
  on.exit(file.h5$close_all(), add = TRUE)
  .write_h5files(file.h5, data)
}


#' @param h5file A new h5 file
#' @param data The data remain to be write
#' @param name Default is NULL. One of the names of a list.
#' @rdname write_h5files
.write_h5files <- function(h5file, data, name = NULL){
  if(is.null(name)) subfile <- h5file else subfile <- h5file$create_group(name)
  if(is.null(data)){
    hdf5r::h5attr(subfile, "class") <- "null"
    subfile[[name]] <- "null"
  } else if(dynutils::is_sparse(data)){
    sparse_matrix <- methods::as(data, "dgCMatrix")
    hdf5r::h5attr(subfile, "class") <- "dgCMatrix"
    subfile[["i"]] <- as.vector(sparse_matrix@i)
    subfile[["p"]] <- sparse_matrix@p
    subfile[["Dim"]] <- sparse_matrix@Dim
    subfile[["x"]] <- sparse_matrix@x
    if(!is.null(sparse_matrix@factors)) subfile[["factors"]] <- "null"
    hdf5r::h5attr(subfile, "rownames") <- sparse_matrix@Dimnames[[1]]
    hdf5r::h5attr(subfile, "colnames") <- sparse_matrix@Dimnames[[2]]
  }else if(is.atomic(data)){
    if(is.matrix(data)){
      hdf5r::h5attr(subfile, "class") <- "matrix"
      if(!is.null(rownames(data))) hdf5r::h5attr(subfile, "rownames") <- rownames(data)
      if(!is.null(colnames(data))) hdf5r::h5attr(subfile, "colnames") <- colnames(data)
    }else if(is.character(data)){
      hdf5r::h5attr(subfile, "class") <- "character"
    }else if(is.integer(data)){
      hdf5r::h5attr(subfile, "class") <- "integer"
    }else if(is.numeric(data)){
      hdf5r::h5attr(subfile, "class") <- "numeric"
    }else if(is.null(data)){
      hdf5r::h5attr(subfile, "class") <- "null"; data <- "null"
    }else if(is.logical(data)){
      hdf5r::h5attr(subfile, "class") <- "logical"
      ifelse(data, "true", "false")
    }
    subfile[["data"]] <- data
  }else if(is.data.frame(data)){
    hdf5r::h5attr(subfile, "class") <- "data.frame"
    if(!is.null(names(data))) subfile[["names"]] <- names(data)
    if(!is.null(rownames(data))) hdf5r::h5attr(subfile, "rownames") <- rownames(data)
    if(!is.null(colnames(data))) hdf5r::h5attr(subfile, "colnames") <- colnames(data)
    subfile2 <- subfile$create_group("data")
    for(col_name in names(data)){
      .write_h5files(h5file = subfile2,
                     data = data[, col_name],
                     name = col_name)
    }
  }else if(is.list(data)){
    hdf5r::h5attr(subfile, "class") <- "list"
    if(!is.null(names(data))) subfile[["names"]] <- names(data)
    subfile2 <- subfile$create_group("data")
    for(name in names(data)){
      .write_h5files(h5file = subfile2,
                     data = data[[name]],
                     name = name)
    }
  }
}
