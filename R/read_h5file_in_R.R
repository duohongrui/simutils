#' Read An H5 File into R
#'
#' @param file_path The file path of the candidate .h5 file
#' @importFrom purrr map_dfc map
#' @importFrom Matrix sparseMatrix
#' @return An R object
#' @export
#'
#' @examples
#' # Generate a list
#' x <- list("length" = 2,
#'           "integer" = 2L,
#'           "logical" = TRUE)
#'
#' # Set a path
#' file_path <- tempfile(fileext = ".h5")
#' # Write an h5 file
#' write_h5files(x, file_path = file_path)
#' # Read the h5 file into R
#' r_object <- h5_to_r(file_path)
h5_to_r <- function(file_path){
  h5file <- hdf5r::H5File$new(file_path, mode = "r")
  .h5_to_r(h5file)
}


#' @param h5file An h5 file
#' @rdname h5_to_r
.h5_to_r <- function(h5file){
  class_check <- hdf5r::h5attr(h5file, "class")
  if(class_check == "data.frame"){
    col_names <- hdf5r::h5attr(h5file, "colnames")
    row_names <- hdf5r::h5attr(h5file, "rownames")
    subh5file <- h5file[["data"]]
    r_output <- base::suppressMessages(as.data.frame(purrr::map_dfc(col_names, ~subh5file[[.]][["data"]][])))
    if(!is.null(row_names)) rownames(r_output) <- row_names
    if(!is.null(col_names)) colnames(r_output) <- col_names
    r_output
  }else if(class_check == "dgCMatrix"){
    col_names <- if(!is.null(hdf5r::h5attr(h5file, "colnames"))) hdf5r::h5attr(h5file, "colnames") else NULL
    row_names <- if(!is.null(hdf5r::h5attr(h5file, "rownames"))) hdf5r::h5attr(h5file, "rownames") else NULL
    Matrix::sparseMatrix(i = h5file[["i"]][],
                         p = h5file[["p"]][],
                         x = h5file[["x"]][],
                         dims = h5file[["Dim"]][],
                         dimnames = list(row_names, col_names),
                         index1 = FALSE)
  }else if(class_check == "character" | class_check == "numeric"){
    h5file[["data"]][]
  }else if(class_check == "matrix"){
    as.integer(h5file[["data"]][,])
  }else if(class_check == "integer"){
    as.integer(h5file[["data"]][])
  }else if(class_check == "null"){
    NULL
  }else if(class_check == "logical"){
    ifelse(h5file[["data"]][] == "true", TRUE, FALSE)
  }else if(class_check == "list"){
    if("names" %in% names(h5file)) object_names <- h5file[["names"]][]
    subh5file <- h5file[["data"]]
    r_output <- purrr::map(object_names, ~.h5_to_r(subh5file[[.]])) %>% stats::setNames(object_names)
    r_output
  }
}
