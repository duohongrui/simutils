# x <- simmethods::get_method(method = 'PROSSTT')
# x <- x["splat"]
#
# params <- get_default_value(x)
#
# get_default_value <- function(x){
#
#   method_name <- names(x)
#
#   param_list_name <- paste0(method_name, '_parameters')
#
#   default_value <- list()
#
#   for(i in seq_len(length(method_name))){
#
#     sublist <- x[[method_name[i]]][[param_list_name[i]]]
#
#     param_names <- names(sublist)
#
#     default_value[[method_name[i]]] <- purrr::map(param_names, function(x){
#
#       sublist[[x]][["default"]]
#
#     }) %>% stats::setNames(param_names)
#
#   }
#   return(default_value)
#
# }
#
#
#
# set.seed(1)
# ref_data <- matrix(rpois(n = 10^6, lambda = 0.5), nrow = 2000)
# colnames(ref_data) <- paste0("cell_", 1:ncol(ref_data))
# rownames(ref_data) <- paste0("gene_", 1:nrow(ref_data))
# sparse_matrix <- Matrix::Matrix(data = ref_data, sparse = T)
#
# df <- as.data.frame(sparse_matrix)
#
# file_name <- tempfile(fileext = ".h5")
# file.h5 <- hdf5r::H5File$new(file_name, mode = "w")
# # on.exit(file.h5$close_all(), add = TRUE)
#
#
# write_h5files <- function(h5file, data, name = NULL){
#
#   if(is.null(name)) subfile <- h5file else subfile <- h5file$create_group(name)
#
#   if(is.null(data)){
#     hdf5r::h5attr(subfile, "class") <- "null"
#     subfile[[name]] <- "null"
#   } else if(dynutils::is_sparse(data)){
#
#     sparse_matrix <- as(data, "dgCMatrix")
#
#     hdf5r::h5attr(subfile, "class") <- "dgCMatrix"
#
#     subfile[["i"]] <- as.vector(sparse_matrix@i)
#
#     subfile[["p"]] <- sparse_matrix@p
#
#     subfile[["Dim"]] <- sparse_matrix@Dim
#
#     subfile[["x"]] <- sparse_matrix@x
#
#     if(!is.null(sparse_matrix@factors)) subfile[["factors"]] <- "null"
#
#     hdf5r::h5attr(subfile, "rownames") <- sparse_matrix@Dimnames[[1]]
#
#     hdf5r::h5attr(subfile, "colnames") <- sparse_matrix@Dimnames[[2]]
#
#   }else if(is.atomic(data)){
#
#     hdf5r::h5attr(subfile, "class") <- "vector"
#
#     if(is.null(data)){
#       data <- "null"
#     }
#
#     if(is.logical(data)){
#       ifelse(data, "1", "0")
#     }
#
#     subfile[["data"]] <- data
#
#   }else if(is.data.frame(data)){
#
#     hdf5r::h5attr(subfile, "class") <- "data.frame"
#
#     if(!is.null(names(data))) subfile[["names"]] <- names(data)
#
#     if(!is.null(data)) hdf5r::h5attr(subfile, "rownames") <- rownames(data)
#
#     if(!is.null(data)) hdf5r::h5attr(subfile, "colnames") <- colnames(data)
#
#     for(col_name in names(data)){
#       write_h5files(h5file = subfile,
#                     data = data[, col_name],
#                     name = col_name)
#     }
#   }else if(is.matrix(data)){
#
#     hdf5r::h5attr(subfile, "class") <- "matrix"
#
#     subfile[["data"]] <- data
#
#   }else if(is.list(data)){
#
#     hdf5r::h5attr(subfile, "class") <- class(data)
#
#     if(!is.null(names(data))) subfile[["names"]] <- names(data)
#
#     subfile2 <- subfile$create_group("data")
#
#     for(name in names(data)){
#
#       write_h5files(h5file = subfile2,
#                     data = data[[name]],
#                     name = name)
#
#     }
#   }
# }
#
#
# x <- list("length" = 2,
#           "reference_datasets" = sparse_matrix,
#           "reference_df" = df,
#           "list" = params,
#           "reference_matrix" = ref_data)
#
# write_h5files(h5file = file.h5, data = x)
#
# a <- hdf5r::H5File$new(file_name, mode = "r")
#
#
# m <- a[["data"]][["reference_matrix"]][["data"]]
