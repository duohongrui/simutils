x <- simmethods::get_method()

params <- get_default_value(x)

get_default_value <- function(x){

  method_name <- names(x)

  param_list_name <- paste0(method_name, '_parameters')

  default_value <- list()

  for(i in seq_len(length(method_name))){

    sublist <- x[[method_name[i]]][[param_list_name[i]]]

    param_names <- names(sublist)

    default_value[[method_name[i]]] <- purrr::map(param_names, function(x){

      sublist[[x]][["default"]]

    }) %>% stats::setNames(param_names)

  }
  return(default_value)

}


file_name <- tempfile(fileext = ".h5")
file.h5 <- hdf5r::H5File$new(file_name, mode = "w")
# on.exit(file.h5$close_all(), add = TRUE)
x <- params


write_h5files <- function(h5file, data, name = NULL){

  if(is.null(name)) subfile <- h5file else subfile <- h5file$create_group(name)

  if(is.null(data)){
    hdf5r::h5attr(subfile, "class") <- "NULL"
    subfile[[name]] <- "NULL"
  }

  if(is.list(data)){

    hdf5r::h5attr(subfile, "class") <- class(data)

    subfile[["names"]] <- names(data)

    subfile2 <- subfile$create_group("data")

    for(name in names(data)){

      write_h5files(h5file = subfile2,
                    data = data[[name]],
                    name = name)

    }
  }

  if(is.vector(data)){
    hdf5r::h5attr(subfile, "class") <- class(data)
    if(is.logical(data)){
      ifelse(data, "true", "false")
    }
    subfile[[name]] <- data
  }

  if(class(data) == "SplatParams"){
    hdf5r::h5attr(subfile, "class") <- "SplatParams"
    subfile[[name]] <- data
  }

}


write_h5files(h5file = file.h5, data = x)

a <- hdf5r::H5File$new(file_name, mode = "r")



