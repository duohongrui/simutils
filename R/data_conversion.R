#' Convert Data Format
#'
#' @param SCE_object A SingleCellExperiment object
#' @param return_format Which data format you want, list or Seurat or h5ad or
#' SingleCellExperiment
#' @param in_docker Logical. Do you process simulation step in Docker?
#' @param local_path If \code{in_docker} is TRUE, you must input the directory
#' path on your local device
#' @importFrom SeuratDisk SaveH5Seurat Convert
#' @importFrom glue glue
#' @return The object you want, list or Seurat or h5ad or SingleCellExperiment.
#' @export
#'
data_conversion <- function(
    SCE_object,
    return_format = "SingleCellExperiment",
    in_docker = FALSE,
    local_path = NULL
){
  simulate_result <- SCE_object
  if(return_format == "SingleCellExperiment"){
    simulate_result <- simulate_result
  }
  if(return_format == "list"){
    count_data <- SingleCellExperiment::counts(simulate_result)
    col_meta <- as.data.frame(SingleCellExperiment::colData(simulate_result))
    row_meta <- as.data.frame(SingleCellExperiment::rowData(simulate_result))
    simulate_result <- dplyr::lst(count_data,
                                  col_meta,
                                  row_meta)
  }
  if(return_format == "Seurat"){
    simulate_result <- Seurat::as.Seurat(simulate_result,
                                         counts = "counts",
                                         data = NULL)
  }
  if(return_format == "h5ad"){
    # Convert to Seurat object
    simulate_result <- Seurat::as.Seurat(simulate_result,
                                         counts = "counts",
                                         data = NULL)
    if(in_docker){
      if(is.null(local_path)){
        stop("Please input the directory path on your local device!")
      }
      # data path in Docker container
      tmp_path <- "/home/admin/docker_path"
      data_name <- paste0(simutils::time_string(), ".h5Seurat")
      data_save_name <- file.path(tmp_path, data_name) %>% simutils::fix_path()
    }else{
      # Tmp file
      if(is.null(local_path)){
        tmp_path <- tempdir()
      }else{
        tmp_path <- local_path
      }
      data_save_name <- file.path(tmp_path, paste0(simutils::time_string(), ".h5Seurat")) %>%
        simutils::fix_path()
    }
    SeuratDisk::SaveH5Seurat(simulate_result, filename = data_save_name)
    SeuratDisk::Convert(data_save_name, dest = "h5ad")

    # data path
    if(in_docker){
      data_save_name <- file.path(local_path, data_name)
    }
    data_path <- stringr::str_replace(data_save_name,
                                      pattern = "h5Seurat",
                                      replacement = "h5ad")
    cat(glue::glue("Your data has been save to {data_path}", "\n"))
    simulate_result <- list(file_type = "h5ad",
                            save_path = data_path)
  }
  return(simulate_result)
}



#' Convert R Matrix to Python h5ad File
#'
#' @param data A gene expression matrix with cell on columns and gene on rows.
#' @param data_id Tha data name to be output.
#' @param res The clustering resolution. Default is 0.15.
#' @param group The group information of each cells. If not NULL, the res will be none and the clustering step will not be performed.
#' @param save_to_path The save path on local device.
#' @param verbose If messages are returned during the process.
#' @export
scgan_data_conversion <- function(
    data,
    data_id,
    res,
    group = NULL,
    save_to_path,
    verbose
){
  # Process file path-----------------------------------------------------------
  ## 1. Create a temp directory to store .rds file or .h5 file for Python
  temp_dir <- tempdir()
  ## 2. Set a path for the input data
  temp_data_path <- file.path(temp_dir, "data.rds") %>% simutils::fix_path()

  # Process datasets------------------------------------------------------------
  ## 1. Convert ref_data into .rds or .h5 file and save to input path
  # simutils::write_h5files(data = ref_data, file_path = temp_input_path)
  saveRDS(data, file = temp_data_path)

  # Prepare the input parameters-----------------------------------------------
  ## 1. docker image working directory
  wd <- "/runcode"
  ## 2. local directory of the mount point
  local_path <- temp_dir %>% simutils::fix_path()
  ## 3. docker image directory of the mount point
  docker_path <- "/data"
  ## 4. verbose
  verbose <- verbose
  ## 5. args
  args <- c("/runcode/data_conversion.R")
  ## 6. command
  command <- "Rscript"
  ## 7. container id
  ### (1. Check docker installation
  docker_install <- dynwrap::test_docker_installation()
  if(!docker_install){
    stop("Docker has not been installed or started! Please check it!")
  }
  ### (2. Check the installation of simpipe docker image
  images <- babelwhale::list_docker_images() %>%
    tidyr::unite("Repository", "Tag", sep = ":", col = "Image") %>%
    dplyr::pull("Image")

  if(!"duohongrui/simutils_scgan:latest" %in% images){
    # If not, pull duohongrui/simpipe:latest
    babelwhale::pull_container(container_id = "duohongrui/simutils_scgan:latest")
  }
  ### (3. docker container id
  container_id <- "duohongrui/simutils_scgan"

  if(!is.null(group)){
    res <- NULL
  }

  # Save command parameters
  input_meta <- list(container_id = container_id,
                     command = command,
                     args = args,
                     volums = paste0(local_path, ":", docker_path),
                     workspace = wd,
                     verbose = verbose,
                     data_id = data_id,
                     res = res,
                     group = group)
  saveRDS(input_meta, file.path(local_path, "data_info.rds"))

  # Run container---------------------------------------------------------------
  output <- babelwhale::run(container_id = input_meta$container_id,
                            command = input_meta$command,
                            args = input_meta$args,
                            volumes = input_meta$volums,
                            workspace = input_meta$workspace,
                            verbose = input_meta$verbose,
                            debug = FALSE)

  # Get result------------------------------------------------------------------
  file.copy(from = file.path(local_path, paste0(data_id, ".h5ad")),
            to = file.path(save_to_path, paste0(data_id, ".h5ad")))
  file.copy(from = file.path(local_path, "cluster_number.rds"),
            to = file.path(save_to_path, "cluster_number.rds"))
  message("Output is saved to ", file.path(save_to_path, paste0(data_id, ".h5ad")))
  cluster_number <- readRDS(file.path(save_to_path, "cluster_number.rds"))
  return(list(save_path = file.path(save_to_path, paste0(data_id, ".h5ad")),
              cluster_number = cluster_number))
}
