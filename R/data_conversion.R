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
      tmp_path <- tempdir()
      data_save_name <- file.path(tmp_path, paste0(time_string(), ".h5Seurat")) %>%
        simutils::fix_path()
    }
    SeuratDisk::SaveH5Seurat(simulate_result, filename = data_save_name)
    SeuratDisk::Convert(data_save_name, dest = "h5ad")

    # data path
    data_path <- stringr::str_replace(data_save_name,
                                      pattern = "h5Seurat",
                                      replacement = "h5ad")
    cat(glue::glue("Your data has been save to {data_path}", "\n"))
    simulate_result <- list(file_type = "h5ad",
                            save_path = data_path)
  }
  return(simulate_result)
}
