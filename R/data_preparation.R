#' Define Meta Information of A Dataset
#'
#' @param id Dataset ID
#' @param repository Which repository or cloud contains the dataset, e.g. GEO, ArrayExpress.
#' @param accession_number The accession number in the repository if available.
#' @param URL The URL to downloiad the dataset.
#' @param platform Sequencing platform.
#' @param species Organism.
#' @param organ Which source or organ do the samples come from?
#' @param cell_num The cell number of the dataset.
#' @param gene_num The gene number of the dataset.
#' @param data_type Which type is the dataset? e.g. count, FPKM.
#' @param ERCC Whether the dataset contains ERCC spike-in genes?
#' @param dilution_factor If there are spike-in genes, dilution factor is neccessary.
#' @param volume The volume (nanoliter) of mix liquid used in experiment.
#' @param group_condition Vector. The group assignment of cells in the dataset.
#' @param treatment Vector. The treatment of cells in the experiment.
#' @param batch_info Vector. The batch information of cells in the experiment.
#' @param cluster_info Vector. The cluster information of cells.
#' @param spatial_coordinate A data frame of x and y coordinates of spatial transcriptomics data
#' @param start_cell The cell id which is determined as the start_cell of a trajectory
#' @importFrom dplyr lst
#' @return A list of dataset information.
#' @export
meta_info <- function(
  id,
  repository = NULL,
  accession_number = NULL,
  URL = NULL,
  platform,
  species,
  organ = NULL,
  cell_num,
  gene_num,
  data_type = "count",
  ERCC = FALSE,
  dilution_factor = NULL,
  volume = NULL,
  group_condition = NULL,
  treatment = NULL,
  batch_info = NULL,
  cluster_info = NULL,
  spatial_coordinate = NULL,
  start_cell = NULL
){
  if(ERCC){
    if(is.null(dilution_factor) | is.null(volume)){
      stop("please input dilution factor and volume when ERCC spike-in genes are available.")
    }
  }
  dplyr::lst(id,
             repository,
             accession_number,
             URL,
             platform,
             species,
             organ,
             cell_num,
             gene_num,
             data_type,
             ERCC,
             dilution_factor,
             volume,
             group_condition,
             treatment,
             batch_info,
             cluster_info,
             spatial_coordinate,
             start_cell)
}


#' Define A Start Cell of A Trajectory
#'
#' @param meta_data The position information of cells or Spots
#'
#' @return A cell id
#' @export
#'
start_cell <- function(meta_data){
  if("true_y" %in% colnames(meta_data)){
    meta_data <- meta_data %>%
      mutate(y = true_y)
  }
  if("true_x" %in% colnames(meta_data)){
    meta_data <- meta_data %>%
      mutate(x = true_x)
  }
  meta_data <- meta_data %>%
    filter(label == "invasive cancer")
  median_cell <- c(mean(meta_data$x), mean(meta_data$y))
  distance <- sqrt((meta_data$x - median_cell[1]) ^ 2 + (meta_data$y - median_cell[2]) ^ 2)
  min_index <- which(distance == min(distance))
  paste0(meta_data$x[min_index], "x", meta_data$y[min_index])
}
