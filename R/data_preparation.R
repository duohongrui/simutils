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
  spatial_coordinate = NULL
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
             spatial_coordinate)
}
