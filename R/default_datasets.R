#' Default Datasets For Simulations
#'
#' This function is used to get the default datasets for simulations.
#'
#' @param method Method name
#' @import splatter
#' @importFrom scater mockSCE
#'
#' @export
#'
#' @examples
#' datasets <- default_datasets(method = "SCRIP")
default_datasets <- function(method){
  if(method == "SCRIP"){
    set.seed(1)
    datasets <- scater::mockSCE()
  }

  return(datasets)
}
