#' @title Check Prior Information
#' @param method Which method
#' @param step Simulation or estimation
#' @param other_prior Prior information
#' @export
check_prior_info <- function(
    method,
    step,
    other_prior
){
  ## group.condition
  if(step == "estimation"){
    group_con_method <- c("Lun2",
                          "scDD",
                          "zingeR",
                          "SparseDC")
  }
  if(step == "simulation"){
    group_con_method <- c("zinbwaveZinger",
                          "zingeR")
  }

  if(method %in% group_con_method){
    if(is.null(other_prior[["group.condition"]])){
      stop(paste0("Please input group condition information when you use ", method))
    }else{
      other_prior_exec <- other_prior
    }
  }else{
    other_prior_exec <- other_prior
  }
  return(other_prior_exec)
}
