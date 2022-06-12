#' Default Parameters of Simulation Methods
#'
#' This function is used to get the default parameters of a simulation method.
#'
#' @param method Method name
#' @importFrom splatter newSplatParams
#'
#' @export
#'
#' @examples
#' parameters <- default_parameters(method = "splat")
default_parameters <- function(method){
  if(method == "Splat"){
    parameters <- splatter::newSplatParams()
  }
  if(method == "Simple"){
    parameters <- splatter::newSimpleParams()
  }
  if(method == "Kersplat"){
    parameters <- splatter::newKersplatParams()
  }
  if(method == "SplatPop"){
    parameters <- splatter::newSplatPopParams()
  }
  if(method == "Lun"){
    parameters <- splatter::newLunParams()
  }
  if(method == "Lun2"){
    parameters <- splatter::newLun2Params()
  }

  return(parameters)
}
