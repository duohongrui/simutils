#' Default Parameters of Simulation Methods
#'
#' This function is used to get the default parameters of a simulation method.
#'
#' @param method Method name
#' @import splatter
#'
#' @export
#'
#' @examples
#' parameters <- default_parameters(method = "Splat")
default_parameters <- function(method){
  if(method == "Splat" | method == "SCRIP"){
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
  if(method == "scDesign"){
    parameters <- NULL
  }
  if(method == "zinbwave"){
    parameters <- splatter::newZINBParams()
  }
  if(method == "BASiCS"){
    parameters <- splatter::newBASiCSParams()
  }
  if(method == "ESCO"){
    parameters <- ESCO::newescoParams()
  }
  if(method == "SPsimSeq"){
    parameters <- NULL
  }

  return(parameters)
}
