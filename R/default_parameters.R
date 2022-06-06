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
  if(method == "splat"){
    parameters <- splatter::newSplatParams()
  }

  return(parameters)
}
