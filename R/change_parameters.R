#' Change parameters
#'
#' @param function_expr The function expression used for get the parameters.
#' @param other_prior The list of new parameters remained to replace.
#' @param step Which step you are in.
#'
#' @export
change_parameters <- function(
    function_expr,
    other_prior,
    step
){
  if(step == "estimation"){
    estimate_formals <- as.list(formals(eval(parse(text = function_expr))))
    for(param in names(estimate_formals)){
      names_wait_check <- names(other_prior)
      if(param %in% names_wait_check){
        estimate_formals[[param]] <- other_prior[[param]]
      }
    }
    return(estimate_formals)
  }

  if(step == "simulation"){
    simulate_formals <- as.list(formals(eval(parse(text = function_expr))))
    for(param in names(simulate_formals)){
      names_wait_check <- names(other_prior)
      if(param %in% names_wait_check){
        simulate_formals[[param]] <- other_prior[[param]]
      }
    }
    return(simulate_formals)
  }
}



#' Change Parameters About scGAN
#'
#' @param project_name Name your experiment
#' @param ... Other parameters
#' @importFrom rjson fromJSON
#'
#' @export
#'
change_scGAN_parameters <- function(
    project_name,
    ...
){
  scGAN_params <- system.file("scGAN_parameters.json", package = "simutils")
  params <- rjson::fromJSON(file = scGAN_params)
  names(params[["experiments"]]) <- project_name
  params <- change_values_in_list(list = params, ...)
  return(params)
}

#' Change Values From a List
#'
#' @param list The list
#' @param ... Other attributes and new values
#' @importFrom rrapply rrapply
#'
#' @export
#'
change_values_in_list <- function(
    list,
    ...
){
  reset_params <- list(...)
  for(id in names(reset_params)){
    reset_value <- reset_params[[id]]
    list <- rrapply::rrapply(list,
                             condition = function(x, .xname) .xname == id,
                             f = function(x) reset_value,
                             how = "replace")

  }
  return(list)
}
