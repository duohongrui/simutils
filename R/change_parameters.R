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
#' @param new_values New values. A list.
#'
#' @export
#'
change_scGAN_parameters <- function(
    project_name,
    new_values
){
  if(!requireNamespace("jsonlite", quietly = TRUE)){
    stop("Package \"jsonlite\" must be installed by \"install.packages('jsonlite')\" command.")
  }
  scGAN_params <- system.file("scGAN_parameters.json", package = "simutils")
  params <- jsonlite::fromJSON(scGAN_params,
                               simplifyDataFrame = TRUE,
                               simplifyVector = FALSE)
  names(params[["experiments"]]) <- project_name
  params <- change_values_in_list(list = params, new_values)
  return(params)
}

#' Change Values From a List
#'
#' @param list The list.
#' @param new_values Other attr ibutes and new values.
#'
#' @export
#'
change_values_in_list <- function(
    list,
    new_values
){
  if(!requireNamespace("rrapply", quietly = TRUE)){
    stop("Package \"rrapply\" must be installed by \"install.packages('rrapply')\" command.")
  }
  reset_params <- new_values
  for(id in names(reset_params)){
    reset_value <- reset_params[[id]]
    if(id == "GPU"){
      list$exp_param$GPU <- list(reset_value)
    }else{
      list <- rrapply::rrapply(list,
                               condition = function(x, .xname) .xname == id,
                               f = function(x) reset_value,
                               how = "replace")
    }
  }
  return(list)
}
