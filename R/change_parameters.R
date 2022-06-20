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
