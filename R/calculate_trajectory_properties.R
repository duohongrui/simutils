#' Calculate Metrics for Comparison of Two Trajectories
#'
#' @param ref_data A gene matrix or the trajectory object generated by `dynwrap::wrap_expression` from real data
#' @param ref_data_grouping The labels of cells in real data
#' @param sim_data A gene matrix or the trajectory object generated by `dynwrap::wrap_expression` from simulated data
#' @param sim_data_grouping The labels of cells in simulated data
#' @param seed Random seed
#' @param verbose Whether the messages are returned to users when processing
#' @importFrom dynwrap wrap_expression add_grouping infer_trajectory
#' @importFrom tislingshot ti_slingshot
#'
#' @return A list
#' @export
#'
calculate_trajectory_properties <- function(
  ref_data,
  ref_data_grouping = NULL,
  sim_data,
  sim_data_grouping = NULL,
  seed = 1,
  verbose = TRUE
){
  ### Check data
  if(is.data.frame(ref_data)){
    ref_data <- as.matrix(ref_data)
  }
  if(is.data.frame(sim_data)){
    sim_data <- as.matrix(sim_data)
  }

  ### Build standard objects
  if(!dynwrap::is_wrapper_with_expression(ref_data)){
    ref_data <- dynwrap::wrap_expression(counts = t(ref_data),
                                         expression = log2(t(ref_data) + 1))
  }
  if(!dynwrap::is_wrapper_with_expression(sim_data)){
    sim_data <- dynwrap::wrap_expression(counts = t(sim_data),
                                         expression = log2(t(sim_data) + 1))
  }

  ### grouping
  if(!is.null(ref_data_grouping)){
    ref_data <- dynwrap::add_grouping(dataset = ref_data,
                                      grouping = ref_data_grouping)
  }
  if(!is.null(sim_data_grouping)){
    sim_data <- dynwrap::add_grouping(dataset = sim_data,
                                      grouping = sim_data_grouping)
  }

  ### trajectory inference
  if(!dynwrap::is_wrapper_with_trajectory(ref_data)){
    ref_model <- dynwrap::infer_trajectory(dataset = ref_data,
                                           method = tislingshot::ti_slingshot(),
                                           parameters = NULL,
                                           give_priors = NULL,
                                           seed = seed,
                                           verbose = verbose)
  }else{
    ref_model <- ref_data
  }
  if(!dynwrap::is_wrapper_with_trajectory(sim_data)){
    sim_model <- dynwrap::infer_trajectory(dataset = sim_data,
                                           method = tislingshot::ti_slingshot(),
                                           parameters = NULL,
                                           give_priors = NULL,
                                           seed = seed,
                                           verbose = verbose)
  }else{
    sim_model <- sim_data
  }

  ### Calculate metrics

  result <- calculate_traj_metrics(model_ref = ref_model,
                                   model_sim = sim_model)
  return(result)
}