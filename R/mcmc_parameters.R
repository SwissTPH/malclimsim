#' Create Adaptive Proposal and MCMC Control Parameters
#'
#' This function defines the default values for adaptive proposal control parameters
#' and MCMC control parameters, while allowing users to specify custom values.
#'
#' @param stage possible values are "NULL", stage1", "stage2", or "noadapt"
#' @param initial_vcv_weight Weight for the initial variance-covariance matrix (default = 1).
#' @param initial_scaling Scaling factor for the proposal (default = 2).
#' @param initial_scaling_weight Optional weight for initial scaling (default = NULL).
#' @param min_scaling Minimum scaling factor for the proposal (default = 0.1).
#' @param scaling_increment Increment for scaling factor (default = NULL).
#' @param log_scaling_update Logical, should the scaling be updated on a log scale? (default = TRUE).
#' @param acceptance_target The target acceptance rate for the proposal (default = 0.234).
#' @param forget_rate Rate at which the proposal forgets past history (default = 0.6).
#' @param forget_end Time step at which the forgetting stops (default = Inf).
#' @param adapt_end Time step at which the adaptation stops (default = Inf).
#' @param pre_diminish Time steps before the diminishing starts (default = 40000).
#' @param n_steps Total number of MCMC steps (default = 10000).
#' @param n_burnin Number of burn-in steps (default = 0).
#' @param n_chains Number of MCMC chains (default = 4).
#' @param n_workers Number of workers for parallel execution (default = 4).
#' @param n_threads_total Total number of threads for parallel execution (default = 8).
#'
#' @return A list containing the adaptive proposal parameters and MCMC control parameters.
#'
#' @export
create_mcmc_params <- function(stage = "stage1",
    initial_vcv_weight = 1, initial_scaling = 2, initial_scaling_weight = NULL,
    min_scaling = 0.1, scaling_increment = NULL, log_scaling_update = TRUE,
    acceptance_target = 0.234, forget_rate = 0.6, forget_end = Inf,
    adapt_end = Inf, pre_diminish = 40000,
    n_steps = 10000, n_burnin = 0, n_chains = 4,
    n_workers = 4, n_threads_total = 8
) {
  if(is.null(stage)){
    # Define adaptive proposal control
    adaptive_param <- adaptive_proposal_control(
      initial_vcv_weight = initial_vcv_weight,
      initial_scaling = initial_scaling,
      initial_scaling_weight = initial_scaling_weight,
      min_scaling = min_scaling,
      scaling_increment = scaling_increment,
      log_scaling_update = log_scaling_update,
      acceptance_target = acceptance_target,
      forget_rate = forget_rate,
      forget_end = forget_end,
      adapt_end = adapt_end,
      pre_diminish = pre_diminish
    )

    # Define MCMC control parameters
    control_params <- list(
      n_steps = n_steps,
      n_burnin = n_burnin,
      n_chains = n_chains,
      n_workers = n_workers,
      n_threads_total = n_threads_total
    )

  }

  if(stage == "stage1"){
    adaptive_param <- adaptive_proposal_control(
      initial_vcv_weight = 1,
      initial_scaling = 2,
      initial_scaling_weight = NULL,
      min_scaling = 0.1,
      scaling_increment = NULL,
      log_scaling_update = TRUE,
      acceptance_target = 0.234,
      forget_rate = 0.6,
      forget_end = Inf,
      adapt_end = Inf,
      pre_diminish = 40000,
      n_steps = 40000,
      n_burnin = 0,
      n_chains = 4,
      n_workers = 4,
      n_threads_total = 8
    )
  }

  if(stage == "stage2"){
    adaptive_param <- adaptive_proposal_control(
      initial_vcv_weight = 100,
      initial_scaling = 2,
      initial_scaling_weight = NULL,
      min_scaling = 0.1,
      scaling_increment = NULL,
      log_scaling_update = TRUE,
      acceptance_target = 0.234,
      forget_rate = 0.4,
      forget_end = Inf,
      adapt_end = Inf,
      pre_diminish = 20000,
      n_steps = 40000,
      n_burnin = 0,
      n_chains = 4,
      n_workers = 4,
      n_threads_total = 8
    )
  }

  if(stage == "noadapt"){
    adaptive_param <- adaptive_proposal_control(
      initial_vcv_weight = 5000,
      initial_scaling = 1,
      initial_scaling_weight = NULL,
      min_scaling = 0,
      scaling_increment = NULL,
      log_scaling_update = TRUE,
      acceptance_target = 0.234,
      forget_rate = 0.4,
      forget_end = Inf,
      adapt_end = 20000,
      pre_diminish = 0,
      n_steps = 100000,
      n_burnin = 0,
      n_chains = 1,
      n_workers = 1
    )
  }

  # Return a list with the adaptive proposal and MCMC control parameters
  return(list(
    adaptive_params = adaptive_param,
    control_params = control_params
  ))
}
