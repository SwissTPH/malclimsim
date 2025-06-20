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
                               n_workers = 4, n_threads_total = 8)
{
  if(is.null(stage)){
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
      pre_diminish = 40000
    )

    control_params <- list(n_steps = 40000, n_burnin = 0, n_chains = 4, n_workers = 4, n_threads_total = 8)
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
      pre_diminish = 20000
    )

    control_params <- list(n_steps = 40000, n_burnin = 0, n_chains = 4, n_workers = 4, n_threads_total = 8)
  }

  if(stage == "noadapt"){
    adaptive_param <- adaptive_proposal_control(
      initial_vcv_weight = 500,
      initial_scaling = 1,
      initial_scaling_weight = NULL,
      min_scaling = 0,
      scaling_increment = NULL,
      log_scaling_update = TRUE,
      acceptance_target = 0.234,
      forget_rate = 0.4,
      forget_end = Inf,
      adapt_end = 20000,
      pre_diminish = 0
    )

    control_params <- list(n_steps = 100000, n_burnin = 20000, n_chains = 3, n_workers = 3, n_threads_total = 6)
  }

  return(list(
    adaptive_params = adaptive_param,
    control_params = control_params
  ))
}



#' Create Starting Values for Parameters
#'
#' This function generates starting values for parameters based on the specified method. It can either draw random values from a uniform distribution or assign values from a specified interval.
#'
#' @param params_to_estimate A character vector of parameter names that are to be estimated.
#' @param control_params A list of control parameters (must include `n_chains` to specify the number of chains).
#' @param min_max_start_values A named list with minimum and maximum start values for each parameter.
#' @param random A logical value. If `TRUE`, random values are drawn from a uniform distribution; if `FALSE`, values are evenly spaced between the lower and upper bounds.
#' @param seed An integer specifying the seed for random number generation (default is 10).
#' @param model A model object with a `param()` method that returns all parameter names.
#' @param param_inputs A list of parameter values passed to initialize the model, used to extract parameters.
#'
#' @return A matrix of starting values for each parameter and chain.
#'
#' @export
create_start_values <- function(params_to_estimate, control_params, min_max_start_values = NULL,
                                random = TRUE, seed = 10, model, param_inputs) {

  # Retrieve all parameter names from the model
  model_instance <- model$new(pars = param_inputs, time = 0, n_particles = 1)
  param_names <- intersect(names(param_inputs), names(model_instance$param()))
  param_names <- param_names[sapply(param_inputs[param_names], length) == 1]

  # Filter params_to_estimate to include only model parameters
  params_to_estimate <- intersect(params_to_estimate, param_names)

  # Default min_max_start_values if NULL, providing predefined bounds for each parameter
  if (is.null(min_max_start_values)) {
    min_max_start_values <- list(
      ######################################
      ## Parameter estimated in the paper ##
      ######################################
      lag_R = c(15, 30),
      lag_T = c(15,30),
      alpha = c(0.2, 0.8),
      sigma_LT = c(5, 9),
      sigma_RT = c(4, 9),
      R_opt = c(0.5, 0.9),
      k1 = c(0.1, 0.8),
      size_1 = c(5, 30),
      eff_SMC = c(0.2, 0.8),
      s = c(0.01, 10),
      qR = c(0.001, 0.5),

      #################################################
      ## Parameter estimated that could be estimated ##
      #################################################
      size_2 = c(5, 30),
      kappa_C = c(70, 200),
      kappa_A = c(70, 200),
      z = c(0.1, 0.7),
      beta_1 = c(-5, 5),
      beta_2 = c(0, 1),
      T_opt = c(24, 27),
      b = c(22, 32),
      fT_C = c(0.1, 0.7)
      )
  }

  # Set up the number of chains and prepare a matrix for starting values
  n_chains <- control_params$n_chains
  start_values <- matrix(NA, nrow = length(params_to_estimate), ncol = n_chains)
  rownames(start_values) <- params_to_estimate

  # Set seed for reproducibility
  set.seed(seed)

  # Generate the starting values
  for (name in params_to_estimate) {
    # bounds <- min_max_start_values[[name]] %||% c(0, 1)  # Default to (0, 1) if not specified
    bounds <- if (is.null(min_max_start_values[[name]])) c(0, 1) else min_max_start_values[[name]]
    if (random) {
      start_values[name, ] <- runif(n_chains, min = bounds[1], max = bounds[2])
    } else {
      start_values[name, ] <- seq(bounds[1], bounds[2], length.out = n_chains)
    }
  }

  return(start_values)
}


#' Create a Proposal Matrix
#'
#' This function generates a proposal matrix for use in MCMC methods, where each parameter has its own proposal variance.
#'
#' @param params_to_estimate A character vector of parameter names to estimate.
#' @param proposal_variance A named list containing the variance for each parameter. If `NULL`, default values are used.
#' @param model A model object with a `param()` method that returns all parameter names.
#' @param param_inputs A list of parameter values passed to initialize the model, used to extract parameters.
#'
#' @return A square matrix with dimensions equal to the number of parameters, with diagonal elements corresponding to the proposal variances for each parameter.
#'
#' @export
create_proposal_matrix <- function(params_to_estimate, proposal_variance = NULL, correlation_matrix = NULL, model, param_inputs) {

  # Retrieve all parameter names from the model
  model_instance <- model$new(pars = param_inputs, time = 0, n_particles = 1)
  param_names <- intersect(names(param_inputs), names(model_instance$param()))
  param_names <- param_names[sapply(param_inputs[param_names], length) == 1]

  # Filter params_to_estimate to include only model parameters
  params_to_estimate <- intersect(params_to_estimate, param_names)

  # Default proposal variance if NULL, providing predefined variances for each parameter
  if (is.null(proposal_variance)) {
    proposal_variance <- list(
      ######################################
      ## Parameter estimated in the paper ##
      ######################################
      lag_T = 800, lag_R = 800, alpha = 4,
      sigma_LT = 3, sigma_RT = 3, R_opt = 2,
      k1 = 3, size_1 = 15, eff_SMC = 0.5, s = 0.3,
      qR = 1,

      #################################################
      ## Parameter estimated that could be estimated ##
      #################################################
      z = 0.1, size_2 = 5, fT_C = 0.1, T_opt = 5,
      kappa_C = 50, kappa_A = 50, beta_1 = 2,
      beta_2 = 5
    )
  }

  # Initialize the proposal matrix with small default variance (0.01) for all parameters
  proposal_matrix <- diag(0.01, length(param_names))
  rownames(proposal_matrix) <- param_names
  colnames(proposal_matrix) <- param_names

  # Update the diagonal elements with the specified proposal variances
  for (name in params_to_estimate) {
    #proposal_matrix[name, name] <- proposal_variance[[name]] %||% 0.5  # Default to 0.1 if not specified
    proposal_matrix[name, name] <- if (is.null(proposal_variance[[name]])) 0.5 else proposal_variance[[name]]

  }

  # If a correlation matrix is provided, integrate it into the proposal matrix
  if (!is.null(correlation_matrix)) {
    corr_params <- intersect(rownames(correlation_matrix), param_names)

    if (length(corr_params) > 0) {
      # Extract standard deviations for the correlated parameters
      std_devs <- sqrt(diag(proposal_matrix[corr_params, corr_params, drop = FALSE]))

      # Construct the covariance matrix for the provided correlations
      cov_matrix <- diag(std_devs) %*% correlation_matrix[corr_params, corr_params, drop = FALSE] %*% diag(std_devs)

      # Place the covariance matrix into the larger proposal matrix
      proposal_matrix[corr_params, corr_params] <- cov_matrix
    } else {
      warning("No matching parameters found between the correlation matrix and model parameters.")
    }
  }

  return(proposal_matrix)
}


#' Extract Variance-Covariance Matrix and Restart Values
#'
#' This function extracts the variance-covariance matrix and restart values from the results of an MCMC run.
#' It can save the computed proposal matrix and start values to specified file paths.
#'
#' @param results list: Results of the `inf_run` function containing MCMC parameters and chain identifiers.
#' @param S_prev numeric: Number of previous steps in the MCMC chains to include in the calculation.
#' @param save boolean: Whether or not to save the computed results to files (default: TRUE).
#' @param param_names list: List of model parameter names, used for labeling the proposal matrix.
#' @param file_proposal character: File path to save the proposal covariance matrix (if `save = TRUE`).
#' @param file_start character: File path to save the start values (if `save = TRUE`).
#'
#' @return A list containing:
#' \describe{
#'   \item{start_values}{A matrix of restart values (medians of the selected MCMC steps) for each parameter.}
#'   \item{proposal_matrix}{The variance-covariance matrix for use as the proposal distribution.}
#' }
#'
#' @examples
#' # Example usage:
#' results <- list(
#'   mcmc_run = list(
#'     pars = matrix(runif(1000), ncol = 10),
#'     chain = rep(1:5, each = 200)
#'   )
#' )
#' param_names <- paste0("param_", 1:10)
#' S_prev <- 100
#' file_proposal <- "proposal_matrix.rds"
#' file_start <- "start_values.rds"
#'
#' output <- extract_vcv(
#'   results = results,
#'   S_prev = S_prev,
#'   save = FALSE,
#'   param_names = param_names,
#'   file_proposal = file_proposal,
#'   file_start = file_start
#' )
#'
#' print(output)
#'
#' @export
extract_vcv <- function(results, S_prev, save = TRUE, param_names,
                        file_proposal = "", file_start = "") {
  # Extract the MCMC parameters and chain identifiers
  mcmc_pars <- results$mcmc_run$pars
  mcmc_chains <- results$mcmc_run$chain

  # Get the number of chains
  unique_chains <- unique(mcmc_chains)

  # Select S_prev subset from each chain
  selected_pars_list <- lapply(unique_chains, function(chain) {
    chain_indices <- which(mcmc_chains == chain)
    chain_subset <- tail(chain_indices, S_prev)
    mcmc_pars[chain_subset, ]
  })

  selected_pars <- do.call(rbind, selected_pars_list)

  # Compute the median across all chains for each parameter
  restart_values <- apply(selected_pars, 2, median)

  start_values <- matrix(NA, nrow = length(restart_values), ncol = 1)
  start_values[, 1] <- restart_values
  rownames(start_values) <- colnames(selected_pars)

  # Compute the variance-covariance matrix considering all chains
  int_vcv <- cov(selected_pars)

  proposal_matrix <- diag(0.1, length(param_names))
  colnames(proposal_matrix) <- param_names
  rownames(proposal_matrix) <- param_names

  for (row_name in rownames(int_vcv)) {
    for (col_name in colnames(int_vcv)) {
      i <- which(rownames(proposal_matrix) == row_name)
      j <- which(colnames(proposal_matrix) == col_name)
      proposal_matrix[i, j] <- int_vcv[row_name, col_name]
    }
  }

  if (save) {
    saveRDS(proposal_matrix, file_proposal)
    saveRDS(start_values, file_start)
  }

  return(list(start_values, proposal_matrix))
}



#' Reorder MCMC Start Values to Match Prior Specification
#'
#' Ensures that the rows of a matrix of MCMC start values are ordered consistently with the names in the `param_priors` list. This is useful when passing initial values to an MCMC routine that expects a specific parameter order.
#'
#' @param start_values A matrix of initial parameter values for MCMC, where each row corresponds to a parameter and each column to a chain or initialization.
#'                     The row names should be the parameter names.
#' @param param_priors A named list of prior specifications (e.g., as used by `mcstate::pmcmc_parameters$new()`), where names correspond to parameter names.
#'
#' @return A reordered matrix of start values with the row order matching the order of names in `param_priors`. Parameters missing from `start_values` are dropped (with a warning optionally).
#'
#' @details
#' The function matches parameter names from `start_values` to those in `param_priors`, and reorders them accordingly.
#' Any parameters present in `param_priors` but missing in `start_values` are ignored in the returned matrix.
#'
#' @examples
#' start_vals <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, dimnames = list(c("beta", "gamma"), NULL))
#' priors <- list(gamma = list(min = 0, max = 1, prior = dunif),
#'                beta = list(min = 0, max = 1, prior = dunif))
#' reordered <- reorder_start_values(start_vals, priors)
#'
#' @export
reorder_start_values <- function(start_values, param_priors) {
  # Get the order of parameters based on names(param_priors)
  prior_order <- names(param_priors)

  # Extract the names of the parameters from start_values
  start_values_params <- rownames(start_values)

  # Match the parameters to ensure the correct ordering
  ordered_params <- intersect(prior_order, start_values_params)

  # Reorder start_values matrix to match the order in prior_order
  reordered_start_values <- start_values[ordered_params, , drop = FALSE]

  # Optional warning for missing parameters
  # missing_params <- setdiff(prior_order, ordered_params)
  # if(length(missing_params) > 0) {
  #   warning(paste("The following parameters are missing from start_values:", paste(missing_params, collapse = ", ")))
  # }

  return(reordered_start_values)
}


