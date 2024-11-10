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
      mu_EI = c(1/14, 1/8), qR = c(0.001, 0.2), a_R = c(0.4, 0.8), b_R = c(1, 3),
      eff_SMC = c(0.2, 0.8), s = c(0.2, 0.8), phi = c(0.2, 0.6), k_par = c(0.6, 0.8),
      delta_temp = c(-4, 4), mu_RS_C = c(1/350, 1/250), z = c(0.2, 0.6), z_A = c(0.2, 0.6),
      z_C2 = c(0.2, 0.6), rho = c(0.3, 0.9), eta = c(0.1, 0.9), size = c(6, 7),
      phi_C2 = c(0.2, 0.6), phi_A = c(0.2, 0.6), tau = c(0.2, 0.6), p_surv = c(0.89, 0.92),
      mu_IR = c(1/10, 1/2), shift1 = c(1, 30), shift2 = c(1, 30), kappa = c(0.2, 0.6)
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
    bounds <- min_max_start_values[[name]] %||% c(0, 1)  # Default to (0, 1) if not specified
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
create_proposal_matrix <- function(params_to_estimate, proposal_variance = NULL, model, param_inputs) {

  # Retrieve all parameter names from the model
  model_instance <- model$new(pars = param_inputs, time = 0, n_particles = 1)
  param_names <- intersect(names(param_inputs), names(model_instance$param()))
  param_names <- param_names[sapply(param_inputs[param_names], length) == 1]

  # Filter params_to_estimate to include only model parameters
  params_to_estimate <- intersect(params_to_estimate, param_names)

  # Default proposal_variance if NULL, providing predefined variances for each parameter
  if (is.null(proposal_variance)) {
    proposal_variance <- list(
      mu_EI = 0.1, qR = 0.1, a_R = 0.2, b_R = 0.2, eff_SMC = 0.1, s = 0.3,
      phi = 0.1, k_par = 0.1, delta_temp = 0.1, mu_RS_C = 0.1, z = 0.1,
      z_A = 0.1, z_C2 = 0.1, rho = 0.1, eta = 0.1, size = 0.1, phi_C2 = 0.1,
      phi_A = 0.1, tau = 0.1, p_surv = 0.1, mu_IR = 0.1, shift1 = 0.1,
      shift2 = 0.1, kappa = 0.1
    )
  }

  # Initialize the proposal matrix with small default variance (0.01) for all parameters
  proposal_matrix <- diag(0.01, length(param_names))
  rownames(proposal_matrix) <- param_names
  colnames(proposal_matrix) <- param_names

  # Update the diagonal elements with the specified proposal variances
  for (name in params_to_estimate) {
    proposal_matrix[name, name] <- proposal_variance[[name]] %||% 0.1  # Default to 0.1 if not specified
  }

  return(proposal_matrix)
}

reorder_start_values <- function(start_values, param_priors) {
  # Get the order of parameters based on names(param_priors)
  prior_order <- names(param_priors)

  # Extract the names of the parameters from start_values
  start_values_params <- rownames(start_values)

  # Match the parameters to ensure the correct ordering
  ordered_params <- intersect(prior_order, start_values_params)

  # Reorder start_values matrix to match the order in prior_order
  reordered_start_values <- start_values[ordered_params, , drop = FALSE]

  # Ensure that the order matches exactly the appearance of names(param_priors)
  missing_params <- setdiff(prior_order, ordered_params)

  #if(length(missing_params) > 0) {
  #  warning(paste("The following parameters are missing from start_values:", paste(missing_params, collapse = ", ")))
  #}

  return(reordered_start_values)
}

