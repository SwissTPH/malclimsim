#' Return Default Priors Used for Inference Procedure
#'
#' This function returns a list of default priors used in the inference procedure. Each parameter has
#' an associated prior distribution, an initial value, and range constraints.
#'
#' @return A list where each element corresponds to a parameter. Each parameter element is a list containing:
#' \describe{
#'   \item{initial}{Initial value for the parameter.}
#'   \item{min}{Minimum value for the parameter.}
#'   \item{max}{Maximum value for the parameter.}
#'   \item{prior}{A function defining the prior distribution of the parameter, returning the log probability density.}
#' }
#'
#' @examples
#' # Retrieve the default priors
#' priors <- return_default_priors()
#'
#' # List the names of all parameters with priors
#' names(priors)
#'
#' # Access the prior definition for a specific parameter
#' if ("phi" %in% names(priors)) {
#'   phi_prior <- priors$phi
#'   print(phi_prior)
#'
#'   # Evaluate the log-prior density at a specific value
#'   phi_value <- 0.5
#'   log_density <- phi_prior$prior(phi_value)
#'   print(log_density)
#' }
#' @export
return_default_priors <- function(){
  default_priors <- list(
    ######################################
    ## Parameter estimated in the paper ##
    ######################################
    # Lag of rainfall and temperature
    lag_R = list(initial = 0, min = 0, max = 60, integer = TRUE, prior = function(p) dunif(p, min = 0, max = 60, log = TRUE)),
    lag_T = list(initial = 0, min = 0, max = 60, integer = TRUE, prior = function(p) dunif(p, min = 0, max = 60, log = TRUE)),

    # EIR scaling parameter
    alpha = list(initial = 2, min = 0, max = 20, prior = function(p) dunif(p, min = 0, max = 20, log = TRUE)),

    # Left and right of normal distribution descripting temp-EIR relationship
    sigma_LT = list(initial = 5, min = 2, max = 8, prior = function(p) dnorm(p, mean = 3.2, sd = 1, log = TRUE)), # Gamma distribution to ensure positivity
    sigma_RT = list(initial = 3, min = 2, max = 8, prior = function(p) dnorm(p, mean = 2.0, sd = 1, log = TRUE)), # Gamma distribution to ensure positivity

    # "Optimum" for rainfall-EIR relationship
    R_opt = list(initial = 0, min = -10, max = 10, prior = function(p) dnorm(p, mean = 0, sd = 1, log = TRUE)), # Normal distribution symmetric around 0

    # "Slope" for rainfall-EIR relationship
    k1 = list(initial = 1, min = 0.000001, max = 10, prior = function(p) dunif(p, min = 0.000001, max = 10, log = TRUE)),

    # Dispersion parameter of negative binomial observation model
    size_1 = list(initial = 5.5, min = 0.0001, max = 200, prior = function(p) dunif(p, min = 0.0001, max = 300, log = TRUE)),

    # Seasonal malaria chemoprevention effectiveness
    eff_SMC = list(initial = 0.5, min = 0.00, max = 1, prior = function(p) dunif(p, min = 0, max = 1, log = TRUE)),

    # Population scaling factor
    s = list(initial = 0.8, min = 0.0001, max = 20, prior = function(p) dunif(p, min = 0.0001, max = 20, log = TRUE)),

    b = list(initial = 1, min = 0.001, max = 200, prior = function(p) dunif(p, min = 0.001, max = 200, log = TRUE)),

    # Relative infectivity of symptomatic and asymptomatic individuals
    #qR = list(initial = 0.01, min = 1e-7, max = 1, prior = function(p) dnorm(p, mean = 0.001, sd = 0.02, log = TRUE)),
    #qR = list(initial = 0.01, min = 1e-7, max = 1, prior = function(p) dnorm(p, mean = 0.05, sd = 0.1, log = TRUE)),

    qR = list(initial = 0.01, min = 1e-7, max = 1, prior = function(p) dnorm(p, mean = 0.24, sd = 0.5, log = TRUE)),
    ########################################
    ## Parameters that could be estimated ##
    ########################################
    size_2 = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    kappa_C = list(initial = 40, min = 0.01, max = 2000, prior = function(p) dunif(p, min = 0.01, max = 2000, log = TRUE)),
    kappa_A = list(initial = 40, min = 0.01, max = 2000, prior = function(p) dunif(p, min = 0.01, max = 2000, log = TRUE)),
    z = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    beta_1 = list(initial = 0, min = -100, max = 100, prior = function(p) dunif(p, min = -100, max = 100, log = TRUE)),
    beta_2 = list(initial = 0, min = 0, max = 1, prior = function(p) dunif(p, min = 0, max = 1, log = TRUE)),
    T_opt = list(initial = 26.12, min = 0, max = 40, prior = function(p) dnorm(p, mean = 26.12, sd = 3, log = TRUE)),
    c_X = list(initial = 1, min = 0.001, max = 50, prior = function(p) dunif(p, min = 0.001, max = 50, log = TRUE)),
    fT_C = list(initial = 0.27, min = 0.001, max = 1, prior = function(p) dunif(p, min = 0.001, max = 1, log = TRUE)),

    ###############################################################
    ## All parameters to be fixed must also have a prior defined ##
    ###############################################################
    mu_TS = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    mu_IR = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    mu_RS = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    mu_EI = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    delta_b = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    delta_d = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    delta_a = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    N = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    percAdult = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    pi_s_1 = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    c_s = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    clim_SMC_lag = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    r_C_0 = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    r_A_0 = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE))

  )
  return(default_priors)
}

#' Build a list of mcstate::pmcmc_parameter() objects,
#' starting from the return_default_priors() template and then layering on any overrides.
#'
#' @param param_inputs Named list of initial values (so that we know which parameters actually exist).
#' @param proposal_matrix Numeric matrix (rownames must match parameter names).
#' @param params_to_estimate Character vector of names we actually want to estimate.
#' @param override_priors Optional named list of lists: each element must match the structure
#'                        returned by return_default_priors().
#' @return Named list of \code{mcstate::pmcmc_parameter} objects for all parameters
#'          present in both \code{param_inputs} and \code{proposal_matrix}.
#' @export
build_priors <- function(param_inputs,
                         proposal_matrix,
                         params_to_estimate,
                         override_priors = NULL) {
  # 0. Basic argument checks
  if (is.null(param_inputs) || is.null(proposal_matrix) || is.null(params_to_estimate)) {
    stop("param_inputs, proposal_matrix, and params_to_estimate must all be provided.")
  }

  # 1. Load default priors
  base_priors <- return_default_priors()

  # 2. Overlay any user-supplied overrides
  if (!is.null(override_priors)) {
    base_priors <- modifyList(base_priors, override_priors)
  }

  # 3. Find parameters defined in both inputs and proposal matrix
  common_params <- intersect(names(param_inputs), rownames(proposal_matrix))
  if (length(common_params) == 0) {
    stop("No parameters are common between param_inputs and proposal_matrix.")
  }

  # 4. Ensure all params_to_estimate exist in the common set
  missing_est <- setdiff(params_to_estimate, common_params)
  if (length(missing_est) > 0) {
    stop(sprintf(
      "The following params_to_estimate are not found in both param_inputs and proposal_matrix: %s",
      paste(missing_est, collapse = ", ")
    ))
  }

  # 5. Final set of parameters to build priors for: all common params
  valid_names <- common_params

  # 6. Construct pmcmc_parameter() objects
  priors_list <- vector("list", length(valid_names))
  names(priors_list) <- valid_names
  for (nm in valid_names) {
    if (!nm %in% names(base_priors)) {
      stop(sprintf(
        "No default prior defined for parameter '%s'; please add to return_default_priors() or override_priors.",
        nm
      ))
    }
    info <- base_priors[[nm]]
    priors_list[[nm]] <- mcstate::pmcmc_parameter(
      name    = nm,
      initial = info$initial,
      min     = info$min,
      max     = info$max,
      prior   = info$prior,
      integer = if (!is.null(info$integer)) info$integer else FALSE
    )
  }

  return(priors_list)
}



#' View Details of Specified Priors
#'
#' Returns a data frame summarizing the details of the specified priors, including initial values,
#' bounds, and the functional form of the prior distribution.
#'
#' @param param_inputs A named list of base parameter values. Used if \code{priors} is not supplied.
#' @param proposal_matrix Covariance matrix used in proposal distribution. Required if \code{priors} is NULL.
#' @param params_to_estimate A character vector specifying which parameters' priors should be viewed.
#' @param priors Optional list of priors. If \code{NULL}, the function calls \code{initialize_priors()} to generate defaults.
#'
#' @return A data frame with one row per parameter and the following columns:
#' \describe{
#'   \item{Name}{The parameter name.}
#'   \item{Initial}{The initial value used in inference.}
#'   \item{Min}{The lower bound of the prior.}
#'   \item{Max}{The upper bound of the prior.}
#'   \item{Description}{A textual representation of the prior function.}
#' }
#'
#' @examples
#' # Define dummy prior functions
#' dummy_prior <- function(p) dnorm(p, mean = 0, sd = 1, log = TRUE)
#'
#' # Create a minimal list of priors matching the expected structure
#' dummy_priors <- list(
#'   a_R = list(name = "a_R", initial = 0.1, min = -2, max = 2, prior = dummy_prior),
#'   b_R = list(name = "b_R", initial = 1.5, min = 0, max = 3, prior = dummy_prior),
#'   qR  = list(name = "qR",  initial = 0.01, min = 0, max = 1, prior = dummy_prior),
#'   z   = list(name = "z",   initial = 0.2, min = 0, max = 1, prior = dummy_prior),
#'   eff_SMC = list(name = "eff_SMC", initial = 0.5, min = 0, max = 1, prior = dummy_prior),
#'   phi = list(name = "phi", initial = 0.3, min = -1, max = 1, prior = dummy_prior),
#'   size = list(name = "size", initial = 5, min = 0.01, max = 20, prior = dummy_prior)
#' )
#'
#' # View priors for a subset of parameters using a supplied list
#' view_priors(
#'   param_inputs = NULL,
#'   proposal_matrix = NULL,
#'   params_to_estimate = c("a_R", "b_R", "qR", "z", "eff_SMC", "phi", "size"),
#'   priors = dummy_priors
#' )
#' @export
view_priors <- function(param_inputs, proposal_matrix, params_to_estimate, priors = NULL) {
  if (is.null(priors)) {
    priors <- initialize_priors(param_inputs, proposal_matrix, params_to_estimate)
  }

  # Filter the priors to only include those specified in params_to_estimate
  filtered_priors <- priors[names(priors) %in% params_to_estimate]

  # Convert the filtered priors list to a data frame for easy viewing
  priors_df <- do.call(rbind, lapply(filtered_priors, function(prior) {
    data.frame(
      Name = prior$name,
      Initial = prior$initial,
      Min = prior$min,
      Max = prior$max,
      Description = paste(deparse(prior$prior), collapse = " ")
    )
  }))

  return(priors_df)
}


#' Plot Default Priors
#'
#' Plots the default or user-specified prior distributions for the parameters of interest.
#'
#' @param param_inputs A list of model inputs (used only if default priors are initialized).
#' @param proposal_matrix A covariance matrix used to generate default priors if `priors` is not supplied.
#' @param params_to_estimate A character vector specifying which parameters to include in the plot.
#' @param priors A named list of priors. Each prior should be a list containing `min`, `max`, and a `prior` function
#'        that returns the log-density for a given numeric input. If `NULL`, default priors will be computed.
#'
#' @return A `ggplot` object showing the prior distributions of the selected parameters.
plot_priors <- function(param_inputs, proposal_matrix, params_to_estimate, priors = NULL) {
  if (is.null(priors)) {
    priors <- initialize_priors(param_inputs, proposal_matrix, params_to_estimate)
  }

  # Filter the priors to only include those specified in params_to_estimate
  filtered_priors <- priors[names(priors) %in% params_to_estimate]

  # Initialize an empty list to collect plot data
  plot_data_list <- list()

  # Loop over the filtered priors to generate plot data
  for (prior_name in names(filtered_priors)) {
    prior <- filtered_priors[[prior_name]]

    # Check if min and max are finite numbers
    if (is.finite(prior$min) && is.finite(prior$max)) {
      param_range <- seq(prior$min, prior$max, length.out = 1000)
      prior_values <- exp(prior$prior(param_range))  # Ensure prior is on original scale (not log)

      # Collect the data for plotting
      plot_data_list[[prior_name]] <- data.frame(
        Parameter = prior_name,
        Value = param_range,
        Density = prior_values
      )
    } else {
      warning(paste("Skipping prior for", prior_name, "due to non-finite min or max values."))
    }
  }

  # Combine all plot data into a single data frame
  plot_data <- do.call(rbind, plot_data_list)

  # Create the plot using ggplot2
  prior_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Value, y = Density)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ Parameter, scales = "free", ncol = 2) +
    ggplot2::labs(title = "Prior Distributions for Specified Parameters",
                  x = "Parameter Value",
                  y = "Density") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      strip.text = ggplot2::element_text(size = 14, face = "bold")
    )

  return(prior_plot)
}


#' Update Default Priors
#'
#' Updates the default prior specifications for a set of parameters, allowing users to replace or modify
#' elements such as distribution type, bounds, and initial values.
#'
#' @param param_inputs A named list of model parameter inputs used to initialize default priors.
#' @param proposal_matrix A proposal covariance matrix used to set up default priors.
#' @param params_to_estimate A character vector of parameters to estimate and update priors for.
#' @param priors Optional. A list of existing priors. If NULL (default), default priors are generated using \code{initialize_priors()}.
#' @param new_priors A named list of new prior specifications. Each element should be a list with updated fields (e.g., \code{initial}, \code{min}, \code{max}, \code{prior}).
#'
#' @return A list of prior specifications with selected entries updated based on \code{new_priors}.
#'
#' @export
update_priors <- function(param_inputs, proposal_matrix, params_to_estimate, priors = NULL, new_priors) {
  if (is.null(priors)) {
    priors <- initialize_priors(param_inputs, proposal_matrix, params_to_estimate)
  }

  for (param in names(new_priors)) {
    if (param %in% names(priors)) {
      priors[[param]] <- modifyList(priors[[param]], new_priors[[param]])
    } else {
      warning(paste("Parameter", param, "not found in priors. Ignoring."))
    }
  }

  return(priors)
}
