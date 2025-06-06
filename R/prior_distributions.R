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
#' # Get default priors
#' priors <- return_default_priors()
#'
#' # Access the prior for a specific parameter
#' phi_prior <- priors$phi
#' print(phi_prior)
#'
#' # Evaluate the prior for phi at a specific value
#' phi_value <- 0.5
#' log_prior_density <- phi_prior$prior(phi_value)
#' print(log_prior_density)
#'
#' @export
return_default_priors <- function(){
  default_priors <- list(

    # 1. Transmission and Recovery Parameters
    #phi = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 40, 12, log = TRUE)),
    phi = list(initial = 0.2, min = 0.01, max = 2, prior = function(p) dunif(p, min = 0.01, max = 2, log = TRUE)),
    phi_1 = list(initial = 0.2, min = 0.01, max = 0.999, prior = function(p) dnorm(p, mean = 0.75, sd = 0.05, log = TRUE)),
    c_phi = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    #qR = list(initial = 0.2, min = 1e-6, max = 0.5, prior = function(p) dunif(p, min = 1e-7, max = 0.5, log = TRUE)),
    qR = list(initial = 0.2, min = 1e-6, max = 0.5, prior = function(p) dnorm(p, 0.24, 0.5, log = TRUE)),
    qR1 = list(initial = 0.2, min = 1e-6, max = 1, prior = function(p) dnorm(p, 0.24, 0.05, log = TRUE)),
    #qR2 = list(initial = 0.2, min = 1e-6, max = 1, prior = function(p) dunif(p, min = 0.0001, max = 1, log = TRUE)),
    #qR2 = list(initial = 0.001, min = 1e-6, max = 0.05, prior = function(p) dunif(p, min = 1e-6, max = 0.05, log = TRUE)),
    c_qR = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    #qR2 = list(initial = 0.01, min = 1e-7, max = 1, prior = function(p) dnorm(p, mean = 0.5, sd = 0.02, log = TRUE)),
    qR2 = list(initial = 0.01, min = 1e-7, max = 1, prior = function(p) dnorm(p, mean = 0.001, sd = 0.02, log = TRUE)),
    #mu_RS_C = list(initial = 1/200, min = 1/400, max = 1/120, prior = function(p) dgamma(p, shape = 2, rate = 480, log = TRUE)),
    mu_RS_C = list(initial = 1/150, min = 1/199, max = 1, prior = function(p) dunif(p, min = 1/199, max = 1, log = TRUE)),
    mu_EI = list(initial = 1/8, min = 1/15, max = 1/6, prior = function(p) dgamma(p, shape = 6, rate = 63, log = TRUE)),
    phi_C2 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 40, 12, log = TRUE)),
    phi_A = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 40, 12, log = TRUE)),
    mu_TS = list(initial = 0.2, min = 0, max = 2, prior = function(p) dunif(p, min = 0, max = 2, log = TRUE)),
    mu_IR = list(initial = 1/5, min = 0.001, max = 1, prior = function(p) dunif(p, min = 0.001, max = 1, log = TRUE)),
    fT_C = list(initial = 0.27, min = 0.001, max = 1, prior = function(p) dunif(p, min = 0.001, max = 1, log = TRUE)),
    #s_1 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    s_1 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dnorm(p, mean = 0.875, sd = 0.05, log = TRUE)),
    #c_s = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    c_s = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dnorm(p, mean = 0.85, sd = 0.2, log = TRUE)),
    #c_s = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dnorm(p, mean = 0.5, sd = 0.25, log = TRUE)),
    eta = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),

    # 2. Survival and Population Parameters
    a_R = list(initial = 0.5, min = 0.01, max = 1, prior = function(p) dbeta(p, 6, 12, log = TRUE)),
    b_R = list(initial = 3, min = -Inf, max = Inf, prior = function(p) dnorm(p, mean = 2, sd = 0.25, log = TRUE)),
    p_surv = list(initial = 0.91, min = 0.88, max = 0.97, prior = function(p) dnorm(p, mean = 0.91, sd = 0.01, log = TRUE)),
    #size = list(initial = 5.5, min = 1, max = 50, prior = function(p) dnorm(p, mean = 5.5, sd = 0.5, log = TRUE)),
    size = list(initial = 5.5, min = 1, max = 50, prior = function(p) dnorm(p, mean = 5.5, sd = 0.5, log = TRUE)),
    size_1 = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    size_2 = list(initial = 5.5, min = 0.01, max = 200, prior = function(p) dunif(p, min = 0.01, max = 300, log = TRUE)),
    kappa_C = list(initial = 40, min = 0.01, max = 2000, prior = function(p) dunif(p, min = 0.01, max = 2000, log = TRUE)),
    kappa_A = list(initial = 40, min = 0.01, max = 2000, prior = function(p) dunif(p, min = 0.01, max = 2000, log = TRUE)),
    # 3. Intervention and Effectiveness Parameters

    #eff_SMC = list(initial = 0.6, min = 0.00, max = 1, prior = function(p) dbeta(p, 2, 2, log = TRUE)),
    eff_SMC = list(initial = 0.5, min = 0.00, max = 1, prior = function(p) dunif(p, min = 0, max = 1, log = TRUE)),
    #z = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 4, 2, log = TRUE)),
    z = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    p_MH_C = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    rho = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    w1 = list(initial = 0, min = 0, max = 0.2, prior = function(p) dunif(p, min = 0, max = 0.2, log = TRUE)),
    w2 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),

    # Multiplicative constants (k parameters)
    #k1 = list(initial = 1, min = 0.01, max = 5, prior = function(p) dnorm(p, mean = 1, sd = 0.1, log = TRUE)),
    k2 = list(initial = 1, min = 0.5, max = 1.5, prior = function(p) dnorm(p, mean = 1, sd = 0.1, log = TRUE)),
    k3 = list(initial = 1, min = 0.5, max = 1.5, prior = function(p) dnorm(p, mean = 1, sd = 0.1, log = TRUE)),
    k4 = list(initial = -0.00924, min = -1, max = 1.5, prior = function(p) dnorm(p, mean = -0.00924, sd = 0.01, log = TRUE)),
    k5 = list(initial = 0.453, min = 0.1, max = 2, prior = function(p) dnorm(p, mean = 0.453, sd = 0.1, log = TRUE)),
    k7 = list(initial = 0.000112, min = 0.00001, max = 0.0003, prior = function(p) dnorm(p, mean = 0.000112, sd = 0.00001, log = TRUE)),

    # Additive constants (c parameters)
    c1 = list(initial = 0, min = -0.5, max = 0.5, prior = function(p) dnorm(p, mean = 0, sd = 0.1, log = TRUE)),
    c3 = list(initial = 0, min = -0.5, max = 0.5, prior = function(p) dnorm(p, mean = 0, sd = 0.1, log = TRUE)),
    c4 = list(initial = 0, min = -1, max = 8, prior = function(p) dnorm(p, mean = 4.77, sd = 0.1, log = TRUE)),
    c5 = list(initial = 14.7, min = -0.5, max = 30, prior = function(p) dnorm(p, mean = 14.7, sd = 0.1, log = TRUE)),
    c6 = list(initial = 34, min = 0.8, max = 45, prior = function(p) dnorm(p, mean = 34, sd = 5, log = TRUE)),
    c7 = list(initial = 15.384, min = 5, max = 30, prior = function(p) dnorm(p, mean = 15.384, sd = 4, log = TRUE)),
    c8 = list(initial = 35, min = 20, max = 50, prior = function(p) dnorm(p, mean = 35, sd = 5, log = TRUE)),
    c9 = list(initial = 0.01, min = 0.005, max = 0.02, prior = function(p) dnorm(p, mean = 0.01, sd = 0.005, log = TRUE)),

    # 4. Population Proportions and Initial Conditions
    s = list(initial = 0.8, min = 0.01, max = 20, prior = function(p) dunif(p, min = 0.01, max = 20, log = TRUE)),
    N = list(initial = 0.2, min = 0, max = 500000, prior = function(p) dunif(p, min = 0, max = 500000, log = TRUE)),
    percAdult = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    percC1 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    percC2 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    SC0 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    EC0 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    IC0 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    TC0 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    RC0 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    SA0 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    EA0 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    IA0 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    TA0 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    RA0 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),

    # 5. Temperature and Environmental Parameters
    shift1 = list(initial = 5, min = 0, max = 90, prior = function(p) dunif(p, min = 0, max = 90, log = TRUE)),
    shift2 = list(initial = 5, min = 0, max = 90, prior = function(p) dunif(p, min = 0, max = 90, log = TRUE)),
    scale = list(initial = 0.2, min = 0.01, max = 1.5, prior = function(p) dunif(p, min = 0.01, max = 1.5, log = TRUE)),
    loc = list(initial = 0.2, min = 0, max = 5, prior = function(p) dunif(p, min = 0, max = 5, log = TRUE)),
    mean_t = list(initial = 20, min = 10, max = 40, prior = function(p) dunif(p, min = 10, max = 40, log = TRUE)),
    k_par = list(initial = 0.7, min = 0.4, max = 1, prior = function(p) dunif(p, min = 0.4, max = 1, log = TRUE)),
    delta_temp = list(initial = 0, min = -10, max = 2, prior = function(p) dunif(p, min = -5, max = 2, log = TRUE)),

    # 6. Adjustment and Scaling Parameters
    tau = list(initial = 0.2, min = 0, max = 1, prior = function(p) dbeta(p, 40, 12, log = TRUE)),
    kappa = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 0.125, 0.125, log = TRUE)),
    z_A = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 0.125, 0.125, log = TRUE)),
    z_C2 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 0.125, 0.125, log = TRUE)),
    D <- list(initial = 2, min = 1, max = 200, integer = TRUE, prior = function(p) dunif(p, min = 1, max = 200, log = TRUE)),
    lag_R = list(initial = 0, min = 0, max = 60, integer = TRUE, prior = function(p) dunif(p, min = 0, max = 60, log = TRUE)),
    lag_T = list(initial = 0, min = 0, max = 60, integer = TRUE, prior = function(p) dunif(p, min = 0, max = 60, log = TRUE)),
    lag_SMC = list(initial = 0, min = 0, max = 30, integer = TRUE, prior = function(p) dunif(p, min = 0, max = 30, log = TRUE)),
    # Gamma distribution to ensure positivity
    alpha = list(initial = 0.5, min = 0, max = 1, prior = function(p) dunif(p, min = 0, max = 1, log = TRUE)),
    #beta_1 = list(initial = 0, min = 0, max = 1, prior = function(p) dunif(p, min = 0, max = 1, log = TRUE)),
    beta_1 = list(initial = 0, min = -100, max = 100, prior = function(p) dunif(p, min = -100, max = 100, log = TRUE)),
    beta_2 = list(initial = 0, min = 0, max = 1, prior = function(p) dunif(p, min = 0, max = 1, log = TRUE)),

    # T_opt - Must be greater than 0, likely between 24 and 32
    T_opt = list(initial = 26.12, min = 0, max = 40, prior = function(p) dnorm(p, mean = 26.12, sd = 3, log = TRUE)),

    # R_opt - Can be negative, likely between -5 and 5
    R_opt = list(initial = 0, min = -10, max = 10, prior = function(p) dnorm(p, mean = 0, sd = 1, log = TRUE)), # Normal distribution symmetric around 0

    # k1 - Can be negative but likely positive and less than 2
    k1 = list(initial = 1, min = 0.000001, max = 10, prior = function(p) dunif(p, min = 0.000001, max = 10, log = TRUE)),
    #k1 = list(initial = 1, min = -5, max = 5, prior = function(p) dnorm(p, mean = 1, sd = 1, log = TRUE)), # Normal distribution with bias towards positive values

    # sigma_T - Must be positive, value is uncertain
    sigma_LT = list(initial = 5, min = 3.5, max = 7, prior = function(p) dnorm(p, mean = 5, sd = 0.5, log = TRUE)), # Gamma distribution to ensure positivity
    sigma_RT = list(initial = 3, min = 2, max = 7, prior = function(p) dnorm(p, mean = 3, sd = 2, log = TRUE)), # Gamma distribution to ensure positivity

    #b = list(initial = 1, min = 0.001, max = 50, prior = function(p) dgamma(p, shape = 2, rate = 1, log = TRUE))
    b = list(initial = 1, min = 0.001, max = 50, prior = function(p) dunif(p, min = 0.001, max = 50, log = TRUE))
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
#'                        returned by return_default_priors(). If you supply an element "phi" here,
#'                        it replaces the entire \code{phi} block in the default.
#' @return Named list of \code{mcstate::pmcmc_parameter} objects, one per parameter in \code{params_to_estimate}.
#' @export
build_priors <- function(param_inputs,
                         proposal_matrix,
                         params_to_estimate,
                         override_priors = NULL) {
  # 1. Grab the “package defaults”:
  base_priors <- return_default_priors()

  # 2. If the user supplied overrides, overlay them:
  if (!is.null(override_priors)) {
    # modifyList() merges named elements, replacing any named components.
    base_priors <- modifyList(base_priors, override_priors)
  }

  # 3. Only keep those parameters that are actually in param_inputs ∩ params_to_estimate ∩ proposal_matrix
  all_names   <- intersect(names(param_inputs), rownames(proposal_matrix))
  valid_names <- intersect(all_names, params_to_estimate)

  # 4. Build pmcmc_parameter() for each valid parameter:
  priors_list <- list()
  for (nm in valid_names) {
    if (! nm %in% names(base_priors)) {
      # If someone forgot to specify a prior in the “defaults,” you could throw an error or fallback:
      stop(paste0("No default prior found for parameter ‘", nm, "’."))
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



#' View Default Priors
#'
#' This function displays the default priors and their details, including initial values, min/max bounds, and prior distributions.
#'
#' @param params_to_estimate Character vector specifying which parameters' priors should be viewed.
#' @param priors List of priors to view. If NULL, the default priors will be displayed.
#' @return A data frame containing the details of each prior specified in params_to_estimate.
#' @export
#' @examples
#' view_priors(param_inputs, proposal_matrix, params_to_estimate = c("a_R", "b_R", "qR", "z", "eff_SMC", "phi", "size"))
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
#' This function plots the default priors for the specified parameters, showing their distributions.
#'
#' @param param_inputs List of input parameters for the model.
#' @param proposal_matrix Covariance matrix for the MCMC proposal distribution.
#' @param params_to_estimate Character vector specifying which parameters' priors should be plotted.
#' @param priors List of priors to plot. If NULL, the default priors will be displayed.
#' @return A ggplot object showing the prior distributions for each of the specified parameters.
#' @export
#' @examples
#' plot_priors(param_inputs, proposal_matrix, params_to_estimate = c("a_R", "b_R", "qR", "z", "eff_SMC", "phi", "size"))
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
#' This function allows the user to update the default priors with new values.
#'
#' @param priors List of priors to update. Default is NULL, which means using the default priors.
#' @param new_priors A named list of new prior specifications to replace the defaults.
#' @return A list containing updated priors.
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
