# Define a function to initialize priors based on parameters in the model
initialize_priors <- function(param_inputs, proposal_matrix, params_to_estimate) {

  # Step 1: Extract parameter names from the proposal matrix to ensure we match the model
  param_names <- rownames(proposal_matrix)

  # Step 2: Get the valid parameters (those that exist in both param_inputs and proposal_matrix)
  valid_params_in_model <- names(param_inputs)[names(param_inputs) %in% param_names]

  # Step 3: Ensure we include parameters from params_to_estimate
  valid_params <- unique(c(valid_params_in_model, params_to_estimate))

  # Step 4: Filter out parameters that are vectors (i.e., those whose length is greater than 1)
  valid_params <- valid_params[sapply(param_inputs[valid_params], length) == 1]

  # Define the prior functions and parameter bounds as specified in the model
  default_priors <- list(

    # 1. Transmission and Recovery Parameters
    phi = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 40, 12, log = TRUE)),
    qR = list(initial = 0.2, min = 1e-6, max = 1, prior = function(p) dbeta(p, 1, 150, log = TRUE)),
    mu_RS_C = list(initial = 1/200, min = 1/400, max = 1/120, prior = function(p) dgamma(p, shape = 2, rate = 480, log = TRUE)),
    mu_EI = list(initial = 1/8, min = 1/15, max = 1/6, prior = function(p) dgamma(p, shape = 6, rate = 63, log = TRUE)),
    phi_C2 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 40, 12, log = TRUE)),
    phi_A = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 40, 12, log = TRUE)),
    mu_TS = list(initial = 0.2, min = 0, max = 2, prior = function(p) dunif(p, min = 0, max = 2, log = TRUE)),
    mu_IR = list(initial = 1/5, min = 0.001, max = 1, prior = function(p) dunif(p, min = 0.001, max = 1, log = TRUE)),

    # 2. Survival and Population Parameters
    a_R = list(initial = 0.5, min = 0.01, max = 1, prior = function(p) dbeta(p, 6, 12, log = TRUE)),
    b_R = list(initial = 3, min = -Inf, max = Inf, prior = function(p) dnorm(p, mean = 2, sd = 0.25, log = TRUE)),
    p_surv = list(initial = 0.91, min = 0.88, max = 0.97, prior = function(p) dnorm(p, mean = 0.91, sd = 0.01, log = TRUE)),
    size = list(initial = 5.5, min = 5, max = 10, prior = function(p) dnorm(p, mean = 5.5, sd = 0.5, log = TRUE)),

    # 3. Intervention and Effectiveness Parameters
    eff_SMC = list(initial = 0.6, min = 0.00, max = 1, prior = function(p) dbeta(p, 2, 2, log = TRUE)),
    z = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 4, 2, log = TRUE)),
    p_MH_C = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dunif(p, min = 0.01, max = 1, log = TRUE)),
    rho = list(initial = 0.2, min = 0.01, max = 0.99, prior = function(p) dunif(p, min = 0.01, max = 0.99, log = TRUE)),

    # 4. Population Proportions and Initial Conditions
    s = list(initial = 0.8, min = 0.01, max = 1, prior = function(p) dgamma(p, shape = 20, rate = 20, log = TRUE)),
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
    delta_temp = list(initial = 0, min = -5, max = 2, prior = function(p) dunif(p, min = -5, max = 2, log = TRUE)),

    # 6. Adjustment and Scaling Parameters
    tau = list(initial = 0.2, min = 0, max = 1, prior = function(p) dbeta(p, 40, 12, log = TRUE)),
    kappa = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 0.125, 0.125, log = TRUE)),
    z_A = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 0.125, 0.125, log = TRUE)),
    z_C2 = list(initial = 0.2, min = 0.01, max = 1, prior = function(p) dbeta(p, 0.125, 0.125, log = TRUE))
  )



  # Initialize an empty list to store priors for each parameter in the proposal matrix
  priors <- list()

  # Loop through each parameter in the proposal matrix and assign priors
  for (param in valid_params) {
    # Check if a specific prior is defined for this parameter in default_priors
    if (param %in% names(default_priors)) {
      param_info <- default_priors[[param]]
      prior_function <- param_info$prior
      initial_value <- param_info$initial
      min_value <- param_info$min
      max_value <- param_info$max
    } else {
      # Apply a generic default prior if none specified
      prior_function <- function(p) dunif(p, min = 0, max = 1, log = TRUE)
      initial_value <- 0.2
      min_value <- 0.01
      max_value <- 1
    }

    # Initialize the parameter with the appropriate prior, setting bounds and initial values
    priors[[param]] <- mcstate::pmcmc_parameter(
      name = param,
      initial = initial_value,
      min = min_value,
      max = max_value,
      prior = prior_function
    )
  }

  # Return the list of initialized priors
  return(priors)
}


