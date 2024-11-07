# age_for_inf can be either 'total', 'sep_ages', or 'total + sep_ages'
inf_run <- function(model, param_inputs,
                    control_params, params_to_estimate, proposal_matrix,
                    adaptive_params, start_values, noise = FALSE, seed = 24, month = FALSE,
                    month_unequal_days = FALSE, dates,
                    age_for_inf, synthetic = TRUE, incidence_df = NULL, clim_model = FALSE,
                    save_trajectories = TRUE, rerun_n = Inf, rerun_random = FALSE, dir = ""){

  ##############################################################################
  #### ------------------ GENERATING SYNTHETIC DATA ------------------------ ###
  ##############################################################################

  if(synthetic){
    # incidence dataframe in specific format from data_sim
    if(month){
      incidence_df <- data_sim(model, param_inputs,
                               start_date = dates[1], end_date = dates[2],
                               month = TRUE, save = FALSE,
                               month_unequal_days = month_unequal_days)
    }else{incidence_df <- data_sim(model, param_inputs, start_date = dates[1],
                                   end_date = dates[2], month = FALSE, save = FALSE,
                                   month_unequal_days = month_unequal_days)}

    if(noise){
      set.seed(seed)
      incidence_df$inc <- rpois(nrow(incidence_df), lambda = incidence_df$inc)
      incidence_df$inc_C <- rpois(nrow(incidence_df), lambda = incidence_df$inc_C)
      incidence_df$inc_A <- rpois(nrow(incidence_df), lambda = incidence_df$inc_A)
    }
  }

  # Loading in SMC and decay vectors
  SMC <- param_inputs$SMC
  decay <- param_inputs$decay
  cov_SMC <- param_inputs$cov_SMC

  temp <- param_inputs$temp
  c_R_D <- param_inputs$c_R_D

  ##############################################################################
  #### -------------- DEFINING COMPARISON FUNCTIONS ------------------------ ###
  ##############################################################################
  # Negative binomial distribution function - see 'observation_functions' for more details

  # Remove the date column from the observed incidence data
  incidence_observed <- incidence_df[-1]  # Removing the first column (date)
  comparison_fn <- generate_incidence_comparison(month, age_for_inf, incidence_observed)


  simulated_result <- data_sim_for_inference(model, param_inputs = param_inputs, dates = dates,
                                             noise = FALSE, month = month)
  initial_time_obs <- simulated_result$month_no[simulated_result$date_ymd == incidence_df$month[1]] # matching first observation to simulation
  if(length(initial_time_obs) > 0){
    incidence_observed$month_no <- seq(initial_time_obs, (nrow(incidence_observed) + initial_time_obs - 1))
  }else{
    initial_time_obs <- 0
  }

  ##############################################################################
  #### -------------- RELATING TIME STEPS AND OBSERVED DATA ---------------- ###
  ##############################################################################
  # Assuming filter_data is now a function in your package or environment
  filt_data <- filter_data(incidence_observed, month = month, initial_time_obs = initial_time_obs)

  ##############################################################################
  #### --------------- DEFINING STATES SAVED FOR SIMULATION ---------------- ###
  ##############################################################################
  # Assuming 'index' is now a function in your package or environment
  filter <- mcstate::particle_deterministic$new(data = filt_data, model = model,
                                                index = index, compare = comparison_fn)
  filter$run(c(param_inputs))

  ##############################################################################
  ## - Function necessary for all input vectors to be treated as parameters - ##
  ##############################################################################
  # Assuming make_transform is now a function in your package or environment
  transform_fn <- make_transform(temp = temp, c_R_D = c_R_D, SMC = SMC, decay = decay, cov_SMC = cov_SMC)

  ##############################################################################
  #### ----------------- DEFINING PRIOR DISTRIBUTIONS ---------------------- ###
  ##############################################################################
  # Initialize the parameters by calling the function
  param_priors <- initialize_priors(param_inputs, proposal_matrix, params_to_estimate)

  # choose proposal matrix - default MVN with given covariance structure
  mcmc_pars <- mcstate::pmcmc_parameters$new(param_priors,
                                             proposal_matrix,
                                             make_transform(temp = temp, c_R_D = c_R_D,
                                                            SMC = SMC, decay = decay,
                                                            cov_SMC = cov_SMC))

  ##############################################################################
  #### ----------------- FIXING PARAMETERS NOT TO BE INFERRED -------------- ###
  ##############################################################################
  # Assuming param_inputs is a named list that contains the values for the parameters
  # and params_to_estimate is the list of parameters you want to estimate

  # First, filter param_inputs to only include parameters found in param_priors
  valid_params <- names(param_inputs)[names(param_inputs) %in% names(param_priors)]

  # Now filter param_inputs to only keep the valid parameters
  paramFix <- param_inputs[valid_params]

  # Further, exclude parameters to be estimated (those in params_to_estimate)
  paramFix <- paramFix[setdiff(names(paramFix), params_to_estimate)]

  # Filter out any parameters that are vectors (length > 1), keeping only scalars (length == 1)
  paramFix <- paramFix[sapply(paramFix, function(x) length(x) == 1)]

  # Convert to a named numeric vector (removes the list structure)
  paramFix <- unlist(paramFix)

  mcmc_pars <- mcmc_pars$fix(paramFix)

  ##############################################################################
  #### ----------------- DEFINING PARAMETERS OF MCMC RUN ------------------- ###
  ##############################################################################
  n_chains = control_params$n_chains
  n_burnin = control_params$n_burnin
  n_steps = control_params$n_steps
  n_workers = control_params$n_workers
  n_threads_total = control_params$n_threads_total

  start_values <- reorder_start_values(start_values, param_priors)

  control1 <- mcstate::pmcmc_control(n_steps = 10, n_burnin = 0, progress = TRUE,
                                     n_chains = n_chains)

  control2 <- mcstate::pmcmc_control(n_steps = n_steps, n_burnin = n_burnin, progress = TRUE, # for some reason this doesn't recognize temp???
                                     n_chains = n_chains, n_workers = n_workers,
                                     n_threads_total = n_threads_total,
                                     adaptive_proposal = adaptive_params,
                                     save_trajectories = save_trajectories,
                                     rerun_every = rerun_n, rerun_random = rerun_random)

  ##############################################################################
  #### ------------------------ RUNNING MCMC ------------------------------- ###
  ##############################################################################
  mcmc_run <- mcstate::pmcmc(mcmc_pars, filter, initial = start_values, control = control1) # I have to do this first for some weird reason
  mcmc_run <- mcstate::pmcmc(mcmc_pars, filter, initial = start_values, control = control2)
  coda_pars <- as.mcmc(cbind(mcmc_run$probabilities, mcmc_run$pars))
  results <- list(mcmc_run = mcmc_run, coda_pars = coda_pars, paramFix = paramFix,
                  param_inputs = param_inputs, incidence_df = incidence_df,
                  model = model, paramList = paramList,
                  n_chains = n_chains)
  return(results)
}
