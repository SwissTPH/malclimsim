library(lubridate)
library(ggplot2)
library(ecmwfr)
library(ncdf4)
library(sf)
library(chirps)
library(dplyr)
library(roxygen2)
library(pak)
library(readxl)
library(odin.dust)
library(scales)  # For better date formatting on x-axis
library(gridExtra)  # For arranging multiple ggplot objects if needed
library(mcstate)

################################################################################
### ---------------- CLIMATE DATA DOWNLOADING AND PROCESSING --------------- ###
################################################################################
library(malclimsim)

# Years used for the analysis
years <- 2014:2023

# Latitude and longitude where climate data (rainfall and temperature) is to be saved
# Rainfall data is from CHIRPS and temperature data is from ERA5
lat <- 8.3
lon <- 17.9
path_to_data <- "C:/Users/putnni/Documents/r-packages/data/"

# Saving the rainfall data locally based on specified path
#save_climate_data(lon, lat, years, path_to_data)

# Reading in climate data saved by `save_climate_data`
temp_path <- paste0(path_to_data, "era5_11051417.nc")
rain_path <- paste0(path_to_data, "chirps_11051418.rds")

# Processing climate data to be used in the model
# `D' controls the number of previous days used for the rainfall rolling average
# calculation and `months_30_days` determines if years are 360 or 365 days
met_360 <- process_climate_data(lon, lat, years, D = 60, temp_path = temp_path,
                            rain_path = rain_path, path_to_data, months_30_days = TRUE)
colnames(met_360)[1] <- "date"
met <- process_climate_data(lon, lat, years, D = 30, temp_path = temp_path,
                            rain_path = rain_path, path_to_data, months_30_days = FALSE)

################################################################################
### ------------------ SMC SCHEDULE FORMATTING FROM DATA ------------------- ###
################################################################################
path_to_SMC <- "C:/Users/putnni/switchdrive/Chad/Data/data-raw/CPS/CPS_coverage_imput_2018.xlsx"

# Format SMC data so that the first column corresponds to the start of the first
# round of SMC and the second column corresponds to the SMC coverage for that round
SMC_data <- readxl::read_excel(path_to_SMC)
View(SMC_data)
SMC_clean <- SMC_data[c("date_start", "smc_couv_tot")]
colnames(SMC_clean) <- c("date_start", "coverage")

SMC_clean <- SMC_clean %>%
  mutate(YearMonth = format(date_start, "%Y-%m")) %>%  # Extract year-month as a new column
  group_by(YearMonth) %>%                        # Group by year-month
  slice(1) %>%                                   # Keep only the first row in each group
  ungroup() %>%                                  # Ungroup to return to a regular data frame
  select(-YearMonth)

smc_schedule <- smc_schedule_from_data(smc_cov = SMC_clean, months_30_days = TRUE, years = years)

################################################################################
### ----------------------- SIMULATING FROM THE MODEL ---------------------- ###
################################################################################
climate_model_1 <- load_model("model_det_1")  # Load the deterministic climate model
climate_model_2 <- load_model("model_det_2")  # Load the deterministic climate model
climate_model_3 <- load_model("model_det_3")  # Load the deterministic climate model
climate_model_4 <- load_model("model_det_4")  # Load the deterministic climate model

# Extract climate vectors
rain <- met_360$anom  # Rainfall anomaly data
temp <- met_360$temp  # Temperature data
# temp <- rep(mean(temp), length(rain))

# Extract key SMC schedule information
SMC <- smc_schedule$SMC  # Indicator for days when an SMC round started (1s and 0s)
decay <- smc_schedule$decay  # Efficacy decay of SMC over time
cov <- smc_schedule$cov  # SMC coverage over time

# Define parameter inputs for the malaria model simulation
param_inputs <- list(
  mu_TS = 1/30, mu_IR = 1/5, eta = 1, mu_RS_C = 1/130, size = 1,
  mu_EI = 1/8, delta_b = 1/(21*365), delta_d = 1/(21*365),
  delta_a = 1/(5 * 365), phi = 1, p_HM = 0.125, p_MH_C = 0.5,
  rho = 1, fT_C = 0.27, z = 1, qR = 0.17, a_R = 0.4, b_R = 3, N = 5e5,
  s = 0.9, p_surv = 0.934, percAdult = 0.81, steps_per_day = 1, SC0 = 0.301,
  EC0 = 0.071, IC0 = 0.042, TC0 = 0.054, RC0 = 0.533, SA0 = 0.281,
  EA0 = 0.067, IA0 = 0.040, TA0 = 0.051, RA0 = 0.562, size_inv = 0,
  eff_SMC = 0, decay = decay, SMC = SMC, cov_SMC = cov,
  c_R_D = rain_rwa, temp = temp_rwa  # Inputs for rainfall and temperature
)


# Run the climate-malaria model simulation
start_date <- ymd("2014-01-01")  # Start date of the simulation
end_date <- ymd("2022-12-31")  # End date of the simulation
results <- data_sim(
  climate_model_4, param_inputs = param_inputs, start_date, end_date,
  month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE
)

################################################################################
### ----------------- PLOTTING SIMULATIONS W/ CLIMATE ---------------------- ###
################################################################################
# Example 1: Plotting only malaria incidence data with only "total" incidence type
plot_time_series(results = results, plot_title = "Malaria Incidence in Rwanda",
                 select_incidence = "<5")

# Example 2: Plotting both malaria incidence (only "<5" and ">=5") and climate data (only "temp")
plot_time_series(results = results, met = met, plot_title = "Malaria Incidence and Climate Data",
                 select_incidence = c("<5", ">=5"), select_climate = "temp", climate_facet = TRUE)

# Example 3: Plotting both malaria incidence (only "total") and climate data (both "temp" and "rollmean")
plot_time_series(results = results, met = met, plot_title = "Malaria Incidence and Climate Data",
                 select_incidence = "total", select_climate = c("temp", "rollmean"), climate_facet = TRUE)

################################################################################
### ------------------------ LOADING OBSERVED DATA ------------------------- ###
################################################################################
path_to_cases <- "C:/Users/putnni/switchdrive/Chad/Data/cases-data/cases_MOISSALA.rds"
obs_cases <- readRDS(path_to_cases)
colnames(obs_cases)[1] <- "date_ymd"
View(cases_moiss)
# First column - `month` - date object corresponding to month of recorded cases
# Second column - `month_no` - numeric where first month is 0
# Third column - `inc_A` - incidence in those >=5
# Fourth column - `inc_C` - incidence in those <5
# Fifth column  - `inc` - total incidence

################################################################################
### --------------------- INFERRING PARAMETER VALUES ----------------------- ###
################################################################################
# Defining parameters to estimate
params_to_estimate <- c(a_R = "a_R", b_R = "b_R", s = "s",
                              qR = "qR", z = "z", eff_SMC = "eff_SMC", "phi",
                              size = "size", p_surv = "p_surv")

params_to_estimate <- c(a_R = "a_R", b_R = "b_R",
                        qR = "qR", z = "z", eff_SMC = "eff_SMC", phi = "phi",
                        size = "size", s = "s")

params_to_estimate <- c(a_R = "a_R", b_R = "b_R",
                        qR = "qR", z = "z", eff_SMC = "eff_SMC", phi = "phi",
                        size = "size", s = "s", c1 = "c1", c3 = "c3", c4 = "c4",
                        c5 = "c5", c6 = "c6", c7 = "c7", c8 = "c8", c9 = "c9",
                        k1 = "k1", k2 = "k2", k3 = "k3", k4 = "k4", k5 = "k5", k7 = "k7")

params_to_estimate <- c(a_R = "a_R", b_R = "b_R",
                        qR = "qR", z = "z", eff_SMC = "eff_SMC", phi = "phi",
                        size = "size", c1 = "c1", k1 = "k1",
                        k2 = "k2", delta_temp = "delta_temp", c4 = "c4",
                        k4 = "k4", k5 = "k5", c5 = "c5", c6 = "c6", k7 = "k7",
                        c7 = "c7", c8 = "c8", p_surv = "p_surv", k3 = "k3", c3 = "c3",
                        s = "s")

params_to_estimate <- c(a_R = "a_R", b_R = "b_R",
                        qR = "qR", z = "z", eff_SMC = "eff_SMC", phi = "phi",
                        size = "size", c6 = "c6", c7 = "c7", c8 = "c8",
                        p_surv = "p_surv", s = "s", k1 = "k1", c1 = "c1")

################################################################################
### ----------------------- DEFINING MCMC PARAMETERS ----------------------- ###
################################################################################
# Call function with default values
params_default <- create_mcmc_params(n_steps = 1000)

# Call function with custom values for n_steps and n_chains
params_custom <- create_mcmc_params(n_steps = 20000, n_chains = 6)

# Access the parameters
adaptive_params_1 <- params_default$adaptive_params
control_params_1 <- params_default$control_params

adaptive_params_2 <- params_custom$adaptive_params
control_params_2 <- params_custom$control_params

################################################################################
### ---------------- DEFINING PROPOSAL AND STARTING VALUES ----------------- ###
################################################################################
# defining starting values for chains
start_values <- create_start_values(params_to_estimate, params_default$control_params,
                                    model = climate_model_4, param_inputs = param_inputs, random = TRUE, seed = NULL)

# choose proposal matrix
proposal_matrix <- create_proposal_matrix(params_to_estimate = params_to_estimate,
                                          model = climate_model_4, param_inputs = param_inputs)


# Default behavior with predefined param_names, min_max_start_values, and proposal_variance
# start_values <- create_start_values(
#  params_to_estimate = c("phi", "qR", "mu_IR"),
#  control_params = list(n_chains = 4)
# )

################################################################################
### -------------- EXAMINING AND CHANGING PRIOR DISTRIBUTIONS -------------- ###
################################################################################
viewed_priors <- view_priors(param_inputs, proposal_matrix, params_to_estimate)
plot_priors(param_inputs, proposal_matrix, params_to_estimate)
# Update default priors
updated_priors <- update_priors(param_inputs, proposal_matrix, params_to_estimate,
                                new_priors = list(
  phi = list(initial = 0.3, min = 0.02, max = 0.9),
  qR = list(prior = function(p) dbeta(p, 5, 10, log = TRUE),)
))


################################################################################
### --------------- RUNNING MCMC (Random Walk Metropolis) ------------------ ###
################################################################################
dates_for_inf <- c("2014-01-01", "2022-12-31") # this isn't really doing anything
results <- inf_run(model = climate_model_4, param_inputs = param_inputs,
                           control_params = params_default$control_params,
                           params_to_estimate = params_to_estimate,
                           proposal_matrix = proposal_matrix,
                           adaptive_params = params_default$adaptive_params,
                           start_values = start_values, month = TRUE,
                           dates = dates_for_inf, age_for_inf = 'u5',
                           synthetic = FALSE, incidence_df = obs_cases,
                           save_trajectories = FALSE, noise = FALSE,
                           rerun_n = 1000, rerun_random = TRUE,
                   param_priors = NULL)

################################################################################
### ----------------------- SOME MCMC DIAGNOSTICS ------------------------- ####
################################################################################
MCMC_diag(results)

################################################################################
### -------------------- EVALUATING MODEL ERROR --------------------------- ####
################################################################################
# Maximum likelihood/posterior
max_ll_post(results)

# Posterior predictive check
plot_observed_vs_simulated(results = results, obs_cases,
                           start_date = "2014-01-01", end_date = "2022-12-31",
                           model = climate_model_4, add_ribbon = TRUE, n_samples = 5,
                           groups = c("inc_C"), met = met_360, climate_facet = TRUE)


# Extract the parameters that had the highest log posterior from your MCMC results
max_posterior_params <- extract_max_posterior_params(results)

# Update the parameter list with the maximum posterior parameters
param_inputs_updated <- update_param_list(results$param_inputs, max_posterior_params)

# Run a model simulation using the updated parameters
start_date <- "2014-01-01"  # Replace with your actual start date
end_date <- "2022-12-31"    # Replace with your actual end date
simulated_df <- simulate_with_max_posterior_params(results, start_date, end_date, model = climate_model_1)

error_metrics <- calculate_model_errors(obs_cases, simulated_df, date_column = "date_ymd")
residual_plot <- plot_residuals(obs_cases, simulated_df, date_column = "date_ymd", groups = c("inc_A", "inc_C"))
residual_plot

################################################################################
### --------- ESTIMATING CASES AVERTED (THIS IS FINALLY INFERENCE) -------- ####
################################################################################
# Step 1: Simulate the scenarios
current_smc_df <- simulate_with_smc(results, start_date = "2014-01-01", end_date = "2022-12-31", model = climate_model)
no_smc_df <- simulate_without_smc(results, start_date = "2014-01-01", end_date = "2022-12-31", model = climate_model)
full_smc_df <- simulate_with_full_coverage(results, start_date = "2014-01-01", end_date = "2022-12-31", model = climate_model, coverage_level = 1)

# Step 2: Calculate the number and percent of cases averted
cases_averted_full_coverage <- calculate_cases_averted(full_smc_df, no_smc_df)
print(cases_averted_full_coverage)

# Step 3: Optionally, compare full coverage to current SMC coverage
cases_averted_current_vs_full <- calculate_cases_averted(full_smc_df, current_smc_df)
print(cases_averted_current_vs_full)

# Assuming the scenarios have been simulated and stored as `observed_df`, `no_smc_df`, `current_smc_df`, `full_smc_df`
plot_smc_scenarios(
  observed_df = obs_cases,
  no_smc_df = no_smc_df,
  current_smc_df = current_smc_df,
  full_smc_df = full_smc_df,
  groups = c("inc_C"),
  include_observed = FALSE
)
