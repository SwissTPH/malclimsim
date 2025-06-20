## code to prepare `DATASET` dataset goes here
library(malclimsim)
library(stringr)
out_dir <- "C:/Users/putnni/Documents/git-hub-repositories/vignette-storage/malclimsim/data/"

path_to_data <- "C:/Users/putnni/switchdrive/Chad/Data" # Replace this with location of chad switchdrive folder
########################
# Simulating SMC Data  #
########################
path_to_SMC <- file.path(path_to_data,"data-raw/CPS/CPS_coverage_imput_2018_with_2023.xlsx")
smc_data_orig <- readxl::read_excel(path_to_SMC)

smc_data_raw <- smc_data_orig

# SMC starts on the 15th of each month
smc_data_raw$date_start <- as.Date(str_replace(smc_data_orig$date_start, "..$", "15"))
smc_data_raw$smc_couv_card <- 0.60
smc_data_raw$smc_couv_card_lower <- 0.55
smc_data_raw$smc_couv_card_upper <- 0.65
smc_data_raw$smc_couv_tot <- 0.85
smc_data_raw$smc_couv_tot_lower <- 0.80
smc_data_raw$smc_couv_tot_upper <- 0.90

usethis::use_data(smc_data_raw, overwrite = TRUE)

smc_data <- load_clean_smc_data(data = smc_data_raw)
smc_schedule <- smc_schedule_from_data(
  smc_cov = smc_data,
  months_30_days = FALSE,
  years = 2018:2023
)

#####################################
# Simulating Weekly Malaria Cases   #
#####################################
malaria_model <- load_model("model_new_R_with_FOI")  # Load the deterministic climate model

mcmc_results <- readRDS("C:/Users/putnni/Documents/git-hub-repositories/vignette-storage/malclimsim/mcmc-results/mcmc_results_stage_3.rds")

# Extract MAP (max posterior) parameter set
max_posterior_params <- extract_max_posterior_params(mcmc_results)
param_inputs_best <- update_param_list(mcmc_results$param_inputs, max_posterior_params)
param_inputs_best$eff_SMC <- 0.7
param_inputs_best$s <- 10
param_inputs_best$qR <- 0.05
param_inputs_best$alpha <- 0.15
param_inputs_best$SMC <- smc_schedule$SMC
param_inputs_best$cov_SMC <- smc_schedule$cov
param_inputs_best$decay <- smc_schedule$decay
param_inputs_best$N <- 150000

# Create covariate matrix (population growth)
r_df <- get_population_scaling(
  n             = nrow(mcmc_results$incidence_df),
  month         = FALSE,
  growth_rate_C = 1.000071,
  growth_rate_A = 1.000092
)
pop_growth_df <- data.frame(
  date_ymd = mcmc_results$incidence_df$date_ymd,
  r_C = r_df$r_C
)

# Define transformation for SMC-adjusted incidence
mu_transform_C <- function(inc_df, param_inputs) {
  inc_df$inc_C * inc_df$r_C
}

# Simulate cases using the best parameter set (deterministic, no noise)
sim_best <- data_sim(
  model         = malaria_model,
  param_inputs  = param_inputs_best,
  start_date    = "2018-01-01",
  end_date      = "2023-12-31",
  noise         = FALSE,
  mu_transform_C = mu_transform_C,
  mu_transform_A = NULL,
  covariate_matrix = pop_growth_df,
  month = FALSE,
  save = FALSE
)

sim_best$inc_C <- sim_best$inc_C_transformed
sim_best <- sim_best %>% select(!c(inc, r_C, inc_C_transformed))

set.seed(20)
sim_best$inc_A <- rnbinom(nrow(sim_best), size = 10, mu = sim_best$inc_A)
sim_best$inc_C <- rnbinom(nrow(sim_best), size = 10, mu = sim_best$inc_C)
sim_best$week_no <- 0:(nrow(sim_best) - 1)

plot(sim_best$date_ymd, sim_best$inc_C, type = "l")

obs_cases <- sim_best

usethis::use_data(obs_cases, overwrite = TRUE)
