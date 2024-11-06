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
save_climate_data(lon, lat, years, path_to_data)

# Reading in climate data saved by `save_climate_data`
temp_path <- paste0(path_to_data, "era5_11051417.nc")
rain_path <- paste0(path_to_data, "chirps_11051418.rds")

# Processing climate data to be used in the model
# `D' controls the number of previous days used for the rainfall rolling average
# calculation and `months_30_days` determines if years are 360 or 365 days
met <- process_climate_data(lon, lat, years, D = 30, temp_path = temp_path,
                            rain_path = rain_path, path_to_data, months_30_days = TRUE)

################################################################################
### ------------------ SMC SCHEDULE FORMATTING FROM DATA ------------------- ###
################################################################################
path_to_SMC <- "C:/Users/putnni/switchdrive/Chad/Data/data-raw/CPS/CPS_coverage_imput_2018.xlsx"

# Format SMC data so that the first column corresponds to the start of the first
# round of SMC and the second column corresponds to the SMC coverage for that round
SMC_data <- read_excel(path_to_SMC)
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

load_model("model_det_1.R")

climate_model <- odin.dust::odin_dust("climate_model_deterministic.R")  # Load the deterministic climate model


# Extract climate vectors
rain <- met$anom  # Rainfall anomaly data
temp <- met$temp  # Temperature data

# Extract key SMC schedule information
SMC <- smc_schedule$SMC  # Indicator for days when an SMC round started (1s and 0s)
decay <- smc_schedule$decay  # Efficacy decay of SMC over time
cov <- smc_schedule$cov  # SMC coverage over time

param_inputs <- list(
  mu_TS = 1/30, mu_IR = 1/5, eta = 1, mu_RS_C = 1/130, size = 1,
  mu_EI = 1/8, delta_b = 1/(21*365), delta_d = 1/(21*365),
  delta_a = 1/(5 * 365), phi = 1, p_HM = 0.125, p_MH_C = 0.5,
  rho = 1, fT_C = 0.27, z = 1, qR = 0.17, a_R = 0.4, b_R = 3, N = 5e5,
  s = 0.9, p_surv = 0.934, percAdult = 0.81, steps_per_day = 1, SC0 = 0.301,
  EC0 = 0.071, IC0 = 0.042, TC0 = 0.054, RC0 = 0.533, SA0 = 0.281,
  EA0 = 0.067, IA0 = 0.040, TA0 = 0.051, RA0 = 0.562, size_inv = 0,
  eff_SMC = 0, decay = decay, SMC = SMC, cov_SMC = cov,
  c_R_D = rain, temp = temp  # Inputs for rainfall and temperature
)


