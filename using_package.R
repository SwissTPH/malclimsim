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

################################################################################
### ---------------- CLIMATE DATA DOWNLOADING AND PROCESSING --------------- ###
################################################################################
library(malclimsim)
years <- 2014:2023
lat <- 8.3
lon <- 17.9
path_to_data <- "C:/Users/putnni/Documents/r-packages/data/"

save_climate_data(lon, lat, years, path_to_data)

temp_path <- paste0(path_to_data, "era5_11051417.nc")
rain_path <- paste0(path_to_data, "chirps_11051418.rds")

met <- process_climate_data(lon, lat, years, D = 30, temp_path = temp_path,
                            rain_path = rain_path, path_to_data)

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

smc_cov <- SMC_clean

smc_df <- smc_schedule_from_data(smc_cov = SMC_clean, months_30_days = TRUE, years = years)
