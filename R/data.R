#' Observed Malaria Incidence Used for Model Calibration
#'
#' Weekly simulated malaria incidence data used to represent observed cases for calibration
#' and validation of the transmission model. Generated deterministically using maximum a posteriori (MAP) parameter values.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{date_ymd}{Date corresponding to each week (class: Date)}
#'   \item{inc_A}{Simulated incidence in adults}
#'   \item{inc_C}{Simulated incidence in children}
#'   \item{week_no}{Week index relative to simulation start}
#' }
#'
#' @source Generated using \code{data_sim()} with MAP-calibrated model parameters.
"obs_cases"


#' Raw SMC Coverage Data from CPS
#'
#' Raw seasonal malaria chemoprevention (SMC) campaign data from a cleaned CPS Excel file covering 2018 to 2023.
#' Includes total and card-confirmed coverage with associated uncertainty bounds.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{date_start}{Date the SMC round began (class: Date)}
#'   \item{epiweek_start}{Epidemiological week the round began}
#'   \item{smc_couv_card}{Card-confirmed SMC coverage}
#'   \item{smc_couv_card_lower}{Lower bound of card-confirmed coverage}
#'   \item{smc_couv_card_upper}{Upper bound of card-confirmed coverage}
#'   \item{smc_couv_tot}{Total SMC coverage}
#'   \item{smc_couv_tot_lower}{Lower bound of total coverage}
#'   \item{smc_couv_tot_upper}{Upper bound of total coverage}
#'   \item{smc_round}{SMC round number (1 through 4 or 5)}
#' }
#'
#' @source Extracted from a cleaned CPS Excel dataset for 2018 to 2023.
"smc_data_raw"
