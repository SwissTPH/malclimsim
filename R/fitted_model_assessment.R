#' Extract and Print Maximum Log-Likelihood and Log-Posterior
#'
#' This function extracts and optionally prints the maximum values of the log-likelihood
#' and log-posterior from the results of an MCMC run. It is useful for quickly evaluating
#' the fit of a model.
#'
#' @param results A list containing the MCMC results from an inference run,
#' including the `coda_pars` matrix with columns for `log_likelihood` and `log_posterior`.
#' @param print_ll Boolean. If `TRUE`, the maximum log-likelihood and log-posterior values
#' will be printed to the console (default: `TRUE`).
#'
#' @return A list with two named elements:
#' \describe{
#'   \item{max_log_likelihood}{The maximum value of the log-likelihood.}
#'   \item{max_log_posterior}{The maximum value of the log-posterior.}
#' }
#'
#' @examples
#' # Simulate results with coda_pars
#' results <- list(
#'   coda_pars = data.frame(
#'     log_likelihood = rnorm(1000, mean = -100, sd = 10),
#'     log_posterior = rnorm(1000, mean = -110, sd = 12)
#'   )
#' )
#'
#' # Extract and print maximum log-likelihood and log-posterior
#' max_ll_post(results)
#'
#' # Suppress printing and only return the values
#' max_ll_post(results, print_ll = FALSE)
#'
#' @export
max_ll_post <- function(results, print_ll = TRUE) {
  # Check if 'coda_pars' exists and contains required columns
  if (!"coda_pars" %in% names(results)) {
    stop("The results list must contain 'coda_pars'.")
  }

  if (!all(c("log_likelihood", "log_posterior") %in% colnames(results$coda_pars))) {
    stop("'coda_pars' must contain 'log_likelihood' and 'log_posterior' columns.")
  }

  # Extract MCMC posterior results
  coda_pars <- results$coda_pars

  # Find maximum values for log_likelihood and log_posterior
  max_log_likelihood <- max(coda_pars[ , "log_likelihood"], na.rm = TRUE)
  max_log_posterior <- max(coda_pars[ , "log_posterior"], na.rm = TRUE)

  # Print the results
  if(print_ll){
    cat("Maximum Log-Likelihood:", max_log_likelihood, "\n")
    cat("Maximum Log-Posterior:", max_log_posterior, "\n")
  }
  return(list(max_log_likelihood = max_log_likelihood,
              max_log_posterior = max_log_posterior))
}
