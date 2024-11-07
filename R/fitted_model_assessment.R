#' Print Maximum Log-Likelihood and Log-Posterior
#'
#' This function extracts and prints the maximum values of the log-likelihood and
#' log-posterior from the results of an MCMC run. Useful for quick evaluation of model fit.
#'
#' @param results A list containing the MCMC results from an inference run,
#' including the `coda_pars` matrix with columns for `log_likelihood` and `log_posterior`.
#'
#' @return NULL. The function prints the maximum values of log-likelihood and log-posterior.
#' @export
#'
#' @examples
#' # Assuming `results` is the result of an inference run with MCMC sampling:
#' max_ll_post(results)
max_ll_post <- function(results) {
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
  cat("Maximum Log-Likelihood:", max_log_likelihood, "\n")
  cat("Maximum Log-Posterior:", max_log_posterior, "\n")
}
