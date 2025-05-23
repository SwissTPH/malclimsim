% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_distributions.R
\name{view_priors}
\alias{view_priors}
\title{View Default Priors}
\usage{
view_priors(param_inputs, proposal_matrix, params_to_estimate, priors = NULL)
}
\arguments{
\item{params_to_estimate}{Character vector specifying which parameters' priors should be viewed.}

\item{priors}{List of priors to view. If NULL, the default priors will be displayed.}
}
\value{
A data frame containing the details of each prior specified in params_to_estimate.
}
\description{
This function displays the default priors and their details, including initial values, min/max bounds, and prior distributions.
}
\examples{
view_priors(param_inputs, proposal_matrix, params_to_estimate = c("a_R", "b_R", "qR", "z", "eff_SMC", "phi", "size"))
}
