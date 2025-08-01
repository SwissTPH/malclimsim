% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference_aux_functions.R
\name{update_inf_stage}
\alias{update_inf_stage}
\title{Update proposal matrix, start values, and MCMC settings for next stage}
\usage{
update_inf_stage(
  results_obj,
  proposal_matrix,
  param_names,
  S_prev = 3000,
  draw_n = 4,
  shrink = 0.8,
  stage,
  n_steps = NULL
)
}
\arguments{
\item{results_obj}{Result object returned by inf_run (with MCMC traces)}

\item{proposal_matrix}{Current proposal variance matrix}

\item{param_names}{Character vector of parameter names (rownames of proposal_matrix)}

\item{S_prev}{Integer: number of samples to use for variance-covariance extraction}

\item{draw_n}{Integer: number of random draws for new start values per chain}

\item{shrink}{Numeric: shrinkage factor for bounds (e.g. 0.8 retains 80\% of interval)}

\item{stage}{Character: which stage to pull from create_mcmc_params()}

\item{n_steps}{Integer or NULL: if not NULL, overrides control_params$n_steps}
}
\value{
A list with elements:
\itemize{
\item proposal_matrix : updated proposal variance matrix
\item start_values    : matrix of new start values (chains x parameters)
\item adaptive_params : list of adaptive settings for this stage
\item control_params  : list of control settings for this stage (with n_steps possibly overridden)
}
}
\description{
Update proposal matrix, start values, and MCMC settings for next stage
}
