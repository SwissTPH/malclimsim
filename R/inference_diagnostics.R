#' Convert MCMC Output to List of Chains for Visualization
#'
#' Converts the output of a multi-chain MCMC run into a `coda::mcmc.list` object, separating chains
#' for easier visualization and diagnostics (e.g., trace plots or Gelman-Rubin diagnostics).
#'
#' @param mcmc_run An object returned by `mcstate::pmcmc()`, containing posterior samples and a `chain` column indicating chain IDs.
#'
#' @return An `mcmc.list` object (from the `coda` package), where each element is an individual chain converted to an `mcmc` object.
#'
#' @details
#' This function assumes that `mcmc_run$pars` contains the posterior samples and that there is a
#' `chain` vector in `mcmc_run$chain` specifying the chain each row belongs to.
#' @export
plot_chains <- function(mcmc_run) {
  chains <- mcmc_run$chain
  df <- data.frame(as.matrix(mcmc_run$pars), chains) # Combine parameter values with chain labels
  chains_sep <- split(df, df$chains) # Separate data by chain
  chains_list <- coda::mcmc.list(lapply(chains_sep, function(x) coda::as.mcmc(x[-ncol(x)]))) # Convert to 'mcmc'
  return(chains_list)
}

#' Plot Correlation of MCMC Samples
#'
#' This function creates a correlation plot for MCMC samples, showing the relationships between sampled parameters.
#' The upper panels display correlation values, and the lower panels show smoothed plots of parameter relationships.
#'
#' @param results list: Results from the `inf_run` function, containing MCMC samples.
#' @param title character: Optional title for the plot.
#'
#' @return A `GGally` ggpairs object visualizing the correlation structure of MCMC samples.
#' @export
plot_corr <- function(results, title = NULL) {
  suppressWarnings(suppressMessages({
    # Convert the relevant portion of results to a dataframe and sample 1000 rows
    df <- as.data.frame(results[[2]][, -(1:3)])
    sampled_df <- df[sample(1:nrow(df), 1000), ]

    # Handle single-parameter cases
    if(is.null(dim(sampled_df))) {
      sampled_df <- data.frame(eff_SMC = sampled_df)
    }

    # Create a correlation plot with upper panels showing correlation values and lower panels showing smoothed plots
    p <- GGally::ggpairs(
      sampled_df,
      upper = list(continuous = GGally::wrap("cor", size = 3)),
      lower = list(continuous = GGally::wrap("smooth", alpha = 0.3, size = 0.1))
    )

    # Apply additional styling to make plot more readable
    p <- p + theme_bw() +
      theme(strip.text = element_text(size = 6),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(angle = 0, hjust = 1),
            axis.ticks.length = unit(0.25, "cm")
      ) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 2))

    # Add title if provided
    if (!is.null(title)) {
      p <- p + ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    }
  }))

  return(p)
}

#' MCMC Diagnostic and Summary Plots with Save Options
#'
#' @param results MCMC result object.
#' @param params Diagnostics to compute.
#' @param thin Thinning interval.
#' @param save Logical; whether to save output (default: FALSE).
#' @param output_dir Directory to save outputs (default: "mcmc_diagnostics").
#' @param file_prefix Optional prefix for filenames.
#' @param plot_width Width of plots in inches (default: 7).
#' @param plot_height Height of plots in inches (default: 5).
#'
#' @return Named list of diagnostics.
#'
#' @importFrom xtable xtable
#' @export
MCMC_diag <- function(results,
                      params = c("trace", "gelman", "corr", "ess", "acf", "quantiles", "acceptance"),
                      thin = 1,
                      save = FALSE,
                      output_dir = "mcmc_diagnostics",
                      file_prefix = "mcmc",
                      plot_width = 7,
                      plot_height = 5) {

  suppressWarnings(suppressMessages({
    if (save && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    pars <- results[[1]]$pars
    n_chains <- results$n_chains
    n_samples <- nrow(pars) / n_chains

    chain_list <- split(pars, rep(1:n_chains, each = n_samples))
    chains <- lapply(chain_list, function(x) {
      as.mcmc(matrix(x, nrow = n_samples, ncol = ncol(pars), dimnames = list(NULL, colnames(pars))))
    })
    mcmc_chains <- as.mcmc.list(chains)

    if (thin > 1) {
      mcmc_chains <- lapply(mcmc_chains, function(chain) window(chain, thin = thin))
      mcmc_chains <- as.mcmc.list(mcmc_chains)
    }

    results_list <- list()

    # Trace plot
    if ("trace" %in% params) {
      cat("\n TRACE PLOT \n")
      if (save) {
        pdf(file.path(output_dir, paste0(file_prefix, "_trace.pdf")),
            width = plot_width, height = plot_height)
        plot(mcmc_chains)
        dev.off()
      } else {
        plot(mcmc_chains)
      }
      results_list$trace <- mcmc_chains
    }

    # Gelman-Rubin
    if ("gelman" %in% params && n_chains > 1) {
      cat("\n GELMAN-RUBIN STATISTIC \n")
      gelman_result <- gelman.diag(mcmc_chains)
      print(gelman_result)
      results_list$gelman <- gelman_result
      if (save) {
        tex <- xtable::xtable(gelman_result$psrf, caption = "Gelman-Rubin Diagnostic")
        print(tex, file = file.path(output_dir, paste0(file_prefix, "_gelman.tex")))
      }
    }

    # Correlation plot
    if ("corr" %in% params) {
      cat("\n CORRELATION PLOT \n")
      results_list$corr <- "Correlation plot generated"
      if (save) {
        pdf(file.path(output_dir, paste0(file_prefix, "_corr.pdf")),
            width = plot_width, height = plot_height)
        plot_corr(results)
        dev.off()
      } else {
        plot_corr(results)
      }
    }

    # Effective Sample Size
    if ("ess" %in% params) {
      cat("\n EFFECTIVE SAMPLE SIZE \n")
      ess_result <- effectiveSize(mcmc_chains)
      print(ess_result)
      results_list$ess <- ess_result
      if (save) {
        tex <- xtable::xtable(as.matrix(ess_result), caption = "Effective Sample Sizes")
        print(tex, file = file.path(output_dir, paste0(file_prefix, "_ess.tex")))
      }
    }

    # Autocorrelation
    if ("acf" %in% params) {
      cat("\n AUTOCORRELATION FUNCTION \n")
      acf_matrix <- as.matrix(do.call(rbind, mcmc_chains))
      if (save) {
        pdf(file.path(output_dir, paste0(file_prefix, "_acf.pdf")),
            width = plot_width, height = plot_height)
        acf(acf_matrix, main = "Autocorrelation")
        dev.off()
      } else {
        acf(acf_matrix, main = "Autocorrelation")
      }
      results_list$acf <- "ACF plotted"
    }

    # Posterior Quantiles
    if ("quantiles" %in% params) {
      cat("\n POSTERIOR QUANTILES \n")
      quantiles <- apply(as.matrix(do.call(rbind, mcmc_chains)), 2, quantile,
                         probs = c(0.025, 0.5, 0.975))
      print(quantiles)
      results_list$quantiles <- quantiles
      if (save) {
        tex <- xtable::xtable(quantiles, caption = "Posterior Quantiles")
        print(tex, file = file.path(output_dir, paste0(file_prefix, "_quantiles.tex")))
      }
    }

    # Acceptance Rate
    if ("acceptance" %in% params) {
      cat("\n ACCEPTANCE RATE \n")
      acceptance_rates <- 1 - rejectionRate(mcmc_chains)
      print(acceptance_rates)
      results_list$acceptance <- acceptance_rates
      if (save) {
        tex <- xtable::xtable(as.matrix(acceptance_rates), caption = "Acceptance Rates")
        print(tex, file = file.path(output_dir, paste0(file_prefix, "_acceptance.tex")))
      }
    }

    return(results_list)
  }))
}

#' Plot Posterior Distributions of Estimated Parameters
#'
#' Generates density or histogram plots for each estimated parameter from MCMC output, with annotated quantiles
#' and optional true values if available. It helps users visualize the posterior distribution of model parameters.
#'
#' @param results_list A list containing MCMC results from multiple runs. Each element should contain posterior samples.
#' @param params_to_estimate A character vector of parameter names to include in the plot.
#' @param dim_plot A numeric vector of length 2 specifying the number of rows and columns for arranging plots.
#' @param show_true Logical; if `TRUE`, true parameter values will be indicated on the plots if provided.
#' @param true_value A named numeric vector of true values for parameters (optional). Must match `params_to_estimate` if used.
#' @param show_prior Logical; if `TRUE`, prior distributions will be overlaid on the plots.
#' @param prior_n Integer; the number of points to sample from the prior distribution.
#' @param title A character string specifying the title for the entire plot layout (optional).
#' @param run_labels A character vector specifying labels for each MCMC run (optional). Must match the length of `results_list`.
#' @param plot_type A character string specifying the type of plot: `"histogram"` or `"density"`.
#'
#' @return A combined plot object displaying posterior distributions for each parameter.
#' The plot includes optional lines for true values if `show_true = TRUE`, and prior distributions if `show_prior = TRUE`.
#' @export
post_plot <- function(results_list, params_to_estimate, dim_plot, show_true = TRUE, true_value = NULL,
                      show_prior = FALSE, prior_n = 10000, title = "", run_labels = NULL, plot_type = "histogram") {
  # Check if the input is a valid list of MCMC results
  if (!is.list(results_list) || length(results_list) < 1) {
    stop("Error: 'results_list' should be a list containing MCMC results from multiple runs.")
  }

  # Validate plot type
  if (!(plot_type %in% c("histogram", "density"))) {
    stop("Error: 'plot_type' must be either 'histogram' or 'density'.")
  }

  # Assign default run labels if not provided
  if (is.null(run_labels)) {
    run_labels <- paste("Run", seq_along(results_list))
  } else if (length(run_labels) != length(results_list)) {
    stop("Error: Length of 'run_labels' must match the number of MCMC runs in 'results_list'.")
  }

  # Combine posterior samples from all MCMC runs
  all_posterior <- do.call(rbind, lapply(seq_along(results_list), function(i) {
    run_data <- data.frame(results_list[[i]][[1]]$pars)
    run_data$Run <- run_labels[i]  # Assign custom run label
    return(run_data)
  }))

  # Transform posterior data to long format for visualization
  post_melt <- reshape2::melt(all_posterior, id.vars = "Run")
  colnames(post_melt) <- c("Run", "variable", "value")
  post_melt$Type <- "posterior"
  post_melt$value.x <- post_melt$value

  # Remove any NA values from the posterior data
  post_melt <- na.omit(post_melt)

  # Handle prior distributions if enabled
  if (show_prior) {
    default_priors <- return_default_priors()
    default_priors <- default_priors[names(default_priors) %in% params_to_estimate]

    prior_melt <- data.frame(Run = character(), variable = character(), value = numeric(), Type = character(), value.x = numeric())

    for (param in names(params_to_estimate)) {
      prior <- default_priors[[param]]
      prior_fn <- prior$prior

      # Determine prior range based on percentiles of posterior values
      param_values <- post_melt$value[post_melt$variable == param]
      posterior_min <- quantile(param_values, 0.000001, na.rm = TRUE)
      posterior_max <- quantile(param_values, 0.999999, na.rm = TRUE)

      prior_x <- seq(posterior_min, posterior_max, length.out = prior_n)
      prior_y <- prior_fn(prior_x)

      # Remove infinite and NA values
      prior_y[is.infinite(prior_y) | is.na(prior_y)] <- 0

      # Add small constant to avoid losing flat priors in visualization

      # Estimate the density of posterior to scale prior accordingly
      density_est <- density(param_values, na.rm = TRUE)
      max_posterior_density <- max(density_est$y)

      # Ensure the prior is non-negative and scale relative to posterior distribution
      prior_y <- prior_y - min(prior_y, na.rm = TRUE) # Shift to non-negative space
      if(max(prior_y) == min(prior_y)){
        prior_y = rep(max_posterior_density * 0.8, length(prior_y))
      }else{
        prior_y <- ((prior_y) / (max(prior_y, na.rm = TRUE) - min(prior_y, na.rm = TRUE))) *
          max_posterior_density * 0.8
      }

      # Apply a safety threshold to prevent negative values
      prior_y <- pmax(prior_y, 0)


      prior_temp <- data.frame(Run = "Prior", variable = rep(param, length(prior_x)), value = prior_y, Type = "prior", value.x = prior_x)
      prior_melt <- rbind(prior_melt, prior_temp)
    }

    # Combine prior and posterior data for plotting
    prior_post <- rbind(prior_melt, post_melt)
  } else {
    prior_post <- post_melt
  }

  # Define color palette, adding red for the prior
  colors <- c(setNames(RColorBrewer::brewer.pal(min(length(run_labels), 8), "Set2"), run_labels), Prior = "red")

  # Generate plots for each parameter with LaTeX-formatted titles
  plot_list <- lapply(names(params_to_estimate), function(param) {
    if (!(param %in% prior_post$variable)) {
      warning(sprintf("Parameter '%s' not found in posterior samples.", param))
      return(NULL)
    }

    # Define LaTeX-formatted titles using latex2exp syntax
    latex_param_titles <- c(
      "lag_R" = "$j_2$",
      "lag_T" = "$j_1$",
      "qR" = "$q_{R}$",
      "phi" = "$\\phi$",
      "z" = "$z$",
      "sigma_LT" = "$\\sigma_{LT}$",
      "sigma_RT" = "$\\sigma_{RT}$",
      "k1" = "$\\tau$",
      "fT_C" = "$\\nu_1$",
      "size_1" = "$size_{1}$",
      "size_2" = "$size_{2}$"
    )

    # Get LaTeX-formatted title or default to param if not found
    plot_title <- ifelse(param %in% names(latex_param_titles),
                         latex_param_titles[[param]],
                         param)

    # Convert to LaTeX expression using TeX() function
    plot_title_latex <- latex2exp::TeX(plot_title)

    # Filter data for the current parameter
    post_data <- subset(prior_post, variable == param &
                          Type == "posterior")
    prior_data <- subset(prior_post, variable == param &
                           Type == "prior")

    p <- ggplot(post_data, aes(x = value, fill = Run, color = Run))

    # Choose between histogram and density plot
    if (plot_type == "histogram") {
      p <- p + geom_histogram(
        aes(y = ..density..),
        bins = 100,
        alpha = 0.5,
        color = "black",
        position = "identity"
      )
    } else if (plot_type == "density") {
      p <- p + geom_density(alpha = 0.5,
                            size = 1,
                            na.rm = TRUE)
    }

    p <- p +
      theme_bw() +
      labs(y = "Density", title = plot_title_latex) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "right"
      ) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      xlim(
        quantile(post_data$value, 0.01, na.rm = TRUE),
        quantile(post_data$value, 0.99, na.rm = TRUE)
      )

    # Add prior if enabled
    if (show_prior) {
      p <- p +
        geom_line(
          data = prior_data,
          aes(x = value.x, y = value, color = "Prior"),
          size = 1.2,
          alpha = 0.8
        ) +
        labs(title = plot_title_latex)
    }

    p
  })

  # Remove NULL plots and combine them
  plot_list <- Filter(Negate(is.null), plot_list)

  combined_plot <- patchwork::wrap_plots(plot_list, ncol = dim_plot[2]) +
    patchwork::plot_layout(guides = 'collect') +
    patchwork::plot_annotation(title = title,
                               theme = theme(plot.title = element_text(hjust = 0.5)))

  return(combined_plot)

}


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
