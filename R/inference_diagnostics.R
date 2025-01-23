# Helper function to create a list of MCMC chains
# - Splits the MCMC output by chain for easier visualization and diagnostics
plot_chains <- function(mcmc_run) {
  chains <- mcmc_run$chain
  df <- data.frame(as.matrix(mcmc_run$pars), chains) # Combine parameter values with chain labels
  chains_sep <- split(df, df$chains) # Separate data by chain
  chains_list <- mcmc.list(lapply(chains_sep, function(x) as.mcmc(x[-ncol(x)]))) # Convert each chain to an 'mcmc' object
  return(chains_list)
}

# Helper function to create a correlation plot for MCMC samples
# - Used within MCMC_diag to show correlation structure of posterior samples
#' Title
#'
#' @param results results from inf_run function
#' @param title title of plot
#'
#' @return
#' @export
#'
#' @examples
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


#' MCMC Diagnostics
#'
#' Provides diagnostic information for MCMC results, including trace plots,
#' chain convergence statistics, correlation structure, effective sample size,
#' autocorrelation function, posterior quantiles, and acceptance rates.
#' Designed for users to assess MCMC chain convergence and parameter mixing.
#'
#' @param results A list containing the MCMC output, typically with at least
#' `results[[1]]` for chain data and `results[[2]]` for trace plot data.
#' Should include an element `n_chains` for the number of chains.
#' @param params A character vector specifying which diagnostics to display.
#' Options include: "trace", "gelman", "corr", "ess", "acf", "quantiles", "acceptance".
#' @param thin An integer specifying thinning interval for the chains.
#' Thinning reduces the number of samples by keeping every `thin`th sample.
#'
#' @return Diagnostic outputs based on the selected `params` are printed or plotted.
#'
#' @export
#'
#' @examples
#' # Assuming 'results' is a valid MCMC result list with required structure:
#' MCMC_diag(results, params = c("trace", "ess"))
MCMC_diag <- function(results, params = c("trace", "gelman", "corr", "ess", "acf", "quantiles", "acceptance"), thin = 1) {
  suppressWarnings(suppressMessages({
    # Extract parameter samples and split into chains
    pars <- results[[1]]$pars
    n_chains <- results$n_chains
    n_samples <- nrow(pars) / n_chains

    # Split parameters into individual chains
    chain_list <- split(pars, rep(1:n_chains, each = n_samples))
    chains <- lapply(chain_list, function(x) {
      as.mcmc(matrix(x, nrow = n_samples, ncol = ncol(pars), dimnames = list(NULL, colnames(pars))))
    })
    mcmc_chains <- as.mcmc.list(chains)

    # Apply thinning if requested
    if (thin > 1) {
      mcmc_chains <- lapply(mcmc_chains, function(chain) window(chain, thin = thin))
      mcmc_chains <- as.mcmc.list(mcmc_chains)
    }

    # Trace plot
    if ("trace" %in% params) {
      cat("\n TRACE PLOT \n")
      plot(mcmc_chains)  # Trace plot for convergence assessment
    }

    # Gelman-Rubin statistic
    if ("gelman" %in% params && n_chains > 1) {
      cat("\n GELMAN-RUBIN STATISTIC \n")
      print(gelman.diag(mcmc_chains))  # Numerical convergence assessment
    }

    # Correlation plot
    if ("corr" %in% params) {
      cat("\n CORRELATION PLOT \n")
      print(plot_corr(results))  # Retain the previous large and clear plot
    }

    # Effective Sample Size
    if ("ess" %in% params) {
      cat("\n EFFECTIVE SAMPLE SIZE \n")
      print(effectiveSize(mcmc_chains))
    }

    # Autocorrelation Function (ACF)
    if ("acf" %in% params) {
      cat("\n AUTOCORRELATION FUNCTION \n")
      acf(as.matrix(do.call(rbind, mcmc_chains)), main = "Autocorrelation")
    }

    # Posterior Quantiles
    if ("quantiles" %in% params) {
      cat("\n POSTERIOR QUANTILES \n")
      quantiles <- apply(as.matrix(do.call(rbind, mcmc_chains)), 2, quantile, probs = c(0.025, 0.5, 0.975))
      print(quantiles)
    }

    # Acceptance Rate
    if ("acceptance" %in% params) {
      cat("\n ACCEPTANCE RATE \n")
      acceptance_rates <- 1 - rejectionRate(mcmc_chains)
      print(acceptance_rates)
    }
  }))
}

MCMC_diag <- function(results, params = c("trace", "gelman", "corr", "ess", "acf", "quantiles", "acceptance"), thin = 1) {
  suppressWarnings(suppressMessages({
    # Extract parameter samples and split into chains
    pars <- results[[1]]$pars
    n_chains <- results$n_chains
    n_samples <- nrow(pars) / n_chains

    # Split parameters into individual chains
    chain_list <- split(pars, rep(1:n_chains, each = n_samples))
    chains <- lapply(chain_list, function(x) {
      as.mcmc(matrix(x, nrow = n_samples, ncol = ncol(pars), dimnames = list(NULL, colnames(pars))))
    })
    mcmc_chains <- as.mcmc.list(chains)

    # Apply thinning if requested
    if (thin > 1) {
      mcmc_chains <- lapply(mcmc_chains, function(chain) window(chain, thin = thin))
      mcmc_chains <- as.mcmc.list(mcmc_chains)
    }

    # Initialize a list to store results
    results_list <- list()

    # Trace plot
    if ("trace" %in% params) {
      cat("\n TRACE PLOT \n")
      plot(mcmc_chains)  # Trace plot for convergence assessment
      results_list$trace <- mcmc_chains
    }

    # Gelman-Rubin statistic
    if ("gelman" %in% params && n_chains > 1) {
      cat("\n GELMAN-RUBIN STATISTIC \n")
      gelman_result <- gelman.diag(mcmc_chains)  # Numerical convergence assessment
      print(gelman_result)
      results_list$gelman <- gelman_result
    }

    # Correlation plot
    if ("corr" %in% params) {
      cat("\n CORRELATION PLOT \n")
      plot_corr(results)  # Retain the previous large and clear plot
      results_list$corr <- "Correlation plot generated"
    }

    # Effective Sample Size
    if ("ess" %in% params) {
      cat("\n EFFECTIVE SAMPLE SIZE \n")
      ess_result <- effectiveSize(mcmc_chains)
      print(ess_result)
      results_list$ess <- ess_result
    }

    # Autocorrelation Function (ACF)
    if ("acf" %in% params) {
      cat("\n AUTOCORRELATION FUNCTION \n")
      acf_result <- acf(as.matrix(do.call(rbind, mcmc_chains)), main = "Autocorrelation")
      results_list$acf <- "Autocorrelation function plotted"
    }

    # Posterior Quantiles
    if ("quantiles" %in% params) {
      cat("\n POSTERIOR QUANTILES \n")
      quantiles <- apply(as.matrix(do.call(rbind, mcmc_chains)), 2, quantile, probs = c(0.025, 0.5, 0.975))
      results_list$quantiles <- quantiles
    }

    # Acceptance Rate
    if ("acceptance" %in% params) {
      cat("\n ACCEPTANCE RATE \n")
      acceptance_rates <- 1 - rejectionRate(mcmc_chains)
      print(acceptance_rates)
      results_list$acceptance <- acceptance_rates
    }

    # Return all results
    return(results_list)
  }))
}




#' Plot Posterior Distributions of Estimated Parameters
#'
#' Generates histogram plots for each estimated parameter in the MCMC output, with annotated quantiles and
#' optional true values if available. Designed for users to visualize the posterior distribution of model parameters.
#'
#' @param results A list containing MCMC results, where `results[[1]]` includes sampled parameter values.
#' @param params_to_estimate A character vector of parameter names to include in the plot.
#' @param dim_plot A numeric vector of length 2 specifying the number of rows and columns for arranging plots.
#' @param show_true Logical; if `TRUE`, true parameter values will be indicated on the plots if provided.
#' @param true_value A named vector of true values for parameters (optional). Must match `params_to_estimate` if used.
#' @param title A character string specifying the title for the entire plot layout (optional).
#' @param show_prior
#' @param prior_n
#'
#' @return A combined plot object displaying histograms of posterior distributions for each parameter,
#' with quantiles annotated. The plot includes optional lines for true values if `show_true = TRUE`
#' and `true_value` is provided.
#'
#' @export
#'
#' @examples
#' # Assuming 'results' contains MCMC posterior samples and parameter names are specified:
#' params_to_estimate <- c("param1", "param2")
#' dim_plot <- c(1, 2) # Arrange in 1 row, 2 columns
#' post_plot(results, params_to_estimate, dim_plot, show_true = TRUE, true_value = c(param1 = 0.5, param2 = -0.2), title = "Posterior Distributions")
post_plot <- function(results, params_to_estimate, dim_plot, show_true = TRUE, true_value = NULL, show_prior = FALSE, prior_n = 10000, title = "") {
  # Check if params_to_estimate is a named vector
  if (is.null(names(params_to_estimate)) || any(names(params_to_estimate) == "")) {
    stop("Error: 'params_to_estimate' must be a named vector where each parameter has an associated name.")
  }

  # Extract posterior samples
  posterior <- data.frame(results[[1]]$pars)
  post_melt <- reshape2::melt(posterior)
  colnames(post_melt) <- c("variable", "value")
  post_melt$Type <- "posterior"
  post_melt$value.x <- post_melt$value

  # Create prior samples if show_prior is TRUE
  if (show_prior) {
    default_priors <- return_default_priors()
    default_priors <- default_priors[names(default_priors) %in% params_to_estimate]

    prior_melt <- data.frame(variable = character(),
                             value = numeric(),
                             Type = character(),
                             value.x = numeric())

    for (param in names(params_to_estimate)) {
      prior <- default_priors[[param]]
      prior_fn <- prior$prior

      # Restrict prior range based on 0.5th and 99.5th percentiles of posterior
      param_values <- post_melt$value[post_melt$variable == param]
      posterior_min <- quantile(param_values, 0.000001)
      posterior_max <- quantile(param_values, 0.999999)

      prior_x <- seq(posterior_min, posterior_max, length.out = prior_n)
      prior_y <- prior_fn(prior_x)

      # Fix potential negative values and infinities
      prior_y[is.infinite(prior_y)] <- 0

      # Estimate posterior density to help with scaling
      density_est <- density(param_values)
      max_posterior_density <- max(density_est$y)
      posterior_mode <- density_est$x[which.max(density_est$y)]

      # Adjust scaling of prior using posterior density estimate
      prior_shift <- min(prior_y) - 0.2
      prior_y <- prior_y - prior_shift  # Shift to positive if necessary
      prior_y <- (prior_y / max(prior_y)) * max_posterior_density * 0.8  # Scale based on density

      # Restrict the prior range to match posterior range
      prior_x <- prior_x[prior_x >= posterior_min & prior_x <= posterior_max]
      prior_y <- prior_y[1:length(prior_x)]

      prior_temp <- data.frame(variable = rep(param, length(prior_x)),
                               value = prior_y,
                               Type = "prior",
                               value.x = prior_x)
      prior_melt <- rbind(prior_melt, prior_temp)
    }

    # Combine prior and posterior data correctly
    prior_post <- rbind(prior_melt, post_melt)
  } else {
    prior_post <- post_melt
  }

  # Calculate quantiles for posterior distributions
  quantiles <- apply(posterior, 2, function(x) quantile(x, probs = c(0.005, 0.025, 0.50, 0.975, 0.995)))
  quantiles_df <- as.data.frame(t(quantiles))
  colnames(quantiles_df) <- c("0.5%", "2.5%", "50%", "97.5%", "99.5%")
  rownames(quantiles_df) <- colnames(posterior)

  # Prepare true values if provided
  if (show_true && !is.null(true_value)) {
    vertical_lines <- data.frame(variable = names(true_value), line_position = as.vector(true_value))
    vertical_lines <- vertical_lines[vertical_lines$variable %in% names(params_to_estimate), ]
  }

  # Generate a plot for each parameter
  plot_list <- lapply(names(params_to_estimate), function(param) {
    if (!(param %in% colnames(posterior))) {
      warning(sprintf("Parameter '%s' not found in posterior samples.", param))
      return(NULL)
    }

    # Filter data for current parameter
    post_data <- subset(prior_post, Type == "posterior" & variable == param)

    p <- ggplot() +
      geom_histogram(data = post_data, aes(x = value, y = ..density.., fill = Type),
                     bins = 100, alpha = 0.5, color = "black") +
      theme_bw() +
      labs(x = param, y = "Density", title = param) +
      theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
      scale_fill_manual(values = c("posterior" = "blue")) +
      xlim(quantile(post_data$value, 0.00001), quantile(post_data$value, 0.99999)) +  # Exclude outliers
      scale_x_continuous(breaks = scales::extended_breaks(n = 3)) +
      scale_y_continuous(breaks = scales::extended_breaks(n = 3)) +
      geom_vline(xintercept = quantiles_df[param, "2.5%"], color = "blue", linetype = "dashed", size = 0.5) +
      geom_vline(xintercept = quantiles_df[param, "50%"], color = "blue", linetype = "solid", size = 0.5) +
      geom_vline(xintercept = quantiles_df[param, "97.5%"], color = "blue", linetype = "dashed", size = 0.5)

    # Add prior as a shaded area if show_prior is TRUE
    if (show_prior) {
      prior_data <- subset(prior_post, Type == "prior" & variable == param)
      p <- p +
        geom_ribbon(data = prior_data, aes(x = value.x, ymin = 0, ymax = value, fill = Type), alpha = 0.3) +
        scale_fill_manual(values = c("prior" = "red", "posterior" = "blue")) +
        xlim(quantile(post_data$value, 0.00001), quantile(post_data$value, 0.99999)) +  # Exclude outliers
        labs(title = param)
    }

    # Add true value line if specified
    if (show_true && !is.null(true_value)) {
      p <- p + geom_vline(data = subset(vertical_lines, variable == param), aes(xintercept = line_position),
                          color = "black", linetype = "solid", size = 0.8)
    }

    p
  })

  # Remove NULL plots (if parameters are missing)
  plot_list <- Filter(Negate(is.null), plot_list)

  # Combine individual parameter plots into a grid
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = dim_plot[2]) +
    plot_annotation(title = title, theme = theme(plot.title = element_text(hjust = 0.5)))

  return(combined_plot)
}


