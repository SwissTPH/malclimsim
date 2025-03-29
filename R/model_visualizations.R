#' Plot Time Series Data of Malaria Incidence and Climate Variables
#'
#' This function creates time series plots for malaria incidence and climate data,
#' with options to customize which incidence types (under 5 years, 5 years and older,
#' or total) and climate variables (temperature or rolling mean) are displayed.
#' The user can also choose to facet the plots or overlay them for better visualization.
#'
#' @param results A data frame containing malaria incidence data with columns for
#'                date (`date_ymd`), and incidence counts for different age groups (e.g.,
#'                `inc_A`, `inc_C`, `inc` for `<5`, `>=5`, and `total` respectively).
#' @param met A data frame containing climate data with columns for date (`date`),
#'            rainfall, temperature, and other relevant climate variables.
#'            The date column in this data must match the date format used in `results`.
#' @param plot_title A string for the title of the plot. Default is "Time Series Data".
#' @param incidence_y_label A string for the y-axis label of the incidence plot.
#'                           Default is "Monthly Malaria Cases".
#' @param climate_y_label A string for the y-axis label of the climate data plot.
#'                        Default is "Temperature (°C)".
#' @param climate_facet A boolean value indicating whether to facet the plots
#'                      (separate plots for incidence and climate data) or not.
#'                      Default is `FALSE` (plots will be overlaid).
#' @param select_incidence A vector of incidence types to include in the plot.
#'                         Options are `">=5"`, `"<5"`, and `"total"`. Default is all three.
#' @param select_climate A vector of climate variables to include in the plot.
#'                       Options are `"temp"` (temperature) and `"rollrain"` (rolling mean).
#'                       Default is both `"temp"` and `"rollrain"`.
#' @param incidence_colors A named vector of colors for the selected incidence types.
#'                         Default is `">=5" = "blue"`, `"<5" = "red"`, `"total" = "green"`.
#' @param climate_colors A named vector of colors for the selected climate variables.
#'                       Default is `"temp" = "orange"`, `"rollrain" = "purple"`.
#' @param climate_alpha A numeric value controlling the transparency of the climate data
#'                      lines. Default is 0.7 (semi-transparent).
#' @param base_size A numeric value controlling the base font size of the plot.
#'                  Default is 15.
#'
#' @return A ggplot object containing the time series plot(s) of malaria incidence and climate data.
#' @export
#'
#' @examples
#' # Example 1: Plotting only the total incidence of malaria
#' plot_time_series(results = results, plot_title = "Total Malaria Incidence", select_incidence = "total")
#'
#' # Example 2: Plotting malaria incidence and climate data (temperature only)
#' plot_time_series(results = results, met = met, plot_title = "Malaria and Temperature",
#'                  select_incidence = c(">=5", "<5"), select_climate = "temp")
#'
#' # Example 3: Plotting both incidence and climate data with faceting
#' plot_time_series(results = results, met = met, plot_title = "Malaria and Climate Data",
#'                  select_incidence = "total", select_climate = c("temp", "rollrain"),
#'                  climate_facet = TRUE)

plot_time_series <- function(results, met = NULL,
                             plot_title = "Time Series Data",
                             incidence_y_label = "Monthly Malaria Cases",
                             climate_y_label = "Temperature (°C)",
                             climate_facet = FALSE,
                             select_incidence = c(">=5", "<5", "total"),
                             select_climate = c("temp", "cumrain"),
                             incidence_colors = c(">=5" = "blue", "<5" = "red", "total" = "green"),
                             climate_colors = c("temp" = "orange", "cumrain" = "purple"),
                             climate_alpha = 0.7,
                             base_size = 15) {

  # Prepare the malaria data (assumed to be monthly)
  # Reshape the malaria incidence data into a long format, so that each row represents
  # a date and an incidence type, and we can easily filter the desired types.
  results_long <- results %>%
    select(date_ymd, inc_A, inc_C, inc) %>%
    pivot_longer(cols = c(inc_A, inc_C, inc),
                 names_to = "Type",   # Create a column 'Type' for age group types
                 values_to = "Incidence") %>%
    mutate(Type = recode(Type, inc_A = ">=5", inc_C = "<5", inc = "total")) %>%
    filter(Type %in% select_incidence)  # Filter based on the selected incidence types

  # Prepare the climate data, if provided
  if (!is.null(met)) {
    # Ensure `met` data's `zoo` columns are converted to numeric for plotting
    met <- met %>%
      mutate(
        cumrain = as.numeric(cumrain),  # Convert rollrain to numeric
        anom = as.numeric(anom),          # Convert anomaly to numeric (if necessary)
        dates = as.Date(date)  # Convert date to proper format for merging
      )

    # Filter the climate data to match the dates from `results` and the selected variables
    met_monthly <- met %>%
      filter(dates %in% results$date_ymd) %>%
      select(dates, all_of(select_climate)) %>%
      pivot_longer(cols = all_of(select_climate),
                   names_to = "Climate_Var",  # Create a column 'Climate_Var' for climate types
                   values_to = "Value")  # Column for climate values
  }

  # Plot the malaria incidence data
  p_incidence <- ggplot(data = results_long, aes(x = date_ymd, y = Incidence, color = Type)) +
    geom_line(size = 1.2) +  # Create the line plot for incidence
    scale_color_manual(values = incidence_colors) +  # Custom colors for incidence types
    labs(title = plot_title, x = NULL, y = incidence_y_label, color = "Incidence Type") +
    theme_minimal(base_size = base_size) +  # Set the minimal theme with custom base size
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
          plot.title = element_text(hjust = 0.5)) +  # Center the title
    scale_x_date(date_breaks = "6 months", date_labels = "%b %Y")  # Format x-axis as dates with 6-month breaks

  # Plot the climate data if provided
  if (!is.null(met)) {
    p_climate <- ggplot(data = met_monthly, aes(x = dates, y = Value, color = Climate_Var)) +
      geom_line(size = 1.2, alpha = climate_alpha) +  # Create the climate line plot
      scale_color_manual(values = climate_colors) +  # Custom colors for climate variables
      labs(x = NULL, y = climate_y_label, color = "Climate Variable") +
      theme_minimal(base_size = base_size) +  # Set the minimal theme for the climate plot
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +  # Rotate x-axis labels and center title
      scale_x_date(date_breaks = "6 months", date_labels = "%b %Y")  # Date formatting for x-axis

    # Arrange the plots if faceting is required
    if (climate_facet) {
      gridExtra::grid.arrange(p_incidence, p_climate, ncol = 1)  # Arrange as separate facets
    } else {
      p_incidence + geom_line(data = met_monthly, aes(x = dates, y = Value, color = Climate_Var))  # Overlay the plots
    }
  } else {
    print(p_incidence)  # If no climate data is provided, just print the incidence plot
  }
}

post_pred_plot <- function(results, sim_df, model, dates_sim, dates_obs, show_clim = TRUE,
                           met = NULL, title = "", ages = NULL, theme = NULL) {
  obs_data <- results$incidence_df # observed

  # observed data dataframe for plotting
  obs_df <- data.frame(
    date = sim_df$date,
    lower = rep(NA, nrow(sim_df)),
    upper = rep(NA, nrow(sim_df)),
    value = c(obs_data$inc_C, obs_data$inc_A, obs_data$inc),
    variable1 = sim_df$variable1,
    variable2 = rep("Observed", nrow(sim_df))
  )
  comb_df <- rbind(sim_df, obs_df)

  comb_df$variable1 <- forcats::fct_relevel(comb_df$variable1, "u5", "o5", "total")
  if(!is.null(ages)){
    comb_df <- subset(comb_df, variable1 %in% ages)
  }

  p1 <- comb_df %>% ggplot(aes(x = date, y = value)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable2), alpha = 0.4) +
    geom_line(data = subset(comb_df, variable2 == "Simulated")) +
    geom_point(data = subset(comb_df, variable2 == "Observed"), size = 2) +
    facet_wrap(~ variable1, ncol = 1, scales = "free") +
    theme_bw() +
    ggtitle(title) +
    ylab("Monthly Cases") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold"), # Facet label size
      legend.position = "none",
      legend.title = element_blank()
    ) + theme

  if(show_clim){
    clim_df <- create_clim_df(met)
    p_clim <- clim_df %>% ggplot(aes(x = date, y = value)) +
      geom_line(color = "blue", linewidth = 1.2) + theme_bw() +
      theme(
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
      )

    p1 <- (p1 + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())) /
      (p_clim + theme +
         theme(axis.title.y = element_blank())) +
      plot_layout(heights = c(3, 1))
  }
  return(p1)
}


# Sample n times the full posterior distribution and use these parameter values
# to simulate different trajectories. The 1st and 99th quantiles of incidence will be taken
# to construct a ribbon and the sample with the maximum likelihood will be used as
# a sample trajectory. Alongside, the real data will be plotted as points.
create_sim_df <- function(results, n, dates_sim, dates_obs, model){
  n_years <- length(seq(year(as.Date(dates_sim[1])), year(as.Date(dates_sim[2]))))
  n_months <- n_years * 12
  coda_pars <- results$coda_pars
  param_inputs <- results$param_inputs
  i_to_keep <- sample(1:nrow(coda_pars), n) # picking which rows to keep
  coda_pars_s <- coda_pars[i_to_keep,-c(1, 2, 3)]
  month_inc_C <- matrix(nrow = n, ncol = n_months)
  month_inc_A <- matrix(nrow = n, ncol = n_months)
  best_params <- coda_pars[which.max(coda_pars_s[, 2]), -c(1, 2, 3)]
  param_inputs[names(best_params)] <- best_params
  best_sim <- data_sim(model, param_inputs,
                       start_date = dates_sim[1], end_date = dates_sim[2],
                       month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)

  for (i in 1:nrow(coda_pars_s)) {
    param_values <- coda_pars_s[i, ]
    param_inputs[colnames(coda_pars_s)] <- param_values
    sim_results <- data_sim(model, param_inputs,
                            start_date = dates_sim[1], end_date = dates_sim[2],
                            month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)
    month_inc_C[i, ] <- sim_results$inc_C
    month_inc_A[i, ] <- sim_results$inc_A
  }
  quantiles_C <- apply(month_inc_C, 2, function(x) quantile(x, probs = c(0.01, 0.05, 0.50, 0.95, 0.99)))
  quantiles_A <- apply(month_inc_A, 2, function(x) quantile(x, probs = c(0.01, 0.05, 0.50, 0.95, 0.99)))
  quantiles_total <- apply(month_inc_A + month_inc_C, 2, function(x) quantile(x, probs = c(0.01, 0.05, 0.50, 0.95, 0.99)))

  # data frame for plotting with simulated data
  # lower is 1st quantile and upper is 99th quantile
  sim_df <- data.frame(date = rep(sim_results$date_ymd, 3),
                       lower = c(quantiles_C[1, ], quantiles_A[1, ], quantiles_total[1, ]),
                       upper = c(quantiles_C[5, ], quantiles_A[5, ], quantiles_total[5, ]),
                       #value = c(best_sim$inc_C, best_sim$inc_A, best_sim$inc),
                       value = c(quantiles_C[3, ], quantiles_A[3, ], quantiles_total[3, ]),
                       variable1 = c(
                         rep("u5", n_months),
                         rep("o5", n_months),
                         rep("total", n_months)
                       ),
                       variable2 = rep("Simulated", 3 * n_months))
  sim_df <- sim_df[which(sim_df$date %in% seq(as.Date(dates_obs[1]), as.Date(dates_obs[2]), by = "day")), ]
  return(sim_df)
}

#' Plot Observed vs Simulated Data with Quantiles Ribbon and Optional Climate Data Facet
#'
#' This function simulates data from a model using parameters that maximize the log posterior,
#' then plots the observed and simulated incidence data for comparison. Optionally, a ribbon
#' representing the 1st and 99th quantiles of additional simulations can be included to visualize uncertainty.
#' It also allows for optional inclusion of climate data as a separate facet.
#'
#' @param results The MCMC results object containing the parameter samples.
#' @param obs_cases A data frame of observed cases with columns for `date_ymd` (date), `inc_A`, `inc_C`, and `inc`.
#' @param start_date The start date for the simulation output as a `Date` object or character string.
#' @param end_date The end date for the simulation as a `Date` object or character string.
#' @param model The model function to simulate from.
#' @param add_ribbon Logical; if TRUE, adds a ribbon to the plot representing the 1st and 99th quantiles of the incidence data from additional simulations.
#' @param n_samples The number of samples to draw from the MCMC results for the quantiles ribbon.
#' @param groups A character vector specifying which groups to include in the plot (`inc_A`, `inc_C`, `inc`).
#' @param met Optional climate data frame to be included as a separate facet.
#' @param climate_facet Logical; if TRUE, adds climate data as a separate facet in the plot.
#' @param prewarm_years Integer; the number of years to simulate before the `start_date` for stabilization (default is 2 years).
#' @param days_per_year Integer; the number of days per year in the simulation (default is 360).
#' @return A ggplot object displaying observed data as points, simulated data as a line, and an optional ribbon for quantiles.
#' @export
#' @examples
#' # Assuming `results` contains MCMC output and `obs_cases` is the observed cases data
#' plot_observed_vs_simulated(results, obs_cases, start_date = "2014-01-01",
#'                            end_date = "2014-12-31", model = data_sim,
#'                            add_ribbon = TRUE, n_samples = 100, groups = c("inc_A", "inc_C", "inc"),
#'                            met = met_data, climate_facet = TRUE, prewarm_years = 2)
plot_observed_vs_simulated <- function(results, obs_cases, start_date, end_date, model,
                                       add_ribbon = TRUE, n_samples = 100, groups = c("inc_A", "inc_C", "inc"),
                                       met = NULL, climate_facet = FALSE, prewarm_years = 2, days_per_year = 360,
                                       plot_title = NULL, multiple_results = TRUE,
                                       mu_transform_A = NULL, mu_transform_C = NULL, covariate_matrix = NULL,
                                       legend_labels = c("Original Model - True Temp",
                                                         "Original Model - Constant Temp",
                                                         "Simplified EIR - True Temp") ) {

  prewarm_start_date <- paste0(year(as.Date(start_date)) - prewarm_years, "-", format(as.Date(start_date), "%m-%d"))

  all_sim_data <- list()

  # results_form_ok <- tryCatch(
  #   {
  #     class(results[[1]][[1]][[1]]) == "mcstate_pmcmc"
  #   },
  #   error = function(e) FALSE
  # )
  #
  # model_form_ok <- tryCatch(
  #   {
  #     all(class(model[[1]]) == c("dust_generator", "R6ClassGenerator"))
  #   },
  #   error = function(e) FALSE
  # )
  #
  # err_message <- paste(
  #   "Results Format Okay:", results_form_ok,
  #   "Model Format Okay:", model_form_ok,
  #   "If single result and/or model, try putting inside a list"
  # )
  #
  # if (!(results_form_ok | model_form_ok)) {
  #   stop(err_message)
  # }
  # # Ensure `model` has the same length as `results`
  # if(is.list(model) & ((length(model) != length(results)))){
  #     stop("The number of model must match the number of elements in results.")
  # }else{
  #   model <- list(model)
  #   results <- list(results)
  #   }

  # Iterate through each result set and model in results and model
  for (i in seq_along(results)) {
    if(multiple_results){
      result <- results[[i]]$results
      curr_model <- try(model[[i]], silent = TRUE)
      if(inherits(curr_model, "try-error")){
        stop("Unable to extract correct model. It is possible that `model' is in the wrong format.")
      }
    }else{
      result <- results[[i]]
      curr_model <- model
      }

    # Run the simulation with parameters that maximize log posterior
    sim_data <- try(sim_data <- simulate_with_max_posterior_params(
      result,
      start_date = start_date,
      end_date = end_date,
      model = curr_model,
      prewarm_years = prewarm_years,
      days_per_year = days_per_year,
      mu_transform_A = mu_transform_A,
      mu_transform_C = mu_transform_C,
      covariate_matrix = covariate_matrix
    )
    , silent = TRUE)
    if(inherits(sim_data, "try-error")){
      stop("Unable to simulate from sim_data function. It is possible that `results' is in the wrong format.")
      return(err_sim)
    }

    #sim_data <- simulate_with_max_posterior_params(result, start_date = start_date,
    #                                               end_date = end_date, model = curr_model,
    #                                               prewarm_years = prewarm_years,
    #                                               days_per_year = days_per_year)

    # Filter simulation data to the desired period
    sim_data <- sim_data[sim_data$date_ymd >= as.Date(start_date), ]
    #sim_data$source <- paste0("Sim_", i) # Label each dataset with a unique source ID
    sim_data$source <- legend_labels[i]
    all_sim_data[[i]] <- sim_data
  }

  # Combine all simulated datasets into one
  sim_data_combined <- do.call(rbind, all_sim_data)

  # Ensure date column compatibility with observed data
  colnames(obs_cases)[1] <- "date_ymd"

  # Merge observed data with combined simulated data
  combined_data <- merge(obs_cases, sim_data_combined, by = "date_ymd", suffixes = c("_obs", "_sim"))

  # If add_ribbon is TRUE, generate quantiles for each result set
  if (add_ribbon) {
    # Initialize a list to store quantile data for all sources
    all_quantiles <- list()

    for (i in seq_along(results)) {
      if(multiple_results){
        result <- results[[i]]$results
        curr_model <- try(model[[i]], silent = TRUE)
        if(inherits(curr_model, "try-error")){
          stop("Unable to extract correct model. It is possible that `model' is in the wrong format.")
        }
      }else{
        result <- results[[i]]
        curr_model <- model
      }
      #result <- results[[i]]$results
      #result <- results[[i]]$results
      #curr_model <- model[[i]]

      # Sample parameters and simulate model
      sampled_params <- sample_params(result, n_samples)
      simulations <- simulate_models(
        curr_model,
        result$param_inputs,
        sampled_params,
        start_date = start_date,
        end_date = end_date
      )

      # Calculate incidence quantiles for each group
      incidence_quantiles <- calculate_incidence_quantiles(simulations)
      incidence_quantiles$source <- legend_labels[i] # Label quantiles with source ID

      # Add quantiles to the list
      all_quantiles[[i]] <- incidence_quantiles
    }

    # Combine all quantiles into one data frame
    all_quantiles_combined <- do.call(rbind, all_quantiles)

    # Merge quantiles with combined data once
    combined_data <- merge(
      combined_data,
      all_quantiles_combined,
      by = c("date_ymd", "source"),
      all.x = TRUE
    )
  }



  # reshape2::melt the combined data for plotting
  combined_data_long <- reshape2::melt(combined_data, id.vars = c("date_ymd", "source"),
                             measure.vars = c("inc_A_obs", "inc_C_obs", "inc_obs",
                                              "inc_A_sim", "inc_C_sim", "inc_sim"))

  # Filter selected groups for both observed and simulated
  selected_groups <- unlist(lapply(groups, function(group) c(paste0(group, "_obs"), paste0(group, "_sim"))))
  combined_data_long <- combined_data_long[combined_data_long$variable %in% selected_groups, ]

  # Create a mapping for facet labels
  facet_labels <- c("inc_A" = "Incidence (>=5 years)", "inc_C" = "Incidence (<5 years)", "inc" = "Total Incidence")

  if (is.null(plot_title)) {
    plot_title <- "Observed vs Simulated Monthly Malaria Cases"
  }

  # Plotting observed and simulated data on the same panels
  p <- ggplot(combined_data_long, aes(x = date_ymd, group = source, color = source)) +
    geom_point(data = subset(combined_data_long, grepl("_obs$", variable)), aes(y = value),
               color = "blue", size = 2, alpha = 0.7) +  # Observed data as points
    geom_line(data = subset(combined_data_long, grepl("_sim$", variable)), aes(y = value),
              size = 1.5) +                              # Simulated data as line
    facet_wrap(~ gsub("(_obs|_sim)", "", variable), scales = "free", ncol = 1, labeller = as_labeller(facet_labels)) +
    labs(title = plot_title,
         x = "Date",
         y = "Number of Monthly Malaria Cases",
         color = "Simulation") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold")
    )

  if (add_ribbon) {
    ribbon_data_list <- list()  # Store ribbon data for each group separately

    for (group in groups) {
      # Extract relevant quantiles for the group
      ribbon_data <- combined_data[, c("date_ymd",
                                       "source",
                                       paste0(group, "_q005"),
                                       paste0(group, "_q995"))]
      colnames(ribbon_data) <- c("date_ymd", "source", "q005", "q995")
      ribbon_data$variable <- group  # Add the group as a variable for faceting

      ribbon_data_list[[group]] <- ribbon_data
    }

    # Combine all ribbon data into one data frame
    ribbon_data_combined <- do.call(rbind, ribbon_data_list)

    # Add ribbons to the plot, ensuring proper filtering for facets
    p <- p + geom_ribbon(
      data = ribbon_data_combined,
      aes(
        x = date_ymd,
        ymin = q005,
        ymax = q995,
        fill = source
      ),
      alpha = 0.3,
      inherit.aes = FALSE
    )
  }
  return(p)
}

plot_observed_vs_simulated <- function(results, obs_cases, start_date, end_date, model,
                                       add_ribbon = TRUE, n_samples = 100,
                                       groups = c("inc_A", "inc_C", "inc", "inc_C_transformed"),
                                       met = NULL, climate_facet = FALSE, prewarm_years = 2,
                                       days_per_year = 360, plot_title = NULL, multiple_results = TRUE,
                                       mu_transform_A = NULL, mu_transform_C = NULL, covariate_matrix = NULL,
                                       legend_labels = c("Original Model - True Temp",
                                                         "Original Model - Constant Temp",
                                                         "Simplified EIR - True Temp")) {
  prewarm_start_date <- paste0(year(as.Date(start_date)) - prewarm_years, "-", format(as.Date(start_date), "%m-%d"))
  all_sim_data <- list()

  for (i in seq_along(results)) {
    if (multiple_results) {
      result <- results[[i]]$results
      curr_model <- try(model[[i]], silent = TRUE)
      if (inherits(curr_model, "try-error")) stop("Unable to extract correct model.")
    } else {
      result <- results[[i]]
      curr_model <- model
    }

    sim_data <- try(simulate_with_max_posterior_params(
      result,
      start_date = start_date,
      end_date = end_date,
      model = curr_model,
      prewarm_years = prewarm_years,
      days_per_year = days_per_year,
      mu_transform_A = mu_transform_A,
      mu_transform_C = mu_transform_C,
      covariate_matrix = covariate_matrix
    ), silent = TRUE)

    if (inherits(sim_data, "try-error")) stop("Unable to simulate from sim_data function.")

    sim_data <- sim_data[sim_data$date_ymd >= as.Date(start_date), ]
    sim_data$source <- legend_labels[i]
    all_sim_data[[i]] <- sim_data
  }

  sim_data_combined <- do.call(rbind, all_sim_data)
  colnames(obs_cases)[1] <- "date_ymd"
  combined_data <- merge(obs_cases, sim_data_combined, by = "date_ymd", suffixes = c("_obs", "_sim"))

  # --- Melt all observed/simulated combinations of requested groups
  all_vars <- unlist(lapply(groups, function(g) c(paste0(g, "_obs"), paste0(g, "_sim"))))
  combined_data_long <- reshape2::melt(combined_data, id.vars = c("date_ymd", "source"),
                                       measure.vars = all_vars)

  # --- Set readable facet labels
  facet_labels <- c(
    "inc_A" = "Incidence (>=5 years)",
    "inc_C" = "Incidence (<5 years)",
    "inc"   = "Total Incidence",
    "inc_C_transformed" = "Incidence (<5) Transformed"
  )

  if (is.null(plot_title)) {
    plot_title <- "Observed vs Simulated Monthly Malaria Cases"
  }

  # --- Main plot
  p <- ggplot(combined_data_long, aes(x = date_ymd, group = source, color = source)) +
    geom_point(data = subset(combined_data_long, grepl("_obs$", variable)), aes(y = value),
               color = "blue", size = 2, alpha = 0.7) +
    geom_line(data = subset(combined_data_long, grepl("_sim$", variable)), aes(y = value), size = 1.5) +
    facet_wrap(~ gsub("(_obs|_sim)", "", variable), scales = "free", ncol = 1,
               labeller = as_labeller(facet_labels)) +
    labs(title = plot_title,
         x = "Date",
         y = "Number of Monthly Malaria Cases",
         color = "Simulation") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold")
    )

  # --- Add Ribbon if requested
  if (add_ribbon) {
    ribbon_data_list <- list()

    for (group in groups) {
      qcols <- c(paste0(group, "_q005"), paste0(group, "_q995"))
      if (!all(qcols %in% colnames(combined_data))) next

      ribbon_data <- combined_data[, c("date_ymd", "source", qcols)]
      colnames(ribbon_data) <- c("date_ymd", "source", "q005", "q995")
      ribbon_data$variable <- group
      ribbon_data_list[[group]] <- ribbon_data
    }

    ribbon_data_combined <- do.call(rbind, ribbon_data_list)

    p <- p + geom_ribbon(
      data = ribbon_data_combined,
      aes(x = date_ymd, ymin = q005, ymax = q995, fill = source),
      alpha = 0.3,
      inherit.aes = FALSE
    )
  }

  return(p)
}



# Update create_clim_df function to include temp and rollrain
create_clim_df <- function(clim_df){
  # Repeat the dates twice, once for each climate variable
  month_clim <- rep(as.Date(clim_df$date), 2)

  # Define the variable names and values for each climate variable
  variable <- c(rep("rollrain", nrow(clim_df)), rep("temp", nrow(clim_df)))
  value <- c(clim_df$rollrain, clim_df$temp)

  # Create the data frame with climate data
  clim_plot_df <- data.frame(date = month_clim, variable = variable, value = value)

  return(clim_plot_df)
}




#' Plot Residuals of Observed vs Simulated Data
#'
#' This function creates residual plots for observed and simulated data, allowing the user to select specific groups.
#'
#' @param observed_df Data frame containing the observed data.
#' @param simulated_df Data frame containing the simulated data.
#' @param date_column The name of the date column in both data frames.
#' @param groups A character vector specifying which groups to include in the plot (`inc_A`, `inc_C`, `inc`).
#' @return A ggplot object displaying the residual plots.
#' @export
#' @examples
#' residual_plot <- plot_residuals(obs_cases, simulated_df, date_column = "date_ymd", groups = c("inc_A", "inc_C"))
plot_residuals <- function(observed_df, simulated_df, date_column, groups = c("inc_A", "inc_C", "inc")) {

  # Ensure date compatibility and merge observed and simulated data
  observed_df[[date_column]] <- as.Date(observed_df[[date_column]])
  simulated_df[[date_column]] <- as.Date(simulated_df[[date_column]])
  combined_data <- merge(observed_df, simulated_df, by = date_column, suffixes = c("_obs", "_sim"))

  # Calculate residuals for each group: (observed - simulated)
  combined_data$residual_inc_A <- combined_data$inc_A_obs - combined_data$inc_A_sim
  combined_data$residual_inc_C <- combined_data$inc_C_obs - combined_data$inc_C_sim
  combined_data$residual_inc <- combined_data$inc_obs - combined_data$inc_sim

  # Prepare data for plotting residuals for each selected group
  residuals_long <- data.frame(
    date_ymd = combined_data[[date_column]],
    inc_A = combined_data$residual_inc_A,
    inc_C = combined_data$residual_inc_C,
    inc = combined_data$residual_inc
  ) %>%
    tidyr::pivot_longer(cols = c("inc_A", "inc_C", "inc"),
                        names_to = "Group",
                        values_to = "Residual") %>%
    dplyr::filter(Group %in% groups)  # Filter selected groups

  # Create facet plot for residuals
  residual_plot <- ggplot(residuals_long, aes(x = date_ymd, y = Residual)) +
    geom_line(color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Group, scales = "free_y", ncol = 1,
               labeller = labeller(Group = c(
                 "inc_A" = ">=5 years old (inc_A)",
                 "inc_C" = "<5 years old (inc_C)",
                 "inc" = "Total incidence (inc)"
               ))) +
    labs(
      title = "Residuals: Observed vs Simulated Incidence",
      x = "Date",
      y = "Residuals (Observed - Simulated)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 14, face = "bold")
    )

  return(residual_plot)
}

#' Plot SMC Coverage Scenarios
#'
#' This function plots the simulated malaria cases for different SMC coverage scenarios,
#' including "No SMC", "Current Coverage", and "Full SMC Coverage" scenarios.
#' It optionally includes observed data for comparison.
#'
#' @param observed_df A data frame of observed malaria cases with columns for `date_ymd`, `inc_A`, `inc_C`, and `inc`.
#' @param no_smc_df A data frame of simulated malaria cases with no SMC (`eff_SMC = 0`).
#' @param current_smc_df A data frame of simulated malaria cases with current SMC coverage.
#' @param full_smc_df A data frame of simulated malaria cases with full SMC coverage.
#' @param groups A character vector specifying which groups to include in the plot (`inc_A`, `inc_C`, `inc`).
#' @param include_observed Logical; if `TRUE`, includes observed data in the plot. Default is `TRUE`.
#' @return A ggplot object displaying the simulated data for different SMC scenarios (and optionally observed data).
#' @export
#'
#' @examples
#' # Assuming the scenarios have been simulated and stored as `observed`, `no_smc`, `current_smc`, `full_smc`
#' plot_smc_scenarios(observed, no_smc, current_smc, full_smc, groups = c("inc_A", "inc_C", "inc"), include_observed = TRUE)
plot_smc_scenarios <- function(observed_df, no_smc_df, current_smc_df, full_smc_df, groups = c("inc_A", "inc_C", "inc"), include_observed = TRUE) {

  # Ensure the date column is compatible and has the same name across all data frames
  no_smc_df$date_ymd <- as.Date(no_smc_df$date_ymd)
  current_smc_df$date_ymd <- as.Date(current_smc_df$date_ymd)
  full_smc_df$date_ymd <- as.Date(full_smc_df$date_ymd)

  # Rename columns to add suffixes to avoid conflicts
  colnames(no_smc_df)[-1] <- paste0(colnames(no_smc_df)[-1], "_no_smc")
  colnames(current_smc_df)[-1] <- paste0(colnames(current_smc_df)[-1], "_current_smc")
  colnames(full_smc_df)[-1] <- paste0(colnames(full_smc_df)[-1], "_full_smc")

  # Merge simulated data by date
  combined_data <- merge(no_smc_df, current_smc_df, by = "date_ymd")
  combined_data <- merge(combined_data, full_smc_df, by = "date_ymd")

  # If include_observed is TRUE, merge observed data as well
  if (include_observed) {
    observed_df$date_ymd <- as.Date(observed_df$date_ymd)
    colnames(observed_df)[-1] <- paste0(colnames(observed_df)[-1], "_obs")
    combined_data <- merge(observed_df, combined_data, by = "date_ymd")
  }

  # Convert the combined data into long format for ggplot2
  measure_vars <- unlist(lapply(groups, function(group) {
    c(
      if (include_observed) paste0(group, "_obs"),
      paste0(group, "_no_smc"),
      paste0(group, "_current_smc"),
      paste0(group, "_full_smc")
    )
  }))

  # Ensure the columns exist in the combined data
  measure_vars <- measure_vars[measure_vars %in% colnames(combined_data)]

  combined_data_long <- reshape2::melt(combined_data, id.vars = "date_ymd", measure.vars = measure_vars)

  # Extract the scenario from the variable name for easy labeling
  combined_data_long$Scenario <- gsub(".*_(obs|no_smc|current_smc|full_smc)$", "\\1", combined_data_long$variable)
  combined_data_long$Scenario <- factor(combined_data_long$Scenario, levels = c("obs", "no_smc", "current_smc", "full_smc")) # Ensure levels are consistent
  combined_data_long$Group <- gsub("(_obs|_no_smc|_current_smc|_full_smc)", "", combined_data_long$variable)

  # Filtering to only include the specified groups
  combined_data_long <- combined_data_long[combined_data_long$Group %in% groups, ]

  # Create descriptive labels for facet titles
  facet_labels <- c(
    "inc_A" = "Incidence (>=5 years old)",
    "inc_C" = "Incidence (<5 years old)",
    "inc" = "Total Incidence"
  )

  # Improve scenario labels for clarity
  scenario_labels <- c(
    "obs" = "Observed Data",
    "no_smc" = "No SMC Coverage",
    "current_smc" = "Current SMC Coverage",
    "full_smc" = "Full SMC Coverage (100%)"
  )

  # Plot the simulated malaria cases under different scenarios (and optionally observed data)
  p <- ggplot(combined_data_long, aes(x = date_ymd, y = value, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1.2) +  # Increased line thickness for better visibility
    facet_wrap(~ Group, scales = "free_y", ncol = 1, labeller = as_labeller(facet_labels)) +
    labs(
      title = "Malaria Cases Under Different SMC Scenarios",
      x = "Date",
      y = "Number of Monthly Malaria Cases",
      color = "Scenario",
      linetype = "Scenario"
    ) +
    scale_color_manual(
      values = c("obs" = "black",
                 "no_smc" = "darkred",
                 "current_smc" = "blue",
                 "full_smc" = "green4"),
      labels = scenario_labels
    ) +
    scale_linetype_manual(
      values = c("obs" = "solid",
                 "no_smc" = "dashed",
                 "current_smc" = "twodash",
                 "full_smc" = "dotted"),
      labels = scenario_labels
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12)
    )

  return(p)
}


#' Plot Compartmental Simulation Results
#'
#' This function visualizes the output of a compartmental model simulation with
#' options for line plots, stacked area plots, and faceted visualizations. Users
#' can also plot proportions relative to the total population.
#'
#' @param compart_df Data frame containing compartments simulation output.
#' @param plot_type Character string specifying the type of plot. Options are "line", "stacked", or "facet".
#' @param compartments Optional character vector of compartments to plot (e.g., c("SC", "IC")).
#' @param plot_proportion Logical indicating whether to plot proportions relative to the total population (`P` column).
#' @param log_scale Logical indicating whether to apply a logarithmic scale to the y-axis.
#' @param title Character string specifying the plot title.
#' @param color_palette Character string specifying the color palette from RColorBrewer (default: "Set2").
#' @param y_label Character string specifying the Y-axis label (default: "Population Count").
#'
#' @return A ggplot object visualizing the compartmental data.
#'
#' @examples
#' # Basic line plot
#' plot_compartments(compart_df, plot_type = "line")
#'
#' # Stacked area plot with selected compartments
#' plot_compartments(compart_df, plot_type = "stacked", compartments = c("SC", "IC", "RA"))
#'
#' # Faceted plot with proportions and log scale
#' plot_compartments(compart_df, plot_type = "line", plot_proportion = TRUE, log_scale = TRUE, facet = TRUE)
#'
#' @import ggplot2 reshape2 patchwork
#' @export
plot_compartments <- function(compart_df, plot_type = "line", compartments = NULL,
                              plot_proportion = FALSE, log_scale = FALSE,
                              title = "Epidemic Compartments Over Time",
                              color_palette = "Set2", y_label = "Population Count") {

  # Validate plot_type
  if (!(plot_type %in% c("line", "stacked", "facet"))) {
    stop("Invalid plot_type. Choose from 'line', 'stacked', or 'facet'.")
  }

  # Check if 'P' column exists for proportion calculation
  if (plot_proportion && !("P" %in% colnames(compart_df))) {
    stop("Column 'P' (total population) must be present in the data to plot proportions.")
  }

  # Convert data to long format for ggplot
  compart_long <- reshape2::melt(compart_df, id.vars = "date")

  # Filter for specified compartments if provided
  if (!is.null(compartments)) {
    compart_long <- compart_long[compart_long$variable %in% compartments, ]
  }

  # If plotting proportions, divide each compartment by total population (P)
  if (plot_proportion) {
    compart_long$value <- compart_long$value / compart_df$P[match(compart_long$date, compart_df$date)]
    y_label <- "Proportion of Population"
  }

  # Base ggplot object
  p <- ggplot(compart_long, aes(x = date, y = value, color = variable, fill = variable))

  # Choose plot type
  if (plot_type == "line") {
    p <- p + geom_line(size = 1.2) +
      scale_color_brewer(palette = color_palette) +
      labs(title = title, y = y_label, x = "Date")
  } else if (plot_type == "stacked") {
    p <- p + geom_area(alpha = 0.6, position = "stack") +
      scale_fill_brewer(palette = color_palette) +
      labs(title = title, y = y_label, x = "Date")
  } else if (plot_type == "facet" || facet) {
    p <- p + geom_line(size = 1) +
      facet_wrap(~ variable, ncol = 1, scales = "free_y") +
      scale_color_brewer(palette = color_palette) +
      labs(title = title, y = y_label, x = "Date")
  }

  # Apply logarithmic scale if specified
  if (log_scale) {
    p <- p + scale_y_log10() +
      annotation_logticks(sides = "l")  # Add log ticks on the left
  }

  # Add additional theme modifications
  p <- p + theme_minimal(base_size = 14) +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold")
    )

  return(p)
}

#' Plot Compartmental Simulation Results
#'
#' This function visualizes the output of a compartmental model simulation with
#' options for line plots, stacked area plots, and faceted visualizations. Users
#' can also plot proportions relative to the total population.
#'
#' @param compart_df Data frame containing compartments simulation output.
#' @param plot_type Character string specifying the type of plot. Options are "line", "stacked", or "facet".
#' @param compartments Optional character vector of compartments to plot (e.g., c("SC", "IC")).
#' @param plot_proportion Logical indicating whether to plot proportions relative to the total population (`P` column).
#' @param log_scale Logical indicating whether to apply a logarithmic scale to the y-axis.
#' @param title Character string specifying the plot title.
#' @param color_palette Character string specifying the color palette from RColorBrewer (default: "Set2").
#' @param y_label Character string specifying the Y-axis label (default: "Population Count").
#'
#' @return A ggplot object visualizing the compartmental data.
#'
#' @examples
#' # Basic line plot
#' plot_compartments(compart_df, plot_type = "line")
#'
#' # Stacked area plot with selected compartments
#' plot_compartments(compart_df, plot_type = "stacked", compartments = c("SC", "IC", "RA"))
#'
#' # Faceted plot with proportions and log scale
#' plot_compartments(compart_df, plot_type = "line", plot_proportion = TRUE, log_scale = TRUE, facet = TRUE)
#'
#' @import ggplot2 reshape2 patchwork
#' @export
plot_compartments <- function(compart_df, plot_type = "line", compartments = NULL,
                              plot_proportion = FALSE, log_scale = FALSE,
                              title = "Epidemic Compartments Over Time",
                              color_palette = "Set2", y_label = "Population Count") {

  # Validate plot_type
  if (!(plot_type %in% c("line", "stacked", "facet"))) {
    stop("Invalid plot_type. Choose from 'line', 'stacked', or 'facet'.")
  }

  # Convert data to long format for ggplot
  compart_long <- reshape2::melt(compart_df, id.vars = "date")

  # Filter for specified compartments if provided
  if (!is.null(compartments)) {
    compart_long <- compart_long[compart_long$variable %in% compartments, ]
  }

  # If plotting proportions, divide each compartment by the appropriate population
  if (plot_proportion) {
    compart_long$value <- with(compart_long, ifelse(grepl("C$", variable),
                                                    value / compart_df$PC[match(date, compart_df$date)],
                                                    ifelse(grepl("A$", variable),
                                                           value / compart_df$PA[match(date, compart_df$date)],
                                                           value / compart_df$P[match(date, compart_df$date)])))
    y_label <- "Proportion of Population"
  }

  # Base ggplot object
  p <- ggplot(compart_long, aes(x = date, y = value, color = variable, fill = variable))

  # Choose plot type
  if (plot_type == "line") {
    p <- p + geom_line(size = 1.2) +
      scale_color_brewer(palette = color_palette) +
      labs(title = title, y = y_label, x = "Date")
  } else if (plot_type == "stacked") {
    p <- p + geom_area(alpha = 0.6, position = "stack") +
      scale_fill_brewer(palette = color_palette) +
      labs(title = title, y = y_label, x = "Date")
  } else if (plot_type == "facet") {
    p <- p + geom_line(size = 1) +
      facet_wrap(~ variable, ncol = 1, scales = "free_y") +
      scale_color_brewer(palette = color_palette) +
      labs(title = title, y = y_label, x = "Date")
  }

  # Apply logarithmic scale if specified
  if (log_scale) {
    p <- p + scale_y_log10() +
      annotation_logticks(sides = "l")  # Add log ticks on the left
  }

  # Add additional theme modifications
  p <- p + theme_minimal(base_size = 14) +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold")
    )

  return(p)
}

#' Plot Estimated and True Prevalence Over Time
#'
#' This function creates a time series plot of estimated malaria prevalence.
#' Prevalence types are plotted on the same facet, while districts are shown as separate facets.
#' True prevalence estimates (from MAP data) can also be included and selected (`prev_A`, `prev_C`, `prev`).
#'
#' @param prevalence_data A data frame containing `year`, `District`, `Prevalence_Type`, and `Value` columns.
#' @param prevalence_true (Optional) A data frame containing observed prevalence estimates with columns `District`, `Year`, `prev_A`, `prev_C`, `prev`.
#' @param true_types A character vector specifying which of the true prevalence values to include (`prev_A`, `prev_C`, `prev`). Default is `prev`.
#' @param districts A character vector specifying which districts to include. Default is all.
#' @param prevalence_types A character vector specifying which prevalence types to include. Default is all.
#' @param title A character string specifying the plot title.
#' @param color_palette A character string specifying the color palette from RColorBrewer (default: "Set2").
#'
#' @return A ggplot2 object visualizing the prevalence trends over time.
#'
#' @details
#' - The function allows users to compare model-estimated prevalence against true values from malaria prevalence data (MAP).
#' - True prevalence data can be plotted for **all age groups (`prev`), only adults (`prev_A`), or only children (`prev_C`)**.
#' - If `prevalence_true` is provided, it is plotted as **circles** (`shape = 16`), while model estimates are displayed as **lines + dots**.
#' - If `prevalence_true` is NULL, only model estimates are plotted.
#'
#' @examples
#' # Plot only model estimates
#' plot_prevalence_over_time(prevalence_model_estimated)
#'
#' # Plot with true prevalence values (default = total prevalence)
#' plot_prevalence_over_time(prevalence_model_estimated, prevalence_true)
#'
#' # Compare model estimates to only `prev_A` (adults)
#' plot_prevalence_over_time(prevalence_model_estimated, prevalence_true, true_types = "prev_A")
#'
#' # Compare selected prevalence types
#' plot_prevalence_over_time(prevalence_model_estimated, prevalence_true,
#'                           prevalence_types = c("prev_total_with_R", "prev_A_with_R"),
#'                           true_types = c("prev_A", "prev_C"))
#'
#' @export
plot_prevalence_over_time <- function(prevalence_data,
                                      prevalence_true = NULL,
                                      true_types = "prev",  # Default to total prevalence
                                      districts = NULL,
                                      prevalence_types = NULL,
                                      title = "Estimated and True Prevalence Over Time",
                                      color_palette = "Set2") {

  # Convert year column to numeric for consistency
  prevalence_data$year <- as.numeric(prevalence_data$year)

  # Filter out "with_R" values and rename "no_R" types for merging with true values
  prevalence_data_clean <- prevalence_data %>%
    filter(!stringr::str_detect(Prevalence_Type, "with_R")) %>%  # Remove `with_R` types
    mutate(Prevalence_Type = case_when(
      Prevalence_Type == "prev_A_no_R" ~ "prev_A",
      Prevalence_Type == "prev_C_no_R" ~ "prev_C",
      Prevalence_Type == "prev_total_no_R" ~ "prev",
      TRUE ~ Prevalence_Type  # Keep other types unchanged
    )) %>%
    mutate(Source = "Model")  # Label model estimates

  # Process true prevalence data if provided
  if (!is.null(prevalence_true)) {
    prevalence_true <- prevalence_true %>% dplyr::select(-Value)  # Remove old `Value` column if exists
    prevalence_true$Year <- as.numeric(prevalence_true$Year)

    # Convert true prevalence data to long format
    prevalence_true_long <- prevalence_true %>%
      pivot_longer(cols = c(prev_A, prev_C, prev),
                   names_to = "Prevalence_Type",
                   values_to = "Value") %>%
      rename(year = Year) %>%
      mutate(Source = "Observed (MAP)")  # Label true values

    # Filter only selected true prevalence types
    prevalence_true_long <- prevalence_true_long %>%
      filter(Prevalence_Type %in% true_types)

    # Merge model estimates with true values (keep all model estimates)
    combined_data <- bind_rows(prevalence_data_clean, prevalence_true_long)
  } else {
    combined_data <- prevalence_data_clean  # Use only model estimates if no true data is provided
  }

  # Filter data based on user selection
  if (!is.null(districts)) {
    combined_data <- combined_data %>% filter(District %in% districts)
  }

  if (!is.null(prevalence_types)) {
    combined_data <- combined_data %>% filter(Prevalence_Type %in% prevalence_types)
  }

  # Base ggplot object
  p <- ggplot(combined_data, aes(x = year, y = Value, color = Prevalence_Type, group = interaction(Prevalence_Type, Source))) +
    geom_line(data = combined_data %>% filter(Source == "Model"), size = 1.2) +  # Model estimates as lines
    geom_point(data = combined_data %>% filter(Source == "Model"), size = 2, shape = 16) +  # Model estimates as dots
    geom_point(data = combined_data %>% filter(Source == "Observed (MAP)"), aes(color = Prevalence_Type), size = 3, shape = 16, stroke = 1.5, na.rm = TRUE) +  # True values as filled circles
    facet_wrap(~ District, scales = "free_y") +  # Separate facets by district
    scale_color_brewer(palette = color_palette) +
    labs(
      title = title,
      x = "Year",
      y = "Average Yearly Prevalence Estimate",
      color = "Prevalence Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold")
    )

  return(p)
}

#' Plot Confidence Intervals for Simulated Compartments
#'
#' This function generates a time series plot of median simulated values with confidence intervals,
#' supporting different plot types, proportion scaling, and log scaling, with consistent ribbon and line colors.
#'
#' @param summary_df A data frame containing summarized simulation results with confidence intervals.
#' @param plot_type Character string specifying the type of plot. Options are "line", "stacked", or "facet".
#' @param compartments Optional character vector of compartments to plot (e.g., c("SC", "IC")).
#' @param plot_proportion Logical indicating whether to plot proportions relative to the total population (`P` column).
#' @param log_scale Logical indicating whether to apply a logarithmic scale to the y-axis.
#' @param title Character string specifying the plot title (default: "Confidence Intervals for Simulated Compartments").
#' @param color_palette Character string specifying the color palette from RColorBrewer (default: "Set2").
#' @param y_label Character string specifying the Y-axis label (default: "Value").
#'
#' @return A ggplot object.
#' @examples
#' # Basic line plot
#' plot_simulation_confidence_intervals(summary_df, plot_type = "line")
#'
#' # Stacked area plot with selected compartments
#' plot_simulation_confidence_intervals(summary_df, plot_type = "stacked", compartments = c("SC", "IC"))
#'
#' # Faceted plot with proportions and log scale
#' plot_simulation_confidence_intervals(summary_df, plot_type = "facet", plot_proportion = TRUE, log_scale = TRUE)
#' @export
plot_simulation_confidence_intervals <- function(summary_df, plot_type = "line", compartments = NULL,
                                                 plot_proportion = FALSE, log_scale = FALSE,
                                                 title = "Confidence Intervals for Simulated Compartments",
                                                 color_palette = "Set2", y_label = "Value") {
  # Convert data to long format for flexible plotting
  long_df <- summary_df %>%
    pivot_longer(
      cols = -date,
      names_to = c("compartment", "stat"),
      names_pattern = "(.*)_(median|lower|upper)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(across(c(median, lower, upper), as.numeric))

  # Filter compartments if specified
  if (!is.null(compartments)) {
    long_df <- long_df %>% filter(compartment %in% compartments)
  }

  # Adjust for proportions if requested (assuming "P" exists in the data)
  if (plot_proportion && "P" %in% names(summary_df)) {
    long_df <- long_df %>%
      left_join(summary_df %>% select(date, P), by = "date") %>%
      mutate(across(c(median, lower, upper), ~ . / P)) %>%
      select(-P)
    y_label <- "Proportion of Population"
  }

  # Base ggplot object
  p <- ggplot(long_df, aes(x = date, group = compartment, color = compartment, fill = compartment))

  # Apply plot types
  if (plot_type == "line") {
    p <- p +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +  # Fill uses the same color as the line
      geom_line(aes(y = median), size = 1.2) +
      scale_color_brewer(palette = color_palette) +
      scale_fill_brewer(palette = color_palette)
  } else if (plot_type == "stacked") {
    p <- p +
      geom_area(aes(y = median), alpha = 0.6, position = "stack") +
      scale_fill_brewer(palette = color_palette)
  } else if (plot_type == "facet") {
    p <- p +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
      geom_line(aes(y = median), size = 1.2) +
      facet_wrap(~ compartment, ncol = 1, scales = "free_y") +
      scale_color_brewer(palette = color_palette) +
      scale_fill_brewer(palette = color_palette)
  } else {
    stop("Invalid plot_type. Choose from 'line', 'stacked', or 'facet'.")
  }

  # Apply log scale if requested
  if (log_scale) {
    p <- p + scale_y_log10() + annotation_logticks(sides = "l")
  }

  # Final plot adjustments
  p <- p +
    labs(
      title = title,
      x = "Date",
      y = y_label,
      fill = "Compartment",
      color = "Compartment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold")
    )

  return(p)
}

#' Plot Confidence Intervals for Simulated Compartments with Observed Data (Improved Legend & Colors)
#'
#' This function generates a time series plot of median simulated values with confidence intervals,
#' and overlays observed prevalence data as dots with consistent coloring and a clear legend.
#'
#' @param summary_df A data frame containing summarized simulation results with confidence intervals.
#' @param obs_data A data frame containing observed data with columns `date_ymd`, `prev_A`, and `prev_C`.
#' @param plot_type Character string specifying the type of plot. Options are "line", "stacked", or "facet".
#' @param compartments Optional character vector of compartments to plot (e.g., c("SC", "IC")).
#' @param plot_proportion Logical indicating whether to plot proportions relative to the total population (`P` column).
#' @param log_scale Logical indicating whether to apply a logarithmic scale to the y-axis.
#' @param title Character string specifying the plot title (default: "Confidence Intervals for Simulated Compartments").
#' @param color_palette Character string specifying the color palette from RColorBrewer (default: "Set2").
#' @param y_label Character string specifying the Y-axis label (default: "Value").
#'
#' @return A ggplot object.
#' @export
plot_simulation_confidence_intervals <- function(summary_df, obs_data = NULL, plot_type = "line", compartments = NULL,
                                                 plot_proportion = FALSE, log_scale = FALSE,
                                                 title = "Confidence Intervals for Simulated Compartments",
                                                 color_palette = "Set2", y_label = "Value") {

  library(tidyverse)

  # Convert summary_df to long format for flexible plotting
  long_df <- summary_df %>%
    pivot_longer(
      cols = -date,
      names_to = c("compartment", "stat"),
      names_pattern = "(.*)_(median|lower|upper)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(across(c(median, lower, upper), as.numeric))

  # Filter compartments if specified
  if (!is.null(compartments)) {
    long_df <- long_df %>% filter(compartment %in% compartments)
  }

  # Adjust for proportions if requested (assuming "P" exists in the data)
  if (plot_proportion && "P" %in% names(summary_df)) {
    long_df <- long_df %>%
      left_join(summary_df %>% select(date, P), by = "date") %>%
      mutate(across(c(median, lower, upper), ~ . / P)) %>%
      select(-P)
    y_label <- "Proportion of Population"
  }

  # Define a consistent color palette
  unique_compartments <- unique(long_df$compartment)
  colors <- RColorBrewer::brewer.pal(max(3, length(unique_compartments)), color_palette)
  color_map <- setNames(colors[1:length(unique_compartments)], unique_compartments)

  # Base ggplot object
  p <- ggplot(long_df, aes(x = date, group = compartment)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = compartment), alpha = 0.2) +
    geom_line(aes(y = median, color = compartment), size = 1.2) +
    scale_color_manual(values = color_map) +
    scale_fill_manual(values = color_map)

  # Add observed data as dots with legend integration
  if (!is.null(obs_data)) {
    obs_data <- obs_data %>% mutate(date_ymd = as.Date(date_ymd))

    # Plot observed prev_A
    if ("prev_A" %in% names(obs_data)) {
      p <- p + geom_point(data = obs_data %>% filter(!is.na(prev_A)),
                          aes(x = date_ymd, y = prev_A, shape = "Observed prev_A", color = "prev_A_with_R"),
                          size = 3, inherit.aes = FALSE)
    }

    # Plot observed prev_C
    if ("prev_C" %in% names(obs_data)) {
      p <- p + geom_point(data = obs_data %>% filter(!is.na(prev_C)),
                          aes(x = date_ymd, y = prev_C, shape = "Observed prev_C", color = "prev_C_with_R"),
                          size = 3, inherit.aes = FALSE)
    }

    # Add shapes for observed data in the legend
    p <- p + scale_shape_manual(name = "Data Type", values = c("Observed prev_A" = 16, "Observed prev_C" = 17))
  }

  # Apply plot types
  if (plot_type == "facet") {
    p <- p + facet_wrap(~ compartment, ncol = 1, scales = "free_y")
  }

  # Apply log scale if requested
  if (log_scale) {
    p <- p + scale_y_log10() + annotation_logticks(sides = "l")
  }

  # Final plot adjustments
  p <- p +
    labs(
      title = title,
      x = "Date",
      y = y_label,
      fill = "Simulated",
      color = "Simulated"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )

  return(p)
}

#' Plot Confidence Intervals for Simulated Compartments with Observed Data (Improved Legend & Colors)
#'
#' This function generates a time series plot of median simulated values with confidence intervals,
#' and overlays observed prevalence data as dots with consistent coloring and a clear legend.
#'
#' @param summary_df A data frame containing summarized simulation results with confidence intervals.
#' @param obs_data A data frame containing observed data with columns `date_ymd`, `prev_A`, and `prev_C`.
#' @param plot_type Character string specifying the type of plot. Options are "line", "stacked", or "facet".
#' @param compartments Optional character vector of compartments to plot (e.g., c("SC", "IC")).
#' @param plot_proportion Logical indicating whether to plot proportions relative to the total population (`P` column).
#' @param log_scale Logical indicating whether to apply a logarithmic scale to the y-axis.
#' @param title Character string specifying the plot title (default: "Confidence Intervals for Simulated Compartments").
#' @param color_palette Character string specifying the color palette from RColorBrewer (default: "Set2").
#' @param y_label Character string specifying the Y-axis label (default: "Value").
#'
#' @return A ggplot object.
#' @export
plot_simulation_confidence_intervals <- function(summary_df, obs_data = NULL, plot_type = "line", compartments = NULL,
                                                 plot_proportion = FALSE, log_scale = FALSE,
                                                 title = "Confidence Intervals for Simulated Compartments",
                                                 color_palette = "Set2", y_label = "Value") {

  library(tidyverse)

  # Convert summary_df to long format for flexible plotting
  long_df <- summary_df %>%
    pivot_longer(
      cols = -date,
      names_to = c("compartment", "stat"),
      names_pattern = "(.*)_(median|lower|upper)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(across(c(median, lower, upper), as.numeric))

  # Filter compartments if specified
  if (!is.null(compartments)) {
    long_df <- long_df %>% filter(compartment %in% compartments)
  }

  # Adjust for proportions if requested (assuming "P" exists in the data)
  if (plot_proportion && "P" %in% names(summary_df)) {
    long_df <- long_df %>%
      left_join(summary_df %>% select(date, P), by = "date") %>%
      mutate(across(c(median, lower, upper), ~ . / P)) %>%
      select(-P)
    y_label <- "Proportion of Population"
  }

  # Define a consistent color palette
  unique_compartments <- unique(long_df$compartment)
  colors <- RColorBrewer::brewer.pal(max(3, length(unique_compartments)), color_palette)
  color_map <- setNames(colors[1:length(unique_compartments)], unique_compartments)

  # Base ggplot object
  p <- ggplot(long_df, aes(x = date, group = compartment)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = compartment), alpha = 0.2) +
    geom_line(aes(y = median, color = compartment), size = 1.2) +
    scale_color_manual(values = color_map) +
    scale_fill_manual(values = color_map)

  # Add observed data as dots with correct facet assignment
  if (!is.null(obs_data)) {
    obs_data <- obs_data %>% mutate(date_ymd = as.Date(date_ymd))

    # Reshape observed data for plotting in facets
    obs_long <- obs_data %>%
      pivot_longer(cols = starts_with("prev_"), names_to = "compartment", values_to = "value") %>%
      filter(!is.na(value)) %>%
      mutate(compartment = case_when(
        compartment == "prev_A" ~ "prev_A_with_R",
        compartment == "prev_C" ~ "prev_C_with_R",
        TRUE ~ compartment
      ))

    # Add observed points to the plot, mapped to corresponding facets
    p <- p + geom_point(data = obs_long,
                        aes(x = date_ymd, y = value, shape = "Observed", color = compartment),
                        size = 3, inherit.aes = FALSE)

    # Add shape to legend
    p <- p + scale_shape_manual(name = "Data Type", values = c("Observed" = 16))
  }

  # Apply plot types
  if (plot_type == "facet") {
    p <- p + facet_wrap(~ compartment, ncol = 1, scales = "free_y")
  }

  # Apply log scale if requested
  if (log_scale) {
    p <- p + scale_y_log10() + annotation_logticks(sides = "l")
  }

  # Final plot adjustments
  p <- p +
    labs(
      title = title,
      x = "Date",
      y = y_label,
      fill = "Simulated",
      color = "Simulated"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )

  return(p)
}

#' Prepare Observed and Simulated Data for Plotting
#'
#' This function reshapes and combines simulated and observed data, allowing plotting of multiple
#' simulated variables (e.g., transformed and untransformed incidence) against the same observed variable.
#'
#' @param sim_data Simulated data frame in long format (must include 'date_ymd', 'variable', 'median', etc.).
#' @param obs_data Observed data frame in wide format (e.g., columns: date_ymd, inc_C, ...).
#' @param sim_vars Vector of variable names from simulation to include (e.g., c("inc_C", "inc_C_transformed")).
#' @param obs_var Name of the observed variable to match against (e.g., "inc_C").
#' @param label Label for the simulated data source.
#'
#' @return Combined long-format tibble ready for plotting.
#' @export
prepare_plot_data <- function(sim_data, obs_data, sim_vars, obs_var, label = "Simulated") {
  # --- Simulated data: filter selected variables and tag source
  sim_filtered <- sim_data %>%
    filter(variable %in% sim_vars) %>%
    mutate(source = label)

  # --- Expand observed data to match each sim variable
  repeated_obs <- lapply(sim_vars, function(var_name) {
    obs_data %>%
      select(date_ymd, !!obs_var) %>%
      rename(value = !!obs_var) %>%
      mutate(
        variable = var_name,
        source = "Observed"
      )
  })

  obs_long <- bind_rows(repeated_obs)

  # --- Combine and return
  bind_rows(sim_filtered, obs_long)
}


#' Plot Posterior Predictive Check (PPC)
#'
#' Generates a time series plot comparing observed vs simulated values with optional credible interval ribbons.
#'
#' @param plot_data Output of `prepare_plot_data()`: long-format tibble with `date_ymd`, `value`, `variable`, `source`.
#' @param ci_data Optional output of `summarize_simulations()`, long-format with `date`, `variable`, `lower`, `upper`.
#' @param plot_title Title for the plot.
#' @return A ggplot object.
#' @export
plot_ppc <- function(plot_data, ci_data = NULL, plot_title = "Posterior Predictive Check") {
  # Facet label mapping
  label_map <- c(
    inc_C = "No SMC",
    inc_C_transformed = "With SMC"
  )

  p <- ggplot(plot_data, aes(x = date_ymd, y = value, color = source, group = source)) +
    geom_point(data = subset(plot_data, source == "Observed"), size = 2.5, alpha = 0.8) +
    geom_line(data = subset(plot_data, source != "Observed"), size = 1.3) +
    facet_wrap(~ variable, scales = "free_y", labeller = as_labeller(label_map)) +
    labs(
      title = plot_title,
      x = "Date",
      y = "Monthly Number of Cases",
      color = "Data Source"
    ) +
    scale_color_manual(values = c("Observed" = "black", "Simulated" = "#1f77b4")) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      strip.text = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.position = "top",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    )

  if (!is.null(ci_data)) {
    p <- p + geom_ribbon(
      data = ci_data,
      aes(x = date_ymd, ymin = lower, ymax = upper, fill = variable),
      alpha = 0.25,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }

  return(p)
}


