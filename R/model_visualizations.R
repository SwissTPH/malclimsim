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
#'                       Options are `"temp"` (temperature) and `"rollmean"` (rolling mean).
#'                       Default is both `"temp"` and `"rollmean"`.
#' @param incidence_colors A named vector of colors for the selected incidence types.
#'                         Default is `">=5" = "blue"`, `"<5" = "red"`, `"total" = "green"`.
#' @param climate_colors A named vector of colors for the selected climate variables.
#'                       Default is `"temp" = "orange"`, `"rollmean" = "purple"`.
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
#'                  select_incidence = "total", select_climate = c("temp", "rollmean"),
#'                  climate_facet = TRUE)

plot_time_series <- function(results, met = NULL,
                             plot_title = "Time Series Data",
                             incidence_y_label = "Monthly Malaria Cases",
                             climate_y_label = "Temperature (°C)",
                             climate_facet = FALSE,
                             select_incidence = c(">=5", "<5", "total"),
                             select_climate = c("temp", "rollmean"),
                             incidence_colors = c(">=5" = "blue", "<5" = "red", "total" = "green"),
                             climate_colors = c("temp" = "orange", "rollmean" = "purple"),
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
        rollmean = as.numeric(rollmean),  # Convert rollmean to numeric
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

create_clim_df <- function(clim_df){
  month_clim <- rep(as.Date(clim_df$date), 2)
  variable <- c(rep("rain", nrow(clim_df)))
  value <- c(clim_df$anom)
  variable2 <- rep("Climate", length(value))
  upper <- rep(NA, length(value))
  lower <- rep(NA, length(value))
  clim_plot_df <- data.frame(date = month_clim, variable = variable, value = value,
                             variable2 = variable2, upper = upper, lower = lower)
  clim_plot_df <- data.frame(date = month_clim, lower = lower, upper = upper,
                             value = value, variable1 = variable, variable2 = variable2)
  return(clim_plot_df)
}
