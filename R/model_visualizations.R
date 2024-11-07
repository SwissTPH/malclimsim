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

