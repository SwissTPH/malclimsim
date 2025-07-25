###############################
## General Plotting Functions #
###############################
#' Plot Time Series of Malaria Incidence (Raw and Transformed)
#'
#' This function creates time series plots for malaria incidence, optionally including
#' transformed incidence values (e.g., for observation models). Also supports plotting
#' climate covariates, with faceting and overlay options.
#'
#' @param results A data frame from `data_sim()` containing columns like inc_C, inc_A,
#'                inc_C_transformed, inc_A_transformed, and date_ymd.
#' @param met Optional data frame of climate covariates.
#' @param plot_title Title of the plot.
#' @param incidence_y_label Y-axis label for incidence.
#' @param climate_y_label Y-axis label for climate plot.
#' @param climate_facet Whether to facet incidence and climate plots.
#' @param select_incidence Vector of incidence types to include. Options:
#'                         `"<5"`, `">=5"`, `"total"`, `"transformed_<5"`, `"transformed_>=5"`, `"transformed_total"`.
#' @param select_climate Climate vars to include (e.g., `"temp"`, `"cumrain"`).
#' @param incidence_colors Named vector of colors for incidence types.
#' @param climate_colors Named vector of colors for climate vars.
#' @param climate_alpha Transparency for climate lines.
#' @param base_size Base font size.
#'
#' @return A ggplot or grid of ggplots
#' @export
plot_time_series <- function(results, met = NULL,
                             plot_title = "Malaria Incidence Time Series",
                             incidence_y_label = "Malaria Incidence",
                             climate_y_label = "Climate",
                             climate_facet = FALSE,
                             select_incidence = c(">=5", "<5", "total"),
                             select_climate = c("temp", "cumrain"),
                             incidence_colors = c(">=5" = "blue", "<5" = "red",
                                                  "total" = "green",
                                                  "transformed_>=5" = "skyblue",
                                                  "transformed_<5" = "orange",
                                                  "transformed_total" = "darkgreen"),
                             climate_colors = c("temp" = "orange", "cumrain" = "purple"),
                             climate_alpha = 0.7,
                             base_size = 15) {

  # Prepare long-format incidence data
  incidence_vars <- list(
    ">=5" = "inc_A", "<5" = "inc_C", "total" = "inc",
    "transformed_>=5" = "inc_A_transformed",
    "transformed_<5" = "inc_C_transformed",
    "transformed_total" = "inc_transformed"
  )

  # Add inc_transformed (sum of transformed age groups) if not present
  if ("inc_A_transformed" %in% names(results) &&
      "inc_C_transformed" %in% names(results) &&
      !"inc_transformed" %in% names(results)) {
    results$inc_transformed <- results$inc_A_transformed + results$inc_C_transformed
  }

  valid_vars <- incidence_vars[select_incidence %in% names(incidence_vars)]
  selected_columns <- unique(unlist(valid_vars))
  selected_columns <- selected_columns[selected_columns %in% names(results)]

  results_long <- results %>%
    select(date_ymd, all_of(selected_columns)) %>%
    pivot_longer(-date_ymd, names_to = "var", values_to = "Incidence") %>%
    mutate(Type = recode(var,
                         inc_A = ">=5", inc_C = "<5", inc = "total",
                         inc_A_transformed = "transformed_>=5",
                         inc_C_transformed = "transformed_<5",
                         inc_transformed = "transformed_total")) %>%
    filter(Type %in% select_incidence)

  # Plot incidence
  p_incidence <- ggplot(results_long, aes(x = date_ymd, y = Incidence, color = Type)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = incidence_colors) +
    labs(title = plot_title, x = NULL, y = incidence_y_label, color = "Incidence Type") +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_x_date(date_breaks = "6 months", date_labels = "%b %Y")

  if (!is.null(met)) {
    met <- met %>%
      mutate(across(where(is.factor), as.character)) %>%
      mutate(dates = as.Date(date))

    met_monthly <- met %>%
      filter(dates %in% results$date_ymd) %>%
      select(dates, all_of(select_climate)) %>%
      pivot_longer(-dates, names_to = "Climate_Var", values_to = "Value")

    p_climate <- ggplot(met_monthly, aes(x = dates, y = Value, color = Climate_Var)) +
      geom_line(size = 1.2, alpha = climate_alpha) +
      scale_color_manual(values = climate_colors) +
      labs(x = NULL, y = climate_y_label, color = "Climate Variable") +
      theme_minimal(base_size = base_size) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_x_date(date_breaks = "6 months", date_labels = "%b %Y")

    if (climate_facet) {
      return(gridExtra::grid.arrange(p_incidence, p_climate, ncol = 1))
    } else {
      return(p_incidence + geom_line(data = met_monthly, aes(x = dates, y = Value, color = Climate_Var),
                                     inherit.aes = FALSE, size = 1, alpha = climate_alpha))
    }
  } else {
    return(p_incidence)
  }
}

###################################
## Compartment Plotting Functions #
###################################

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

#########################################
## Posterior Predictive Check Functions #
#########################################

#' Combine Posterior Simulations, Deterministic Prediction, and Observed Data for PPC Plotting
#'
#' @param ci_data Data frame with posterior median and credible intervals. Must include
#'   `date_ymd`, `variable`, `median`, `lower`, and `upper`.
#' @param best_df Long-format data frame from deterministic model run, including columns
#'   `date_ymd`, `variable`, and `value`.
#' @param obs_df Wide-format observed data with columns `date_ymd` and `obs_var`.
#' @param sim_var Character. Name of the variable used for simulated predictions (e.g., "inc_C_transformed").
#' @param obs_var Character. Name of the observed variable in `obs_df` (e.g., "inc_C").
#'
#' @return A tibble combining observed, deterministic, and simulated prediction data, ready for `plot_ppc()`.
#' @export
prepare_ppc_data <- function(ci_data, best_df, obs_df, sim_var, obs_var) {

  # -- Best-fit deterministic run
  best_filtered <- best_df %>%
    dplyr::filter(.data$variable == sim_var) %>%
    dplyr::mutate(label = "Model")

  # -- Observed data
  obs_filtered <- obs_df %>%
    dplyr::select(date_ymd, value = !!rlang::sym(obs_var)) %>%
    dplyr::mutate(variable = sim_var, label = "Observed")

  # -- Combine all three
  dplyr::bind_rows(best_filtered, obs_filtered)
}

#' Plot Posterior Predictive Check (PPC)
#'
#' @param ppc_data Data frame from `prepare_ppc_data()` with columns:
#'   `date_ymd`, `value`, `variable`, and `label`.
#' @param ci_data Optional data frame with columns: `date_ymd`, `lower`, `upper`, `variable`.
#' @param title Plot title.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param show_obs Logical. Whether to show observed data.
#' @param show_best Logical. Whether to show deterministic prediction.
#' @param show_ribbon Logical. Whether to show posterior CI ribbon.
#'
#' @return A `ggplot2` object.
#' @export
plot_ppc <- function(ppc_data,
                     ci_data     = NULL,
                     title       = "Posterior Predictive Check",
                     xlab        = "Date",
                     ylab        = "Value",
                     show_obs    = TRUE,
                     show_best   = TRUE,
                     show_ribbon = TRUE) {

  # Set label levels
  label_levels <- c("Observed", "Model")
  ppc_data$label <- factor(ppc_data$label, levels = label_levels)

  p <- ggplot2::ggplot(ppc_data, ggplot2::aes(x = date_ymd, y = value, color = label))

  # Posterior ribbon
  if (!is.null(ci_data) && show_ribbon) {
    ribbon_data <- ci_data %>%
      dplyr::filter(variable %in% unique(ppc_data$variable))

    p <- p + ggplot2::geom_ribbon(
      data        = ribbon_data,
      ggplot2::aes(x = date_ymd, ymin = lower, ymax = upper),
      fill        = "#D95F02",
      alpha       = 0.25,
      inherit.aes = FALSE
    )
  }

  # Add lines by type
  if (show_obs) {
    p <- p + ggplot2::geom_point(data = ppc_data %>% dplyr::filter(label == "Observed"),
                                 size = 1.5)
  }

  if (show_best) {
    p <- p + ggplot2::geom_line(data = ppc_data %>% dplyr::filter(label == "Model"),
                                linewidth = 1.2)
  }

  # Final plot formatting
  p <- p +
    ggplot2::scale_color_manual(values = c(
      "Observed"           = "#201110",
      "Model" = "#D95F02"
    )) +
    ggplot2::labs(title = title, x = xlab, y = ylab, color = "Source") +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "top"
    )

  return(p)
}


#' Plot Single Time Series Comparison (Uncomplicated or Severe Cases)
#'
#' Displays time series with optional credible intervals, observed points, and
#' consistent legend formatting across scenarios (e.g., With/Without SMC, June/July 2023, etc.)
#'
#' @param plot_data Data frame with columns: date_ymd (Date), value (numeric), label
#' @param ci_data   Optional data frame: date_ymd, lower, upper, label
#' @param obs_data  Optional observed data (must have date_ymd)
#' @param obs_col   Name of obs column in obs_data
#' @param plot_title Title
#' @param xlim      x‐axis limits (Date vector)
#' @param ylim      y‐axis limits (numeric vector)
#' @param severity  Either "Uncomplicated" or "Severe"
#' @param scale_severe Scaling factor for severe cases (used if no per-year mapping is provided)
#' @param scale_severe_by_year Named vector of multipliers by year (optional)
#'
#' @return ggplot2 plot object
#' @export
plot_ppc_single <- function(plot_data,
                            ci_data    = NULL,
                            obs_data   = NULL,
                            obs_col    = NULL,
                            plot_title = "Time Series Comparison",
                            xlim       = NULL,
                            ylim       = NULL,
                            severity   = "Uncomplicated",
                            scale_severe = 1,
                            scale_severe_by_year = NULL) {
  # Define full label set (for consistent legends)
  legend_levels <- c("Observed", "With SMC", "Without SMC", "SMC in 2019", "2023 June", "2023 July",
                     "2023 True SMC Timing")

  # Ensure label column exists
  if (!"label" %in% colnames(plot_data)) {
    plot_data <- plot_data %>% mutate(label = as.character(variable))
  }
  if (!is.null(ci_data) && !"label" %in% colnames(ci_data)) {
    ci_data <- ci_data %>% mutate(label = as.character(variable))
  }

  # Severity scaling helper
  get_mult <- function(d) {
    if (!is.null(scale_severe_by_year)) {
      m <- scale_severe_by_year[as.character(year(d))]
      ifelse(is.na(m), scale_severe, m)
    } else {
      scale_severe
    }
  }

  # Apply scaling if severe
  if (severity == "Severe") {
    plot_data <- plot_data %>% mutate(value = value * get_mult(date_ymd))
    if (!is.null(ci_data)) {
      ci_data <- ci_data %>%
        mutate(lower = lower * get_mult(date_ymd),
               upper = upper * get_mult(date_ymd))
    }
  }

  # Add observed data if present
  if (!is.null(obs_data) && !is.null(obs_col)) {
    obs_df <- obs_data %>%
      select(date_ymd, !!sym(obs_col)) %>%
      rename(value = !!sym(obs_col)) %>%
      mutate(label = "Observed")

    if (severity == "Severe") {
      obs_df <- obs_df %>% mutate(value = value * get_mult(date_ymd))
    }

    plot_data <- bind_rows(plot_data, obs_df)
  }

  # Enforce label as factor with fixed levels
  plot_data <- plot_data %>%
    mutate(label = factor(label, levels = legend_levels))

  if (!is.null(ci_data)) {
    ci_data <- ci_data %>%
      mutate(label = factor(label, levels = legend_levels))
  }

  # Determine which labels are present
  label_levels <- levels(droplevels(plot_data$label))
  present_idx <- legend_levels %in% label_levels

  # Define color palette
  color_palette <- c(
    "Observed"     = "black",
    "With SMC"     = "#6b6363",
    "Without SMC"  = "#201d69",
    "SMC in 2019"  = "#30baff",
    "2023 June"    = "#e41a1c",
    "2023 July"    = "#377eb8",
    "2023 True SMC Timing" = "#e41a1c"

  )

  # Define legend override aesthetics dynamically
  override_aes <- list(
    linetype = ifelse(legend_levels == "Observed", "blank", "solid")[present_idx],
    shape    = ifelse(legend_levels == "Observed", 16, NA)[present_idx],
    size     = ifelse(legend_levels == "Observed", 2.5, 1.2)[present_idx]
  )

  # Start plot
  p <- ggplot(plot_data, aes(x = date_ymd, y = value, color = label)) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle   = element_text(hjust = 0.5, face = "italic"),
      axis.text       = element_text(size = 14),
      axis.title      = element_text(size = 16),
      legend.position = "top",
      legend.title    = element_text(size = 14),
      legend.text     = element_text(size = 13)
    ) +
    labs(
      title    = plot_title,
      subtitle = severity,
      x        = "",
      y        = "Weekly Cases",
      color    = ""
    ) +
    scale_color_manual(
      values = color_palette[label_levels],
      breaks = label_levels,
      drop   = FALSE
    ) +
    guides(color = guide_legend(override.aes = override_aes))

  # Add CI ribbons
  if (!is.null(ci_data)) {
    p <- p + geom_ribbon(
      data = ci_data,
      aes(x = date_ymd, ymin = lower, ymax = upper, fill = label),
      alpha = 0.25,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
      scale_fill_manual(values = color_palette[label_levels], drop = FALSE)
  }

  # Plot lines and points
  p <- p +
    geom_line(data = filter(plot_data, label != "Observed"), size = 1.2, na.rm = TRUE) +
    geom_point(data = filter(plot_data, label == "Observed"), size = 2.5, alpha = 0.8, na.rm = TRUE)

  # Axis limits
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }

  return(p)
}
