#' Plot response curves panel (applicable for any number of time periods)
#'
#' @param model_list Named list. Contains model results for each time period,
#'   e.g., `list(Quarter1 = sdm_q1, Quarter2 = sdm_q2, ...)`. Each element must
#'   contain a data frame `res_cur` with columns: `variable`, `Value`, `Response`,
#'   `lower`, `upper`.
#' @param time_names Character vector. Names of time periods, e.g.,
#'   `c("1st Quarter", "2nd Quarter", "3rd Quarter", "4th Quarter")`.
#' @param colors Optional character vector of colors, length equal to number of
#'   time periods. Defaults to `RColorBrewer::brewer.pal(3, "Set3")` (or rainbow
#'   if more than 12).
#' @param ncol Integer. Number of columns in the panel. Default: 4.
#' @param nrow Integer. Number of rows in the panel. Automatically calculated if `NULL`.
#' @param y_label Character. y‑axis label. Default: `"Habitat suitability"`.
#' @param legend_position Character. Legend position: `"vertical"` (stacked legend
#'   in the last panel) or `"bottom"` (unified legend at the bottom). Default: `"vertical"`.
#' @param width Numeric. Plot width in inches for saving. Default: 8.
#' @param height Numeric. Plot height in inches for saving. Default: 4.
#' @param save_path Character. Path to save the plot (TIFF format). If `NULL`,
#'   the plot is not saved.
#'
#' @return A combined `ggplot` object (patchwork or cowplot).
#' @export
#'
#' @examples
#' \donttest{
#' # Prepare model list (four quarters) containing the model object produced by TB_sdm()
#' model_list <- list(
#'   Q1 = sdm_q1,
#'   Q2 = sdm_q2,
#'   Q3 = sdm_q3,
#'   Q4 = sdm_q4
#' )
#'
#' # Time period names
#' time_names <- c("1st Quarter", "2nd Quarter", "3rd Quarter", "4th Quarter")
#'
#' # Plot response curves panel (vertical stacked legend)
#' res_plot <- TB_rescur(model_list,
#'                       time_names = time_names,
#'                       legend_position = "vertical",
#'                       save_path = "response_curves.tif")
#' }
TB_rescur <- function(model_list,
                      time_names = NULL,
                      colors = NULL,
                      ncol = 4,
                      nrow = NULL,
                      y_label = "Habitat suitability",
                      legend_position = c("vertical", "bottom"),
                      width = 8,
                      height = 4,
                      save_path = NULL) {

  # ---- Argument validation ----
  legend_position <- match.arg(legend_position)

  if (!is.list(model_list) || length(model_list) == 0) {
    stop("model_list must be a non-empty list")
  }

  n_times <- length(model_list)

  # If time names not provided, use default names
  if (is.null(time_names)) {
    time_names <- paste0("Time ", 1:n_times)
  }

  # If colors not provided, use RColorBrewer Set3 palette or rainbow
  if (is.null(colors)) {
    if (n_times <= 12) {
      colors <- RColorBrewer::brewer.pal(max(3, n_times), "Set3")
    } else {
      colors <- grDevices::rainbow(n_times)
    }
  }

  # ---- Combine response curve data from all time periods ----
  df_list <- list()
  for (i in 1:n_times) {
    model <- model_list[[i]]
    time_name <- time_names[i]

    # Extract response curve data
    if (!is.null(model$res_cur)) {
      df_time <- model$res_cur
    } else {
      stop(paste("Model at temporal bin", i, "does not contain res_cur"))
    }

    # Check required columns
    required_cols <- c("variable", "Value", "Response", "lower", "upper")
    if (!all(required_cols %in% colnames(df_time))) {
      stop(paste("res_cur in model", i, "must contain columns:",
                 paste(required_cols, collapse = ", ")))
    }

    df_time <- df_time |>
      dplyr::mutate(Time = factor(time_name, levels = time_names))

    df_list[[i]] <- df_time
  }

  df_all <- dplyr::bind_rows(df_list)

  # Get all unique explanatory variables
  variables <- unique(df_all$variable)
  n_vars <- length(variables)

  # Automatically calculate number of rows
  if (is.null(nrow)) {
    nrow <- ceiling(n_vars / ncol)
  }

  # Create color mapping
  time_colors <- stats::setNames(colors[1:n_times], time_names)

  # ---- Create individual plots for each variable ----
  plot_list <- list()

  for (i in seq_along(variables)) {
    var_name <- variables[i]

    # Filter data for current variable
    df_var <- df_all |> dplyr::filter(variable == var_name)

    # Calculate x-axis range
    x_min <- min(df_var$Value, na.rm = TRUE)
    x_max <- max(df_var$Value, na.rm = TRUE)
    x_range <- c(x_min, x_max)

    # Calculate y-axis range (considering confidence intervals)
    y_min <- min(c(df_var$lower, df_var$Response), na.rm = TRUE)
    y_max <- max(c(df_var$upper, df_var$Response), na.rm = TRUE)
    y_padding <- (y_max - y_min) * 0.05
    y_range <- c(y_min - y_padding, y_max + y_padding)

    # Create plot for single variable
    p <- ggplot2::ggplot(df_var, ggplot2::aes(x = Value, y = Response,
                                              fill = Time, color = Time)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                           alpha = 0.2, color = NA) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::coord_cartesian(xlim = x_range, ylim = y_range) +
      ggplot2::scale_fill_manual(values = time_colors) +
      ggplot2::scale_color_manual(values = time_colors) +
      ggplot2::labs(title = paste0(letters[i], ")"),
                    x = var_name,
                    y = y_label) +
      ggplot2::theme_minimal() +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0, face = "bold", size = 10),
        legend.position = "none",
        axis.title = ggplot2::element_text(size = 8),
        axis.text = ggplot2::element_text(size = 7),
        plot.margin = ggplot2::margin(5, 5, 5, 5)
      )

    plot_list[[i]] <- p
  }

  # ---- Helper: create vertical stacked legend ----
  create_vertical_legend <- function(time_names, colors) {
    n_times <- length(time_names)

    vertical_data <- data.frame()
    for (i in 1:n_times) {
      time_data <- data.frame(
        Time = factor(time_names[i], levels = time_names),
        x = seq(0, 2, length.out = 20),
        y = (n_times - i + 1) + 0.05 * sin(seq(0, 2, length.out = 20))
      )
      vertical_data <- rbind(vertical_data, time_data)
    }

    p <- ggplot2::ggplot(vertical_data, ggplot2::aes(x = x, y = y,
                                                     color = Time, fill = Time)) +
      ggplot2::geom_line(linewidth = 1, ggplot2::aes(group = Time)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = y - 0.1, ymax = y + 0.1, group = Time),
                           alpha = 0.2, color = NA) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_x_continuous(limits = c(0, 3.5), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0.5, n_times + 0.5)) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "none",
        plot.margin = ggplot2::margin(5, 30, 5, 5)
      )

    # Add text labels
    for (i in 1:n_times) {
      p <- p + ggplot2::annotate("text",
                                 x = 2.2,
                                 y = n_times - i + 1,
                                 label = time_names[i],
                                 hjust = 0,
                                 size = 3.5,
                                 color = "black")
    }
    return(p)
  }

  # ---- Create legend panel based on legend_position ----
  if (legend_position == "vertical") {
    # Vertical stacked legend
    legend_panel <- create_vertical_legend(time_names, time_colors)

    # Combine plots: first n_vars positions for plots, last position for legend
    total_plots <- nrow * ncol
    plots_for_combine <- plot_list

    # If total cells > number of variables, fill with empty plots
    if (total_plots > n_vars) {
      for (j in (n_vars + 1):total_plots) {
        if (j == total_plots) {  # Last cell for legend
          plots_for_combine[[j]] <- legend_panel
        } else {
          plots_for_combine[[j]] <- ggplot2::ggplot() + ggplot2::theme_void()
        }
      }
    } else if (total_plots == n_vars) {
      warning("No space for legend. Consider increasing ncol or nrow.")
      plots_for_combine[[n_vars]] <- legend_panel
    }

    combined_plot <- patchwork::wrap_plots(plots_for_combine, nrow = nrow, ncol = ncol)

  } else { # bottom legend
    combined_plot <- patchwork::wrap_plots(plot_list, nrow = nrow, ncol = ncol)
    combined_plot <- combined_plot &
      ggplot2::scale_fill_manual(values = time_colors) &
      ggplot2::scale_color_manual(values = time_colors)

    # Create a standalone legend
    legend_data <- data.frame()
    for (i in seq_along(time_names)) {
      x_seq <- seq(0, 1, length.out = 20)
      y_seq <- i + 0.1 * sin(seq(0, 2 * pi, length.out = 20))
      temp <- data.frame(
        Time = factor(time_names[i], levels = time_names),
        x = x_seq,
        y = y_seq,
        lower = y_seq - 0.15,
        upper = y_seq + 0.15
      )
      legend_data <- rbind(legend_data, temp)
    }

    legend_plot <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y,
                                                             color = Time, fill = Time)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                           alpha = 0.2, color = NA) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_color_manual(values = time_colors, name = "Time Period") +
      ggplot2::scale_fill_manual(values = time_colors, name = "Time Period") +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "center",
        legend.title = ggplot2::element_text(hjust = 0.5, size = 10),
        legend.text = ggplot2::element_text(hjust = 0.5, size = 9),
        legend.key.width = grid::unit(1, "cm"),
        legend.key.height = grid::unit(0.5, "cm"),
        legend.spacing.x = grid::unit(0.3, "cm"),
        plot.margin = ggplot2::margin(0, 0, 0, 0)
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          title = "Time Period",
          nrow = 1,
          override.aes = list(linewidth = 1)
        ),
        fill = ggplot2::guide_legend(
          title = "Time Period",
          nrow = 1,
          override.aes = list(alpha = 0.2)
        )
      )

    legend_grob <- ggplot2::ggplotGrob(legend_plot)
    legend_index <- which(sapply(legend_grob$grobs, function(x) x$name) == "guide-box")
    if (length(legend_index) > 0) {
      legend_only <- legend_grob$grobs[[legend_index[1]]]
    } else {
      legend_only <- legend_grob
    }

    n_items <- length(time_names)
    legend_width <- min(0.8, n_items * 0.12)
    margin_width <- (1 - legend_width) / 2

    combined_plot <- cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot, x = 0, y = 0.12, width = 1, height = 0.88) +
      cowplot::draw_grob(legend_only,
                         x = margin_width,
                         y = 0.02,
                         width = legend_width,
                         height = 0.08)
  }

  # ---- Save plot (if save_path is provided) ----
  if (!is.null(save_path)) {
    grDevices::tiff(save_path, width = width, height = height,
                    units = "in", res = 300, compression = "lzw")
    print(combined_plot)
    grDevices::dev.off()
  }

  return(combined_plot)
}
