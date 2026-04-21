#' Plot variable importance boxplot (applicable for any number of time periods)
#'
#' @param model_list Named list. Contains model results for each time period,
#'   e.g., `list(Quarter1 = sdm_q1, Quarter2 = sdm_q2, ...)`. Each element must
#'   contain a data frame `var_imp` with columns `Var` (variable name) and `Im` (importance).
#' @param time_names Character vector. Names of time periods, e.g.,
#'   `c("1st Quarter", "2nd Quarter", "3rd Quarter", "4th Quarter")`. Default is
#'   `paste0("Time ", 1:n_times)`.
#' @param colors Optional character vector of colors, length equal to number of
#'   time periods. Defaults to `RColorBrewer::brewer.pal(3, "Set3")` (or `rainbow`
#'   if more than 12).
#' @param title Character. Plot title. Default is `NULL`.
#' @param x_label Character. x‑axis label. Default: `"Importance"`.
#' @param y_label Character. y‑axis label. Default: `NULL` (no label).
#' @param width Numeric. Plot width in inches for saving. Default: `5`.
#' @param height Numeric. Plot height in inches for saving. Default: `6`.
#' @param save_path Character. Path to save the plot (TIFF format). If `NULL`,
#'   the plot is not saved.
#'
#' @return A `ggplot` object of the variable importance boxplot.
#' @export
#'
#' @importFrom dplyr mutate bind_rows
#' @importFrom ggplot2 ggplot aes geom_boxplot position_dodge labs scale_fill_manual
#'   theme_minimal theme_classic theme element_text
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices rainbow tiff dev.off
#'
#' @examples
#' \donttest{
#' # Prepare model list (four quarters) containing model object from TB_sdm()
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
#' # Plot variable importance
#' imp_plot <- TB_varimp(model_list,
#'                       time_names = time_names,
#'                       save_path = "./variable_importance.tif")
#' }
TB_varimp <- function(model_list,
                      time_names = NULL,
                      colors = NULL,
                      title = NULL,
                      x_label = "Importance",
                      y_label = NULL,
                      width = 5,
                      height = 6,
                      save_path = NULL) {

  # ---- Input validation ----
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

  # ---- Combine data from all time periods ----
  df_list <- list()
  for (i in 1:n_times) {
    model <- model_list[[i]]
    time_name <- time_names[i]

    # Extract variable importance data
    if (!is.null(model$var_imp)) {
      df_time <- model$var_imp
    } else {
      stop(paste("Model at temporal bin", i, "does not contain var_imp"))
    }

    # Ensure Var and Im columns exist
    if (!all(c("Var", "Im") %in% colnames(df_time))) {
      stop(paste("var_imp in model", i, "must contain 'Var' and 'Im' columns"))
    }

    df_time <- df_time |>
      dplyr::mutate(Time = factor(time_name, levels = time_names))

    df_list[[i]] <- df_time
  }

  df_all <- dplyr::bind_rows(df_list)

  # Create color mapping
  time_colors <- stats::setNames(colors[1:n_times], time_names)

  # ---- Create boxplot ----
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = Im, y = Var, fill = Time)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(0.8),
                          alpha = 0.7,
                          outlier.shape = 21,
                          outlier.size = 2,
                          outlier.stroke = 0.5) +
    ggplot2::labs(x = x_label,
                  y = y_label,
                  title = title,
                  fill = "Time Period") +
    ggplot2::scale_fill_manual(values = time_colors) +
    ggplot2::theme_minimal() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right")

  # ---- Save plot (if save_path is provided) ----
  if (!is.null(save_path)) {
    grDevices::tiff(save_path, width = width, height = height,
                    units = "in", res = 300, compression = "lzw")
    print(p)
    grDevices::dev.off()
  }

  return(p)
}
