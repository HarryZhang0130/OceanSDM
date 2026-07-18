#' Plot variable importance as a grouped barplot across multiple time periods
#'
#' This function reads variable importance CSV files from subfolders of a root
#' directory, where each subfolder corresponds to a time period (e.g., Q1, Q2, etc.).
#' The CSV file must contain columns 'variables', 'corTest', and 'AUCtest' (or other
#' columns as specified by 'method'). Missing variables in any period are filled with 0.
#'
#' @param root_path Character. Path to the root directory containing time period subfolders.
#' @param time_type Character. Type of time periods: "quarter", "month", "season".
#'   Default "quarter". Ignored if `time_folders` is provided.
#' @param time_folders Optional character vector of folder names (e.g., c("Q1","Q2","Q3","Q4")).
#'   If provided, overrides `time_type`.
#' @param time_names Optional character vector of labels for the time periods
#'   (used in legend). Default uses folder names.
#' @param file_name Character. Name of the variable importance CSV file within each
#'   subfolder (e.g., "varimp.csv"). Default: "varimp.csv".
#' @param method Character. Which importance metric to use: either "corTest" or
#'   "AUCtest". Default "AUCtest".
#' @param colors Optional character vector of colors, length equal to number of
#'   time periods. If NULL (default), uses ggplot2's default discrete color palette.
#' @param title Character. Plot title. Default NULL.
#' @param x_label Character. x-axis label. Default: "Importance".
#' @param y_label Character. y-axis label. Default: "Variable".
#' @param group_label Character. Label for the temporal bins (e.g., Quarters, Seasons, Months).
#'   If NULL, automatically set based on `time_type`.
#' @param width Numeric. Plot width in inches for saving. Default: 8.
#' @param height Numeric. Plot height in inches for saving. Default: 6.
#' @param save_path Character. Path to save the plot (TIFF format). If NULL,
#'   the plot is not saved.
#'
#' @return A `ggplot` object of the grouped barplot.
#'
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot aes geom_col position_dodge labs
#'   scale_fill_manual scale_fill_hue theme_minimal theme_classic theme element_text
#' @importFrom grDevices tiff dev.off
#' @importFrom utils read.csv
#' @examples
#' \donttest{
#' # Name each temporal bin
#' time_names <- c("1st Quarter", "2nd Quarter", "3rd Quarter", "4th Quarter")
#' #  Plot variable importance
#' imp_plot <- TB_varimp(root_path = "F:/whaleshark_sdm/result_sdm",
#'                       time_type = "quarter",
#'                       file_name = "VarImportance_Table.csv",
#'                       time_names = time_names,
#'                       time_folders = NULL,
#'                       method = "AUCtest",
#'                       colors = c("1st Quarter" = "#FFB347",
#'                                  "2nd Quarter" = "#4CAF50",
#'                                  "3rd Quarter" = "#FF6B4A",
#'                                  "4th Quarter" = "#5D9BEC"),
#'                       title = NULL,
#'                       x_label = "Importance",
#'                       y_label = NULL,
#'                       group_label= "Quarters",
#'                       width = 8,
#'                       height = 6,
#'                       save_path = "F:/whaleshark_sdm/result_sdm/variable_importance.tif")
#' print(imp_plot)
#' }
#' @export
TB_varimp <- function(root_path,
                      time_type = "quarter",
                      time_folders = NULL,
                      time_names = NULL,
                      file_name = "varimp.csv",
                      method = "AUCtest",
                      colors = NULL,
                      title = NULL,
                      x_label = "Importance",
                      y_label = "Variable",
                      group_label = NULL,
                      width = 8,
                      height = 6,
                      save_path = NULL) {

  # ---- 1. Determine time periods (folders) ----
  if (!is.null(time_folders)) {
    periods <- time_folders
  } else {
    if (time_type == "quarter") {
      periods <- c("Q1", "Q2", "Q3", "Q4")
    } else if (time_type == "month") {
      periods <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    } else if (time_type == "season") {
      periods <- c("Sp", "Sm", "Am", "Wt") # Spring, Summer, Autumn, Winter
    } else {
      stop("time_type must be 'quarter', 'month', or 'season'")
    }
  }

  if (is.null(time_names)) {
    time_names <- periods
  } else {
    if (length(time_names) != length(periods)) {
      stop("length of time_names must equal number of periods")
    }
  }

  # Set default group_label if not provided
  if (is.null(group_label)) {
    if (time_type == "quarter") group_label <- "Quarters"
    else if (time_type == "month") group_label <- "Months"
    else if (time_type == "season") group_label <- "Seasons"
    else group_label <- "Time Periods"
  }

  n_times <- length(periods)

  # Validate method
  method <- match.arg(method, c("corTest", "AUCtest"))

  # ---- 2. Read all variable importance files ----
  var_imp_list <- list()
  all_vars <- c()
  for (i in 1:n_times) {
    period <- periods[i]
    file_path <- file.path(root_path, period, file_name)
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path, " - skipping period ", period)
      next
    }
    df <- utils::read.csv(file_path, stringsAsFactors = FALSE)
    required_cols <- c("variables", "corTest", "AUCtest")
    if (!all(required_cols %in% colnames(df))) {
      warning("File ", file_path, " missing required columns; skipping.")
      next
    }
    df <- df[, c("variables", "corTest", "AUCtest")]
    var_imp_list[[period]] <- df
    all_vars <- union(all_vars, df$variables)
  }

  if (length(var_imp_list) == 0) {
    stop("No valid variable importance files found.")
  }

  # ---- 3. Fill missing variables with 0 for each period ----
  df_all <- data.frame()
  for (i in 1:n_times) {
    period <- periods[i]
    if (is.null(var_imp_list[[period]])) {
      imp_df <- data.frame(variables = all_vars, corTest = NA, AUCtest = NA)
    } else {
      imp_df <- var_imp_list[[period]]
    }
    complete_df <- data.frame(variables = all_vars, stringsAsFactors = FALSE)
    complete_df <- merge(complete_df, imp_df, by = "variables", all.x = TRUE)
    complete_df[is.na(complete_df[[method]]), method] <- 0

    df_plot <- data.frame(
      Variable = complete_df$variables,
      Importance = complete_df[[method]],
      Time = factor(time_names[i], levels = time_names),
      stringsAsFactors = FALSE
    )
    df_all <- rbind(df_all, df_plot)
  }

  # ---- 4. Create grouped barplot (horizontal bars) ----
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = Importance, y = Variable, fill = Time)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.9), width = 0.8) +
    ggplot2::labs(x = x_label, y = y_label, fill = group_label, title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      legend.position = "right"
    )

  # Apply color scale: user-provided colors or ggplot2 default hue
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors[1:n_times])
  } else {
    p <- p + ggplot2::scale_fill_hue()  # ggplot2's default discrete color palette
  }

  # ---- 5. Save plot if requested ----
  if (!is.null(save_path)) {
    grDevices::tiff(save_path, width = width, height = height,
                    units = "in", res = 300, compression = "lzw")
    print(p)
    grDevices::dev.off()
  }

  return(p)
}
