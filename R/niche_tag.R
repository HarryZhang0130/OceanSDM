#' Estimate preference and tolerance ranges for a single species or environmental variable
#'
#' Based on multi-article binned data with different binning schemes.
#' Supports any continuous environmental variable (temperature, dissolved oxygen, etc.).
#'
#' @param data_list List. A list where each element is a data frame from one article,
#' containing columns: \code{bin_lower}, \code{bin_upper}, \code{percentage}, \code{article_weight}.
#' @param global_bins Data frame. A data frame with columns: \code{bin_lower}, \code{bin_upper}
#' defining standardized bins.
#' @param samples_per_article Numeric. Number of simulated samples per article for
#' standardization (default: 10000).
#' @param pref_prob Numeric. Proportion of area on each side of the mode to include in the
#' preference range. Default is 0.5 (50% of each side's area).
#' @param toler_prob Numeric. Proportion of area on each side of the mode to include in the
#' tolerance range. Default is 0.95 (95% of each side's area).
#' @param bw_method Character. Bandwidth selection method for kernel density estimation.
#' Options are:
#' \itemize{
#'   \item \code{"nrd0"}: Silverman's rule-of-thumb (default in \code{stats::density}).
#'         Fast and stable, suitable for unimodal, near-normal distributions, but may
#'         oversmooth multimodal or skewed distributions.
#'   \item \code{"nrd"}: Silverman's rule-of-thumb (variant). Slightly different
#'         scaling factor, similar properties to \code{"nrd0"}.
#'   \item \code{"SJ"}: Sheather–Jones direct plug-in estimator. Adaptive to
#'         data structure, preferred for complex/multimodal distributions, but
#'         computationally more intensive.
#'   \item \code{"bcv"}: Biased cross-validation. Tends to produce smaller bandwidths,
#'         useful for fine-scale features, but can be unstable.
#'   \item \code{"ucv"}: Unbiased cross-validation. Tends to produce larger bandwidths,
#'         more stable than \code{"bcv"} but computationally intensive.
#' }
#' @param variable_name Character. Name of the variable (e.g., \code{"Temperature"},
#' \code{"Dissolved Oxygen"}) for plotting (default: \code{"Environmental variable"}).
#' @param unit Character. Unit of the variable (e.g., \code{"°C"}, \code{"mg/L"}) for plotting
#' @param seed Integer. Random seed for reproducibility (default: 123).
#' @param width Numeric. Plot width in inches (default: 6).
#' @param height Numeric. Plot height in inches (default: 5).
#' @param save_path Character. Directory path to save plot and summary CSV.
#' If \code{NULL}, nothing is saved (default: \code{NULL}).
#'
#' @return A list of class \code{niche_range} containing:
#' \item{standardized_data}{Standardized binned data with weighted percentages}
#' \item{summary}{Data frame with optimal value, preference range, tolerance range}
#' \item{simulated_values}{Vector of simulated values (for KDE)}
#' \item{density_df}{Kernel density estimation results}
#' \item{sample_size}{Total number of simulated samples for KDE}
#' \item{total_individuals}{Total number of original individuals across articles}
#' \item{article_weights}{Normalized article weights}
#' \item{article_standardized_data}{Standardized data per article}
#' \item{variable_name}{User-provided variable name}
#' \item{unit}{User-provided unit}
#' \item{pref_prob_used}{Proportion used for preference range on each side}
#' \item{toler_prob_used}{Proportion used for tolerance range on each side}
#' \item{bw_method}{Bandwidth selection method used}
#'
#' @examples
#' \donttest{
#' # Example data: articles with temperature percentage distributions
#' article1_df <- data.frame(
#'   bin_lower = c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39),
#'   bin_upper = c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42),
#'   percentage = c(0, 0, 0, 0, 0, 0, 0.003, 0.11, 0.29, 0.37, 0.16, 0.06, 0.004, 0.003),
#'   article_weight = 50
#' )
#'
#' article2_df <- data.frame(
#'   bin_lower = c(10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30),
#'   bin_upper = c(12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5),
#'   percentage = c(0.008856, 0.008856, 0.006642, 0.004428, 0.004428, 0.046494,
#'                  0.294465, 0.593358, 0.028782),
#'   article_weight = 8
#' )
#'
#' data_list <- list(article1_df, article2_df)
#' global_bins <- article1_df[, c("bin_lower", "bin_upper")]
#'
#' # Use Sheather-Jones bandwidth (default)
#' temp_fit <- niche_tag(
#'   data_list = data_list,
#'   global_bins = global_bins,
#'   bw_method = "SJ",
#'   pref_prob = 0.5,
#'   toler_prob = 0.95,
#'   variable_name = "Temperature",
#'   unit = "°C"
#' )
#'
#' # Use Silverman's rule-of-thumb for faster computation
#' temp_fit_fast <- niche_tag(
#'   data_list = data_list,
#'   global_bins = global_bins,
#'   bw_method = "nrd0",
#'   pref_prob = 0.5,
#'   toler_prob = 0.95
#' )
#' }
#' @export
niche_tag <- function(
    data_list,
    global_bins,
    samples_per_article = 10000,
    pref_prob = 0.5,
    toler_prob = 0.95,
    bw_method = c( "nrd0", "nrd", "SJ", "bcv", "ucv"),
    variable_name = "Environmental variable",
    unit = NULL,
    seed = 123,
    width = 6,
    height = 5,
    save_path = NULL
) {

  # Match bandwidth method argument
  bw_method <- match.arg(bw_method)

  base::set.seed(seed)

  # Validate probabilities
  if (pref_prob <= 0 || pref_prob >= 1) {
    stop("pref_prob must be between 0 and 1 (exclusive).")
  }
  if (toler_prob <= 0 || toler_prob >= 1) {
    stop("toler_prob must be between 0 and 1 (exclusive).")
  }
  if (pref_prob >= toler_prob) {
    warning("pref_prob (", pref_prob, ") is greater than or equal to toler_prob (",
            toler_prob, "). This is ecologically unusual but allowed.")
  }

  # Create save directory if needed
  if (!is.null(save_path) && !dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  # Check global_bins format
  if (!all(c("bin_lower", "bin_upper") %in% colnames(global_bins))) {
    stop("global_bins must contain columns: bin_lower, bin_upper")
  }

  # Add bin_midpoint and bin_width to global_bins
  global_bins <- global_bins |>
    dplyr::mutate(
      bin_midpoint = (bin_lower + bin_upper) / 2,
      bin_width = bin_upper - bin_lower
    )

  # =========================================================================
  # STEP 1: Simulate raw data per article and re-bin to global standardized bins
  # =========================================================================
  standardized_articles <- list()
  article_weights_raw <- c()

  for (i in seq_along(data_list)) {
    article_data <- data_list[[i]]
    cat("\nProcessing Article", i, "...\n")

    required_cols <- c("bin_lower", "bin_upper", "percentage", "article_weight")
    if (!all(required_cols %in% colnames(article_data))) {
      stop(paste("Article", i, "missing required columns:",
                 paste(setdiff(required_cols, colnames(article_data)), collapse = ", ")))
    }

    weight <- unique(article_data$article_weight)
    if (length(weight) > 1) {
      warning(paste("Article", i, "has multiple weights, using the first one:", weight[1]))
      weight <- weight[1]
    }
    article_weights_raw <- c(article_weights_raw, weight)

    total_pct <- sum(article_data$percentage)
    if (abs(total_pct - 1) > 0.01) {
      warning(paste("Article", i, "percentage sum is", round(total_pct, 3), "- normalizing"))
      article_data$percentage <- article_data$percentage / total_pct
    }

    simulated_raw <- numeric()
    for (j in 1:nrow(article_data)) {
      n_in_bin <- round(article_data$percentage[j] * samples_per_article)
      if (n_in_bin > 0) {
        temp_points <- stats::runif(
          n = n_in_bin,
          min = article_data$bin_lower[j],
          max = article_data$bin_upper[j]
        )
        simulated_raw <- c(simulated_raw, temp_points)
      }
    }

    cat("  Simulated", length(simulated_raw), "raw values from original bins\n")

    if (length(simulated_raw) == 0) {
      warning(paste("Article", i, "generated no simulated data - skipping"))
      next
    }

    standardized_pct <- numeric(nrow(global_bins))
    for (k in 1:nrow(global_bins)) {
      in_bin <- sum(simulated_raw >= global_bins$bin_lower[k] &
                      simulated_raw < global_bins$bin_upper[k])
      standardized_pct[k] <- in_bin / length(simulated_raw)
    }

    if (abs(sum(standardized_pct) - 1) > 0.01) {
      warning(paste("Article", i, "standardized percentages sum to",
                    round(sum(standardized_pct), 3), "- adjusting"))
      standardized_pct <- standardized_pct / sum(standardized_pct)
    }

    standardized_articles[[i]] <- data.frame(
      article = i,
      bin_lower = global_bins$bin_lower,
      bin_upper = global_bins$bin_upper,
      bin_midpoint = global_bins$bin_midpoint,
      percentage = standardized_pct,
      density = standardized_pct / global_bins$bin_width
    )
  }

  # =========================================================================
  # STEP 2: Combine articles using weighted average (weights = sample sizes)
  # =========================================================================
  cat("\n=== Combining articles with weighted average ===\n")
  total_individuals <- sum(article_weights_raw)
  cat("Total number of individuals =", total_individuals, "\n")

  article_weights_norm <- article_weights_raw / sum(article_weights_raw)
  cat("Normalized article weights:", paste(round(article_weights_norm, 3), collapse = ", "), "\n")

  combined_data <- global_bins[, c("bin_lower", "bin_upper", "bin_midpoint", "bin_width")]
  combined_data$percentage <- 0
  combined_data$density <- 0

  for (i in seq_along(standardized_articles)) {
    combined_data$percentage <- combined_data$percentage +
      standardized_articles[[i]]$percentage * article_weights_norm[i]
  }

  total_pct <- sum(combined_data$percentage)
  if (abs(total_pct - 1) > 0.01) {
    cat("Combined percentages sum to", round(total_pct, 3), "- normalizing\n")
    combined_data$percentage <- combined_data$percentage / total_pct
  }
  combined_data$density <- combined_data$percentage / combined_data$bin_width

  cat("\nFinal standardized binned data (weighted average):\n")
  print(combined_data[, c("bin_lower", "bin_upper", "percentage", "density")])

  # =========================================================================
  # STEP 3: Generate final simulated data from standardized weighted bins
  # =========================================================================
  cat("\n=== Generating final simulated data from standardized bins ===\n")
  final_samples <- 50000
  final_simulated <- numeric()

  for (k in 1:nrow(combined_data)) {
    n_in_bin <- round(combined_data$percentage[k] * final_samples)
    if (n_in_bin > 0) {
      values <- stats::runif(
        n = n_in_bin,
        min = combined_data$bin_lower[k],
        max = combined_data$bin_upper[k]
      )
      final_simulated <- c(final_simulated, values)
    }
  }

  cat("Generated", length(final_simulated), "simulated values from standardized bins\n")

  # =========================================================================
  # STEP 4: Kernel density estimation and threshold calculation
  # =========================================================================
  # Use user-specified bandwidth method
  density_est <- stats::density(final_simulated, kernel = "gaussian", bw = bw_method)
  density_df <- data.frame(x = density_est$x, y = density_est$y)

  cat("\nBandwidth method:", bw_method, "| Bandwidth value:", round(density_est$bw, 4), "\n")

  # ---- Helper function: compute CDF from KDE density ----
  compute_cdf_from_kde <- function(density_df) {
    x_vals <- density_df$x
    y_vals <- density_df$y
    n <- length(x_vals)

    # Normalize density to integrate to 1
    integral <- 0
    for (i in 2:n) {
      integral <- integral + (y_vals[i] + y_vals[i-1]) / 2 * (x_vals[i] - x_vals[i-1])
    }
    y_vals_norm <- y_vals / integral

    # Compute cumulative distribution function (CDF)
    cdf <- numeric(n)
    cdf[1] <- 0
    for (i in 2:n) {
      cdf[i] <- cdf[i-1] + (y_vals_norm[i] + y_vals_norm[i-1]) / 2 * (x_vals[i] - x_vals[i-1])
    }
    cdf <- cdf / max(cdf)

    return(list(x = x_vals, cdf = cdf, y_vals_norm = y_vals_norm))
  }

  # ---- Helper function: compute probability mass between two values ----
  compute_prob_between <- function(density_df, x1, x2) {
    cdf_result <- compute_cdf_from_kde(density_df)
    cdf_func <- approxfun(cdf_result$x, cdf_result$cdf, rule = 2)
    prob <- cdf_func(x2) - cdf_func(x1)
    return(max(0, min(1, prob)))
  }

  # ---- Helper function: compute total area on each side of the mode ----
  compute_side_areas <- function(density_df, mode_value) {
    x_vals <- density_df$x
    y_vals <- density_df$y
    n <- length(x_vals)

    # Normalize density
    integral <- 0
    for (i in 2:n) {
      integral <- integral + (y_vals[i] + y_vals[i-1]) / 2 * (x_vals[i] - x_vals[i-1])
    }
    y_vals_norm <- y_vals / integral

    # Find mode index
    mode_idx <- which.min(abs(x_vals - mode_value))

    # Compute left side area (from min(x) to mode)
    left_area <- 0
    for (i in 2:mode_idx) {
      left_area <- left_area + (y_vals_norm[i] + y_vals_norm[i-1]) / 2 * (x_vals[i] - x_vals[i-1])
    }

    # Compute right side area (from mode to max(x))
    right_area <- 0
    for (i in (mode_idx + 1):n) {
      right_area <- right_area + (y_vals_norm[i] + y_vals_norm[i-1]) / 2 * (x_vals[i] - x_vals[i-1])
    }

    return(list(left_area = left_area, right_area = right_area))
  }

  # ---- Helper function: find threshold on one side for a given proportion of that side's area ----
  find_side_threshold <- function(density_df, mode_value, side, target_prob) {
    x_vals <- density_df$x
    y_vals <- density_df$y
    n <- length(x_vals)

    # Normalize density
    integral <- 0
    for (i in 2:n) {
      integral <- integral + (y_vals[i] + y_vals[i-1]) / 2 * (x_vals[i] - x_vals[i-1])
    }
    y_vals_norm <- y_vals / integral

    mode_idx <- which.min(abs(x_vals - mode_value))

    if (side == "left") {
      # Compute total left side area
      total_side_area <- 0
      for (i in 2:mode_idx) {
        total_side_area <- total_side_area + (y_vals_norm[i] + y_vals_norm[i-1]) / 2 * (x_vals[i] - x_vals[i-1])
      }

      target_area <- total_side_area * target_prob

      # Find threshold by expanding leftwards from mode
      threshold <- mode_value
      current_area <- 0
      step_size <- (mode_value - min(x_vals)) / 200

      while (current_area < target_area && threshold > min(x_vals)) {
        threshold <- threshold - step_size
        threshold <- max(threshold, min(x_vals))
        current_area <- compute_prob_between(density_df, threshold, mode_value)
      }

      actual_area <- compute_prob_between(density_df, threshold, mode_value)
      total_area <- total_side_area

    } else { # right side
      # Compute total right side area
      total_side_area <- 0
      for (i in (mode_idx + 1):n) {
        total_side_area <- total_side_area + (y_vals_norm[i] + y_vals_norm[i-1]) / 2 * (x_vals[i] - x_vals[i-1])
      }

      target_area <- total_side_area * target_prob

      # Find threshold by expanding rightwards from mode
      threshold <- mode_value
      current_area <- 0
      step_size <- (max(x_vals) - mode_value) / 200

      while (current_area < target_area && threshold < max(x_vals)) {
        threshold <- threshold + step_size
        threshold <- min(threshold, max(x_vals))
        current_area <- compute_prob_between(density_df, mode_value, threshold)
      }

      actual_area <- compute_prob_between(density_df, mode_value, threshold)
      total_area <- total_side_area
    }

    return(list(
      threshold = threshold,
      actual_area = actual_area,
      total_side_area = total_side_area,
      target_prob = target_prob,
      target_area = target_area,
      proportion_achieved = actual_area / total_side_area
    ))
  }

  # ---- Helper function: compute full asymmetric interval ----
  find_asymmetric_interval <- function(density_df, target_prob) {
    x_vals <- density_df$x
    y_vals <- density_df$y
    n <- length(x_vals)

    # Normalize density
    integral <- 0
    for (i in 2:n) {
      integral <- integral + (y_vals[i] + y_vals[i-1]) / 2 * (x_vals[i] - x_vals[i-1])
    }
    y_vals_norm <- y_vals / integral

    # Find mode
    mode_idx <- which.max(y_vals_norm)
    mode_value <- x_vals[mode_idx]

    # Compute left and right side areas
    side_areas <- compute_side_areas(density_df, mode_value)
    left_total <- side_areas$left_area
    right_total <- side_areas$right_area

    # Find left threshold: covers target_prob of left side area
    left_result <- find_side_threshold(density_df, mode_value, "left", target_prob)
    lower_bound <- left_result$threshold
    left_actual_prob <- left_result$proportion_achieved
    left_actual_area <- left_result$actual_area

    # Find right threshold: covers target_prob of right side area
    right_result <- find_side_threshold(density_df, mode_value, "right", target_prob)
    upper_bound <- right_result$threshold
    right_actual_prob <- right_result$proportion_achieved
    right_actual_area <- right_result$actual_area

    # Total actual area covered (as proportion of total density)
    total_area_left <- compute_prob_between(density_df, lower_bound, mode_value)
    total_area_right <- compute_prob_between(density_df, mode_value, upper_bound)
    total_actual_prob <- total_area_left + total_area_right

    return(list(
      lower = lower_bound,
      upper = upper_bound,
      mode = mode_value,
      left_total_area = left_total,
      right_total_area = right_total,
      left_actual_area = left_actual_area,
      right_actual_area = right_actual_area,
      left_actual_proportion = left_actual_prob,
      right_actual_proportion = right_actual_prob,
      total_actual_prob = total_actual_prob,
      target_prob = target_prob
    ))
  }

  # ---- Optimal value (mode) ----
  optimal_value <- density_df$x[which.max(density_df$y)]

  # ---- Preference Range (asymmetric interval) ----
  pref_interval <- find_asymmetric_interval(density_df, pref_prob)
  lower_pref <- pref_interval$lower
  upper_pref <- pref_interval$upper

  cat("\n=== Preference Range (", pref_prob*100, "% of each side's area) ===\n", sep = "")
  cat("  Mode value:", round(pref_interval$mode, 2), "\n")
  cat("  Left side total area:", round(pref_interval$left_total_area * 100, 1), "%\n")
  cat("    Left threshold:", round(lower_pref, 2), " (covers ",
      round(pref_interval$left_actual_proportion * 100, 1), "% of left side)\n", sep = "")
  cat("  Right side total area:", round(pref_interval$right_total_area * 100, 1), "%\n")
  cat("    Right threshold:", round(upper_pref, 2), " (covers ",
      round(pref_interval$right_actual_proportion * 100, 1), "% of right side)\n", sep = "")
  cat("  Total actual probability covered:", round(pref_interval$total_actual_prob * 100, 1), "%\n")

  # ---- Tolerance Range (asymmetric interval) ----
  toler_interval <- find_asymmetric_interval(density_df, toler_prob)
  lower_tol <- toler_interval$lower
  upper_tol <- toler_interval$upper

  cat("\n=== Tolerance Range (", toler_prob*100, "% of each side's area) ===\n", sep = "")
  cat("  Mode value:", round(toler_interval$mode, 2), "\n")
  cat("  Left side total area:", round(toler_interval$left_total_area * 100, 1), "%\n")
  cat("    Left threshold:", round(lower_tol, 2), " (covers ",
      round(toler_interval$left_actual_proportion * 100, 1), "% of left side)\n", sep = "")
  cat("  Right side total area:", round(toler_interval$right_total_area * 100, 1), "%\n")
  cat("    Right threshold:", round(upper_tol, 2), " (covers ",
      round(toler_interval$right_actual_proportion * 100, 1), "% of right side)\n", sep = "")
  cat("  Total actual probability covered:", round(toler_interval$total_actual_prob * 100, 1), "%\n")

  # =========================================================================
  # STEP 5: Create summary data frame
  # =========================================================================
  summary_df <- data.frame(
    parameter = c("optimal_value",
                  "lower_preferred", "upper_preferred",
                  "lower_tolerance", "upper_tolerance"),
    value = c(optimal_value, lower_pref, upper_pref, lower_tol, upper_tol),
    unit = rep(unit, 5)
  )

  # =========================================================================
  # STEP 6: Return results with S3 class
  # =========================================================================
  result <- list(
    standardized_data = combined_data,
    summary = summary_df,
    simulated_values = final_simulated,
    density_df = density_df,
    sample_size = length(final_simulated),
    total_individuals = total_individuals,
    article_weights = article_weights_norm,
    article_standardized_data = standardized_articles,
    variable_name = variable_name,
    unit = unit,
    pref_prob_used = pref_prob,
    toler_prob_used = toler_prob,
    bw_method = bw_method,
    bw_value = density_est$bw,
    pref_interval = pref_interval,
    toler_interval = toler_interval
  )

  class(result) <- "niche_range"

  # =========================================================================
  # S3 print method
  # =========================================================================
  print.niche_range <- function(x, ...) {
    cat("\n=== Species Preference and Tolerance Range Estimates for", x$variable_name, " ===\n")
    cat("Method: Standardized bins + weighted average across articles\n")
    cat("Bandwidth method:", x$bw_method, "| Bandwidth value:", round(x$bw_value, 4), "\n")
    cat("Range estimation: Asymmetric interval using mode as dividing line\n")
    cat("  Each side independently covers target proportion of its own area\n")
    cat("  Preference: covers", x$pref_prob_used * 100, "% of each side's area\n")
    cat("  Tolerance:  covers", x$toler_prob_used * 100, "% of each side's area\n")
    cat("Final sample size for KDE:", x$sample_size, "\n")
    cat("Total individuals across articles:", x$total_individuals, "\n")
    cat("Number of articles combined:", length(x$article_weights), "\n\n")

    cat("Optimal", x$variable_name, ":\n")
    cat("  Value:", round(x$summary$value[x$summary$parameter == "optimal_value"], 2), x$unit, "\n\n")

    cat("Preferred Range:\n")
    cat("  Lower:", round(x$summary$value[x$summary$parameter == "lower_preferred"], 2), x$unit, "\n")
    cat("  Upper:", round(x$summary$value[x$summary$parameter == "upper_preferred"], 2), x$unit, "\n\n")

    cat("Tolerance Range:\n")
    cat("  Lower:", round(x$summary$value[x$summary$parameter == "lower_tolerance"], 2), x$unit, "\n")
    cat("  Upper:", round(x$summary$value[x$summary$parameter == "upper_tolerance"], 2), x$unit, "\n")
  }

  # =========================================================================
  # S3 plot method
  # =========================================================================
  plot.niche_range <- function(x, ...) {
    optimal <- x$summary$value[x$summary$parameter == "optimal_value"]
    lower_pref <- x$summary$value[x$summary$parameter == "lower_preferred"]
    upper_pref <- x$summary$value[x$summary$parameter == "upper_preferred"]
    lower_tol <- x$summary$value[x$summary$parameter == "lower_tolerance"]
    upper_tol <- x$summary$value[x$summary$parameter == "upper_tolerance"]

    x_label <- paste(x$variable_name, "(", x$unit, ")", sep = "")

    toler_label <- paste("Tolerance (", x$toler_prob_used * 100, "% of each side)", sep = "")
    pref_label <- paste("Preference (", x$pref_prob_used * 100, "% of each side)", sep = "")

    p <- ggplot2::ggplot() +
      ggplot2::geom_col(data = x$standardized_data,
                        ggplot2::aes(x = bin_midpoint, y = density, width = bin_width),
                        fill = "lightblue", color = "black", alpha = 0.5) +
      ggplot2::geom_line(data = x$density_df, ggplot2::aes(x = x, y = y),
                         linewidth = 1.2, color = "darkblue") +
      # Tolerance range (green shaded)
      ggplot2::annotate("rect", xmin = lower_tol, xmax = upper_tol,
                        ymin = 0, ymax = Inf, fill = "green", alpha = 0.05) +
      ggplot2::geom_vline(xintercept = c(lower_tol, upper_tol),
                          linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
      # Preference range (orange shaded)
      ggplot2::annotate("rect", xmin = lower_pref, xmax = upper_pref,
                        ymin = 0, ymax = Inf, fill = "orange", alpha = 0.08) +
      ggplot2::geom_vline(xintercept = c(lower_pref, upper_pref),
                          linetype = "dashed", color = "orange", linewidth = 1) +
      # Optimal value
      ggplot2::geom_vline(xintercept = optimal, linetype = "dashed",
                          color = "red", linewidth = 1) +
      ggplot2::annotate("text", x = optimal, y = max(x$density_df$y) * 0.9,
                        label = paste("Optimal:", round(optimal, 1), x$unit),
                        hjust = -0.1, color = "red", size = 4) +
      # Legend annotations
      ggplot2::annotate("text", x = max(x$density_df$x) * 0.55,
                        y = max(x$density_df$y) * 0.15,
                        label = toler_label,
                        color = "darkgreen", size = 3.5, hjust = 0) +
      ggplot2::annotate("text", x = max(x$density_df$x) * 0.55,
                        y = max(x$density_df$y) * 0.10,
                        label = pref_label,
                        color = "orange", size = 3.5, hjust = 0) +
      ggplot2::labs(title = paste("Species' Preference and Tolerance Range Estimates for", x$variable_name),
                    subtitle = paste("Based on", x$total_individuals, "tagged individuals from",
                                     length(x$article_weights), "articles with standardized bins"),
                    x = x_label,
                    y = "Density") +
      ggplot2::theme_minimal()

    return(p)
  }

  # Print and plot by default
  print.niche_range(result)
  print(plot.niche_range(result))

  # Save outputs if path provided
  if (!is.null(save_path)) {
    grDevices::tiff(paste0(save_path, "/niche_fit.tif"),
                    width = width, height = height, units = "in", res = 300, compression = "lzw")
    print(plot.niche_range(result))
    grDevices::dev.off()

    utils::write.csv(summary_df[, c("parameter", "value")],
                     paste0(save_path, "/summary_niche_fit.csv"),
                     row.names = FALSE)
  }

  return(result)
}
