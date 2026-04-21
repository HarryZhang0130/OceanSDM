#' Estimate preference and tolerance ranges for a single species or environmental variable
#'
#' Based on multi-article binned data with different binning schemes.
#' Supports any continuous environmental variable (temperature, dissolved oxygen, etc.).
#'
#' @param data_list List. A list where each element is a data frame from one article,
#'   containing columns: `bin_lower`, `bin_upper`, `percentage`, `article_weight`.
#' @param global_bins Data frame. A data frame with columns: `bin_lower`, `bin_upper`
#'   defining standardized bins.
#' @param samples_per_article Numeric. Number of simulated samples per article for
#'   standardization (default: 10000).
#' @param toler_quantiles Numeric vector of length 2. Quantiles for tolerance range
#'   (default: c(0.01, 0.99)).
#' @param pref_method Character. Method for preference range: `"FWHM"` (Full Width at
#'   Half Maximum) or `"quantile"` (default: `"FWHM"`).
#' @param pref_quantiles Numeric vector of length 2. Quantiles for preference range
#'   when `pref_method = "quantile"` (default: c(0.25, 0.75)).
#' @param variable_name Character. Name of the variable (e.g., `"Temperature"`,
#'   `"Dissolved Oxygen"`) for plotting (default: `"Environmental variable"`).
#' @param unit Character. Unit of the variable (e.g., `"°C"`, `"mg/L"`) for plotting
#'   (default: `"℃"`).
#' @param seed Integer. Random seed for reproducibility (default: 123).
#' @param width Numeric. Plot width in inches (default: 6).
#' @param height Numeric. Plot height in inches (default: 5).
#' @param save_path Character. Directory path to save plot and summary CSV.
#'   If `NULL`, nothing is saved (default: `NULL`).
#'
#' @return A list of class `niche_range` containing:
#'   \item{standardized_data}{Standardized binned data with weighted percentages}
#'   \item{summary}{Data frame with optimal value, preference range, tolerance range}
#'   \item{simulated_values}{Vector of simulated values (for KDE)}
#'   \item{density_df}{Kernel density estimation results}
#'   \item{sample_size}{Total number of simulated samples for KDE}
#'   \item{total_individuals}{Total number of original individuals across articles}
#'   \item{article_weights}{Normalized article weights}
#'   \item{article_standardized_data}{Standardized data per article}
#'   \item{variable_name}{User-provided variable name}
#'   \item{unit}{User-provided unit}
#'   \item{pref_method}{Method used for preference range}
#'   \item{pref_quantiles}{Preference quantiles (if applicable)}
#'
#' @examples
#' \donttest{
#' # Example data: 7 articles with temperature percentage distributions
#' # Article 1 (largest sample size, used as reference for global bins)
#' article1_df <- data.frame(
#'   bin_lower = c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39),
#'   bin_upper = c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42),
#'   percentage = c(0, 0, 0, 0, 0, 0, 0.003, 0.11, 0.29, 0.37, 0.16, 0.06, 0.004, 0.003),
#'   article_weight = 50
#' )
#'
#' # Article 2
#' article2_df <- data.frame(
#'   bin_lower = c(10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30),
#'   bin_upper = c(12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5),
#'   percentage = c(0.008856, 0.008856, 0.006642, 0.004428, 0.004428, 0.046494,
#'                  0.294465, 0.593358, 0.028782),
#'   article_weight = 8
#' )
#'
#' # Article 3
#' article3_df <- data.frame(
#'   bin_lower = c(10, 12, 14, 16, 18, 20, 22, 24, 26, 28),
#'   bin_upper = c(12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
#'   percentage = c(0.0065445, 0.0302054, 0.0436468, 0.0436468, 0.0436468,
#'                  0.06207209, 0.1258558, 0.2550846, 0.2601188, 0.1291784),
#'   article_weight = 1
#' )
#'
#' # Article 4
#' article4_df <- data.frame(
#'   bin_lower = c(0, 2, 5, 10, 15, 17.5, 20, 22, 25, 27, 29, 31),
#'   bin_upper = c(2, 5, 10, 15, 17.5, 20, 22, 25, 27, 29, 31, 33),
#'   percentage = c(0.00039048, 0.00002092, 0.00055783, 0.02545097, 0.04554678,
#'                  0.08029258, 0.1846834, 0.4400787, 0.1277708, 0.08686102,
#'                  0.00594088, 0.00240564),
#'   article_weight = 15
#' )
#'
#' # Article 5
#' article5_df <- data.frame(
#'   bin_lower = c(0, 5, 10, 15, 18, 20, 22, 24, 26, 28, 30, 32, 34),
#'   bin_upper = c(5, 10, 15, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36),
#'   percentage = c(0.00002896, 0.00004234, 0.00013009, 0.00012676, 0.00011109,
#'                  0.00018206, 0.00055804, 0.00443287, 0.0041829, 0.00020489,
#'                  0, 0, 0),
#'   article_weight = 27
#' )
#'
#' # Article 6
#' article6_df <- data.frame(
#'   bin_lower = c(4, 8, 12, 16, 20, 22, 24, 26, 28),
#'   bin_upper = c(8, 12, 16, 20, 22, 24, 26, 28, 30),
#'   percentage = c(0.000188, 0.000075, 0.000063, 0.000088, 0.000151,
#'                  0.001217, 0.003111, 0.004027, 0.001116),
#'   article_weight = 17
#' )
#'
#' # Article 7
#' article7_df <- data.frame(
#'   bin_lower = c(10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5),
#'   bin_upper = c(12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35),
#'   percentage = c(0, 0.00029412, 0.00011765, 0.00011765, 0.00019608,
#'                  0.00060785, 0.00372549, 0.0037647, 0.00115686, 0.00001961),
#'   article_weight = 17
#' )
#'
#' # Combine all articles into a list
#' data_list <- list(article1_df, article2_df, article3_df, article4_df,
#'                   article5_df, article6_df, article7_df)
#'
#' # Define unified global binning scheme (based on Article 1)
#' global_bins <- article1_df[, c("bin_lower", "bin_upper")]
#'
#' # Run niche estimation
#' temp_fit <- niche_tag(
#'   data_list = data_list,
#'   global_bins = global_bins,
#'   samples_per_article = 10000,
#'   toler_quantiles = c(0.01, 0.99),
#'   pref_method = "FWHM",
#'   pref_quantiles = NULL,
#'   variable_name = "Temperature",
#'   unit = "°C",
#'   seed = 123,
#'   width = 6,
#'   height = 5,
#'   save_path = NULL
#' )
#'
#' # View summary
#' print(temp_fit)
#' temp_fit$summary
#'
#' # Generate plot
#' plot(temp_fit)
#' }
#' @export
niche_tag <- function(
    data_list,
    global_bins,
    samples_per_article = 10000,
    toler_quantiles = c(0.01, 0.99),
    pref_method = c("FWHM", "quantile"),
    pref_quantiles = c(0.25, 0.75),
    variable_name = "Environmental variable",
    unit = "℃",
    seed = 123,
    width = 6,
    height = 5,
    save_path = NULL
    ) {
  # Match argument for pref_method
  pref_method <- match.arg(pref_method)

  base::set.seed(seed)

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

  # ---------------------------------------------------------------------------
  # STEP 1: Simulate raw data per article and re-bin to global standardized bins
  # ---------------------------------------------------------------------------
  standardized_articles <- list()
  article_weights_raw <- c()  # store original weights (sample sizes)

  for (i in seq_along(data_list)) {
    article_data <- data_list[[i]]
    cat("\nProcessing Article", i, "...\n")

    # Check required columns
    required_cols <- c("bin_lower", "bin_upper", "percentage", "article_weight")
    if (!all(required_cols %in% colnames(article_data))) {
      stop(paste("Article", i, "missing required columns:",
                 paste(setdiff(required_cols, colnames(article_data)), collapse = ", ")))
    }

    # Get article weight (should be single value)
    weight <- unique(article_data$article_weight)
    if (length(weight) > 1) {
      warning(paste("Article", i, "has multiple weights, using the first one:", weight[1]))
      weight <- weight[1]
    }
    article_weights_raw <- c(article_weights_raw, weight)

    # Normalize percentages within article (ensure sum = 1)
    total_pct <- sum(article_data$percentage)
    if (abs(total_pct - 1) > 0.01) {
      warning(paste("Article", i, "percentage sum is", round(total_pct, 3), "- normalizing"))
      article_data$percentage <- article_data$percentage / total_pct
    }

    # Simulate raw data from original bins
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

    # Rebin to global standardized bins
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

  # ---------------------------------------------------------------------------
  # STEP 2: Combine articles using weighted average (weights = sample sizes)
  # ---------------------------------------------------------------------------
  cat("\n=== Combining articles with weighted average ===\n")
  total_individuals <- sum(article_weights_raw)
  cat("Total number of individuals =", total_individuals, "\n")

  # Normalize article weights to sum to 1
  article_weights_norm <- article_weights_raw / sum(article_weights_raw)
  cat("Normalized article weights:", paste(round(article_weights_norm, 3), collapse = ", "), "\n")

  # Initialize combined data
  combined_data <- global_bins[, c("bin_lower", "bin_upper", "bin_midpoint", "bin_width")]
  combined_data$percentage <- 0
  combined_data$density <- 0

  for (i in seq_along(standardized_articles)) {
    combined_data$percentage <- combined_data$percentage +
      standardized_articles[[i]]$percentage * article_weights_norm[i]
  }

  # Re-normalize combined percentages
  total_pct <- sum(combined_data$percentage)
  if (abs(total_pct - 1) > 0.01) {
    cat("Combined percentages sum to", round(total_pct, 3), "- normalizing\n")
    combined_data$percentage <- combined_data$percentage / total_pct
  }
  combined_data$density <- combined_data$percentage / combined_data$bin_width

  cat("\nFinal standardized binned data (weighted average):\n")
  print(combined_data[, c("bin_lower", "bin_upper", "percentage", "density")])

  # ---------------------------------------------------------------------------
  # STEP 3: Generate final simulated data from standardized weighted bins
  # ---------------------------------------------------------------------------
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

  # ---------------------------------------------------------------------------
  # STEP 4: Kernel density estimation and threshold calculation
  # ---------------------------------------------------------------------------
  density_est <- stats::density(final_simulated, kernel = "gaussian", bw = "nrd0")
  density_df <- data.frame(x = density_est$x, y = density_est$y)

  # Optimal value (peak of distribution)
  optimal_value <- density_df$x[which.max(density_df$y)]

  # Preference range
  if (pref_method == "FWHM") {
    # Full Width at Half Maximum
    peak_height <- max(density_df$y) / 2
    f <- stats::approxfun(density_df$x, density_df$y - peak_height)
    lower_pref <- tryCatch(
      stats::uniroot(f, interval = c(min(density_df$x), optimal_value))$root,
      error = function(e) {
        warning("Could not find lower preference bound: ", e$message)
        return(NA)
      }
    )
    upper_pref <- tryCatch(
      stats::uniroot(f, interval = c(optimal_value, max(density_df$x)))$root,
      error = function(e) {
        warning("Could not find upper preference bound: ", e$message)
        return(NA)
      }
    )
  } else { # quantile method
    pref_quant <- stats::quantile(final_simulated, probs = pref_quantiles, na.rm = TRUE)
    lower_pref <- pref_quant[1]
    upper_pref <- pref_quant[2]
  }

  # Tolerance range (quantiles)
  tolerance_range <- stats::quantile(final_simulated, probs = toler_quantiles, na.rm = TRUE)

  # Summary data frame
  summary_df <- data.frame(
    parameter = c("optimal_value",
                  "lower_preferred", "upper_preferred",
                  "lower_tolerance", "upper_tolerance"),
    value = c(optimal_value, lower_pref, upper_pref,
              tolerance_range[1], tolerance_range[2]),
    unit = rep(unit, 5)
  )

  # ---------------------------------------------------------------------------
  # STEP 5: Return results with S3 class
  # ---------------------------------------------------------------------------
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
    pref_method = pref_method,
    pref_quantiles = if (pref_method == "quantile") pref_quantiles else NULL
  )
  class(result) <- "niche_range"

  # ---------------------------------------------------------------------------
  # S3 print method
  # ---------------------------------------------------------------------------
  print.niche_range <- function(x, ...) {
    cat("\n=== Species", x$variable_name, "Preference and Tolerance Range Estimates ===\n")
    cat("Method: Standardized bins + weighted average across articles\n")
    cat("Preference range method:", x$pref_method, "\n")
    if (x$pref_method == "quantile") {
      cat("Preference quantiles:", paste(x$pref_quantiles, collapse = ", "), "\n")
    }
    cat("Final sample size for KDE:", x$sample_size, "\n")
    cat("Total individuals across articles:", x$total_individuals, "\n")
    cat("Number of articles combined:", length(x$article_weights), "\n\n")

    cat("Optimal", x$variable_name, ":\n")
    cat("  Value:", round(x$summary$value[x$summary$parameter == "optimal_value"], 2), x$unit, "\n\n")

    cat("Preferred Range:\n")
    cat("  Lower:", round(x$summary$value[x$summary$parameter == "lower_preferred"], 2), x$unit, "\n")
    cat("  Upper:", round(x$summary$value[x$summary$parameter == "upper_preferred"], 2), x$unit, "\n\n")

    cat("Tolerance Range [", paste(toler_quantiles, collapse = "-"), " percentile]:\n", sep = "")
    cat("  Lower:", round(x$summary$value[x$summary$parameter == "lower_tolerance"], 2), x$unit, "\n")
    cat("  Upper:", round(x$summary$value[x$summary$parameter == "upper_tolerance"], 2), x$unit, "\n")
  }

  # ---------------------------------------------------------------------------
  # S3 plot method
  # ---------------------------------------------------------------------------
  plot.niche_range <- function(x, ...) {
    # Extract parameter values
    optimal <- x$summary$value[x$summary$parameter == "optimal_value"]
    lower_pref <- x$summary$value[x$summary$parameter == "lower_preferred"]
    upper_pref <- x$summary$value[x$summary$parameter == "upper_preferred"]
    lower_tol <- x$summary$value[x$summary$parameter == "lower_tolerance"]
    upper_tol <- x$summary$value[x$summary$parameter == "upper_tolerance"]

    # Create axis label
    x_label <- paste(x$variable_name, "(", x$unit, ")", sep = "")

    # Build plot
    p <- ggplot2::ggplot() +
      ggplot2::geom_col(data = x$standardized_data,
                        ggplot2::aes(x = bin_midpoint, y = density, width = bin_width),
                        fill = "lightblue", color = "black", alpha = 0.5) +
      ggplot2::geom_line(data = x$density_df, ggplot2::aes(x = x, y = y),
                         linewidth = 1.2, color = "darkblue") +
      ggplot2::geom_vline(xintercept = optimal, linetype = "dashed",
                          color = "red", linewidth = 1) +
      ggplot2::annotate("text", x = optimal, y = max(x$density_df$y) * 0.9,
                        label = paste("Optimal:", round(optimal, 1), x$unit),
                        hjust = -0.1, color = "red") +
      ggplot2::geom_vline(xintercept = c(lower_pref, upper_pref),
                          linetype = "dashed", color = "orange", linewidth = 1) +
      ggplot2::annotate("rect", xmin = lower_pref, xmax = upper_pref,
                        ymin = 0, ymax = Inf, fill = "orange", alpha = 0.1) +
      ggplot2::geom_vline(xintercept = c(lower_tol, upper_tol),
                          linetype = "dashed", color = "darkgreen", linewidth = 1) +
      ggplot2::annotate("rect", xmin = lower_tol, xmax = upper_tol,
                        ymin = 0, ymax = Inf, fill = "green", alpha = 0.05) +
      ggplot2::labs(title = paste("Species", x$variable_name, "Preference and Tolerance Range Estimates"),
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
    utils::write.csv(summary_df, paste0(save_path, "/summary_niche_fit.csv"), row.names = FALSE)
  }

  return(result)
}
