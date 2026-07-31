#' Estimate preference and tolerance ranges for a species based on occurrence data
#'
#' Extracts environmental values (e.g., temperature, dissolved oxygen)
#' from raster layers at occurrence points, then estimates optimal value,
#' preference range, and tolerance range using the same asymmetric mode-based
#' approach as \code{niche_tag()}.
#'
#' @param occ_path character. Path to occurrence file (table with columns: long, lat, year, month).
#' @param env_path character. Path to NetCDF file containing environmental data (e.g., temperature, DO).
#' @param pref_prob numeric. Proportion of area on each side of the mode to include in the
#' preference range. Default is 0.5 (50% of each side's area).
#' @param toler_prob numeric. Proportion of area on each side of the mode to include in the
#' tolerance range. Default is 0.95 (95% of each side's area).
#' @param bw_method Character. Bandwidth selection method for kernel density estimation.
#' Options are:
#' \itemize{
#'   \item \code{"SJ"} (default): Sheather–Jones direct plug-in estimator. Adaptive to
#'         data structure, preferred for complex/multimodal distributions, but
#'         computationally more intensive.
#'   \item \code{"nrd0"}: Silverman's rule-of-thumb (default in \code{stats::density}).
#'         Fast and stable, suitable for unimodal, near-normal distributions, but may
#'         oversmooth multimodal or skewed distributions.
#'   \item \code{"nrd"}: Silverman's rule-of-thumb (variant). Slightly different
#'         scaling factor, similar properties to \code{"nrd0"}.
#'   \item \code{"bcv"}: Biased cross-validation. Tends to produce smaller bandwidths,
#'         useful for fine-scale features, but can be unstable.
#'   \item \code{"ucv"}: Unbiased cross-validation. Tends to produce larger bandwidths,
#'         more stable than \code{"bcv"} but computationally intensive.
#' }
#' @param save_path character. Directory to save output plot and summary CSV. If NULL, nothing saved.
#' @param min_year numeric. First year of occurrence data (used to compute time index in NetCDF).
#' @param var_name character. Variable name in the NetCDF file (e.g., "thetao_mean", "o2").
#' @param env_label character. Display name for the environmental variable (e.g., "Temperature", "Dissolved Oxygen"). Defaults to var_name.
#' @param unit character. Unit of the environmental variable (e.g., "°C", "mg/L"). Defaults to empty string.
#' @param lon_range numeric vector of length 2: min and max longitude of the raster extent.
#' @param lat_range numeric vector of length 2: min and max latitude of the raster extent.
#' @param width numeric. Plot width in inches (default: 6).
#' @param height numeric. Plot height in inches (default: 5).
#' @param seed integer. Random seed for reproducibility (default: 123).
#'
#' @return A list containing:
#' \item{occ_values}{Extracted environmental values at occurrence points (NAs removed)}
#' \item{summary}{Data frame with optimal value, preference range, tolerance range}
#' \item{density_df}{Kernel density estimation results}
#' \item{pref_prob_used}{Proportion used for preference range on each side}
#' \item{toler_prob_used}{Proportion used for tolerance range on each side}
#' \item{bw_method}{Bandwidth selection method used}
#' \item{bw_value}{Bandwidth value used}
#' \item{env_label}{Display name of the variable}
#' \item{unit}{Unit of the variable}
#'
#' @examples
#' \donttest{
#' occ_path <- system.file("extdata", "R_typus_q1.txt", package = "OceanSDM")
#' env_path <- system.file("extdata", "cmems/IndoWPac_temp_1m_y2015_y2024.nc", package = "OceanSDM")
#'
#' niche_occ(
#'   occ_path = occ_path,
#'   env_path = env_path,
#'   pref_prob = 0.5,
#'   toler_prob = 0.95,
#'   bw_method = "SJ",
#'   save_path = "F:/whaleshark_sdm/temp_fit/habitat",
#'   min_year = 2015,
#'   var_name = "thetao_mean",
#'   env_label = "Temperature",
#'   unit = "°C",
#'   lon_range = c(90, 180),
#'   lat_range = c(-41, 40)
#' )
#' }
#' @export
niche_occ <- function(
    occ_path,
    env_path,
    pref_prob = 0.5,
    toler_prob = 0.95,
    bw_method = c("SJ", "nrd0", "nrd", "bcv", "ucv"),
    save_path = NULL,
    min_year,
    var_name = "thetao_mean",
    env_label = var_name,
    unit = "",
    lon_range,
    lat_range,
    width = 6,
    height = 5,
    seed = 123
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

  # ---- Read occurrence data ----
  df <- utils::read.table(occ_path, header = TRUE)

  # ---- Open NetCDF file ----
  nc <- ncdf4::nc_open(env_path)
  df$env_value <- NA

  cat("Extracting environmental values from NetCDF...\n")

  for (i in 1:nrow(df)) {
    # Compute time index (monthly timesteps from min_year)
    j <- (df$year[i] - min_year) * 12 + df$month[i]

    # Extract data slice depending on dimensionality
    if (length(dim(ncdf4::ncvar_get(nc, varid = var_name))) == 4) {
      env_slice <- ncdf4::ncvar_get(nc, varid = var_name)[,, 1, j]  # 4D: lon, lat, depth, time
    } else if (length(dim(ncdf4::ncvar_get(nc, varid = var_name))) == 3) {
      env_slice <- ncdf4::ncvar_get(nc, varid = var_name)[,, j]  # 3D: lon, lat, time
    } else {
      stop("Unsupported number of dimensions in NetCDF variable")
    }

    # Create SpatRaster from slice using terra
    env_raster <- terra::rast(
      t(env_slice),  # transpose to align with raster convention
      extent = terra::ext(lon_range[1], lon_range[2], lat_range[1], lat_range[2]),
      crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
    )

    # Flip vertically
    env_raster <- terra::flip(env_raster, direction = "vertical")

    # Extract value at occurrence point (columns: long, lat)
    extracted <- terra::extract(env_raster, df[i, c(3, 2)])
    df$env_value[i] <- extracted[1, 2]
  }

  ncdf4::nc_close(nc)

  # ---- Remove NA values ----
  occ_values <- stats::na.omit(df$env_value)

  cat("Summary of extracted", var_name, "values:\n")
  print(summary(occ_values))

  # ---- Kernel density estimation ----
  density_est <- stats::density(occ_values, kernel = "gaussian", bw = bw_method)
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

  # ---- Summary data frame ----
  summary_df <- data.frame(
    parameter = c("optimal_value",
                  "lower_preferred", "upper_preferred",
                  "lower_tolerance", "upper_tolerance"),
    value = c(optimal_value, lower_pref, upper_pref, lower_tol, upper_tol),
    unit = rep(unit, 5)
  )

  # ---- Build result list ----
  result <- list(
    occ_values = occ_values,
    summary = summary_df,
    density_df = density_df,
    pref_prob_used = pref_prob,
    toler_prob_used = toler_prob,
    bw_method = bw_method,
    bw_value = density_est$bw,
    env_label = env_label,
    unit = unit,
    pref_interval = pref_interval,
    toler_interval = toler_interval
  )

  class(result) <- "niche_occ"

  # ---- S3 plot method ----
  plot.niche_occ <- function(x, ...) {
    optimal <- x$summary$value[x$summary$parameter == "optimal_value"]
    lower_pref <- x$summary$value[x$summary$parameter == "lower_preferred"]
    upper_pref <- x$summary$value[x$summary$parameter == "upper_preferred"]
    lower_tol <- x$summary$value[x$summary$parameter == "lower_tolerance"]
    upper_tol <- x$summary$value[x$summary$parameter == "upper_tolerance"]

    occ_df <- data.frame(value = x$occ_values)
    x_label <- ifelse(x$unit == "", x$env_label, paste(x$env_label, "(", x$unit, ")"))

    pref_label <- paste("Preference (", x$pref_prob_used * 100, "% of each side)", sep = "")
    toler_label <- paste("Tolerance (", x$toler_prob_used * 100, "% of each side)", sep = "")

    p <- ggplot2::ggplot() +
      ggplot2::geom_histogram(data = occ_df,
                              ggplot2::aes(x = value, y = ggplot2::after_stat(density)),
                              binwidth = diff(range(occ_df$value)) / 20,
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
      ggplot2::annotate("text", x = lower_pref,
                        y = max(x$density_df$y) * 0.15,
                        label = toler_label,
                        color = "darkgreen", size = 3.5, hjust = 0) +
      ggplot2::annotate("text", x = lower_pref,
                        y = max(x$density_df$y) * 0.10,
                        label = pref_label,
                        color = "orange", size = 3.5, hjust = 0) +
      ggplot2::labs(title = paste("Species Preference and Tolerance Range Estimates for", x$env_label),
                    subtitle = paste("Based on", length(x$occ_values), "occurrence records"),
                    x = x_label,
                    y = "Density") +
      ggplot2::theme_minimal()

    return(p)
  }

  # ---- S3 print method ----
  print.niche_occ <- function(x, ...) {
    cat("\n=== Species", x$env_label, "Preference and Tolerance Range Estimates ===\n")
    cat("Data source: Occurrence records (", length(x$occ_values), " points)\n", sep = "")
    cat("Bandwidth method:", x$bw_method, "| Bandwidth value:", round(x$bw_value, 4), "\n")
    cat("Range estimation: Asymmetric interval using mode as dividing line\n")
    cat("  Each side independently covers target proportion of its own area\n")
    cat("  Preference: covers", x$pref_prob_used * 100, "% of each side's area\n")
    cat("  Tolerance:  covers", x$toler_prob_used * 100, "% of each side's area\n\n")

    cat("Optimal", x$env_label, ":\n")
    cat("  Value:", round(x$summary$value[x$summary$parameter == "optimal_value"], 2), x$unit, "\n\n")

    cat("Preferred Range:\n")
    cat("  Lower:", round(x$summary$value[x$summary$parameter == "lower_preferred"], 2), x$unit, "\n")
    cat("  Upper:", round(x$summary$value[x$summary$parameter == "upper_preferred"], 2), x$unit, "\n\n")

    cat("Tolerance Range:\n")
    cat("  Lower:", round(x$summary$value[x$summary$parameter == "lower_tolerance"], 2), x$unit, "\n")
    cat("  Upper:", round(x$summary$value[x$summary$parameter == "upper_tolerance"], 2), x$unit, "\n")
  }

  # ---- Print and plot ----
  print.niche_occ(result)
  print(plot.niche_occ(result))

  # ---- Save outputs if path provided ----
  if (!is.null(save_path)) {
    cat("\nSaving niche plot and data...\n")
    grDevices::tiff(paste0(save_path, "/Niche_fit.tif"),
                    width = width, height = height, units = "in", res = 300, compression = "lzw")
    print(plot.niche_occ(result))
    grDevices::dev.off()

    utils::write.csv(summary_df[, c("parameter", "value")],
                     paste0(save_path, "/summary_niche_fit.csv"),
                     row.names = FALSE)
  }

  return(result)
}
