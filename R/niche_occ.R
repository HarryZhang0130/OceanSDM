#' Estimate preference and tolerance ranges for a species based on occurrence data
#'
#' Extracts environmental values (e.g., temperature, dissolved oxygen) from raster layers
#' at occurrence points, then estimates optimal value, preference range, and tolerance range.
#'
#' @param occ_path character. Path to occurrence file (table with columns: long, lat, year, month).
#' @param env_path character. Path to NetCDF file containing environmental data (e.g., temperature, DO).
#' @param toler_quantiles vector. Quantiles for tolerance range (default: c(0.01, 0.99)).
#' @param pref_method character. Method for preference range: "FWHM" (Full Width at Half Maximum) or "quantile" (default: "FWHM").
#' @param pref_quantiles vector. Quantiles for preference range when pref_method = "quantile" (default: c(0.25, 0.75)).
#' @param save_path character. Directory to save output plot and summary CSV. If NULL, nothing saved.
#' @param min_year numeric. First year of occurrence data (used to compute time index in NetCDF).
#' @param var_name character. Variable name in the NetCDF file (e.g., "thetao_mean", "o2").
#' @param env_label character. Display name for the environmental variable (e.g., "Temperature", "Dissolved Oxygen"). Defaults to var_name.
#' @param unit character. Unit of the environmental variable (e.g., "°C", "mg/L"). Defaults to empty string.
#' @param lon_range numeric vector of length 2: min and max longitude of the raster extent.
#' @param lat_range numeric vector of length 2: min and max latitude of the raster extent.
#' @param width numeric. Plot width in inches (default: 6).
#' @param height numeric. Plot height in inches (default: 5).
#'
#' @return A list containing:
#'   \item{occ_values}{Extracted environmental values at occurrence points (NAs removed)}
#'   \item{summary}{Data frame with optimal value, preference range, tolerance range}
#'   \item{density_df}{Kernel density estimation results}
#'   \item{pref_method}{Method used for preference range}
#'   \item{pref_quantiles}{Preference quantiles (if applicable)}
#'   \item{env_label}{Display name of the variable}
#'   \item{unit}{Unit of the variable}
#'
#' @examples
#' \donttest{
#' niche_occ(
#'   occ_path = "F:/whaleshark_sdm/occ/R_typus_q1.txt",
#'   env_path = "F:/whaleshark_sdm/cmems/IndoWPac_temp_1m_y2015_y2024.nc",
#'   toler_quantiles = c(0.01, 0.99),
#'   pref_method = "FWHM",
#'   pref_quantiles =NULL,
#'   save_path = "F:/whaleshark_sdm/temp_fit/habitat"
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
    toler_quantiles = c(0.01, 0.99),
    pref_method = c("FWHM", "quantile"),
    pref_quantiles = c(0.25, 0.75),
    save_path = NULL,
    min_year,
    var_name = "thetao_mean",
    env_label = var_name,
    unit = "",
    lon_range,
    lat_range,
    width = 6,
    height = 5
) {
  # Match preference method argument
  pref_method <- match.arg(pref_method)

  # Read occurrence data
  df <- utils::read.table(occ_path, header = TRUE)

  nc <- ncdf4::nc_open(env_path)
  df$env_value <- NA

  for (i in 1:nrow(df)) {
    # Compute time index (monthly timesteps from min_year)
    j <- (df$year[i] - min_year) * 12 + df$month[i]
    # Extract data slice depending on dimensionality
    if (length(dim(ncdf4::ncvar_get(nc, varid = var_name))) == 4) {
      env_slice <- ncdf4::ncvar_get(nc, varid = var_name)[,, 1, j]  # 4D: lon, lat, depth, time
    } else if (length(dim(ncdf4::ncvar_get(nc, varid = var_name))) == 3) {
      env_slice <- ncdf4::ncvar_get(nc, varid = var_name)[,, j]     # 3D: lon, lat, time
    } else {
      stop("Unsupported number of dimensions in NetCDF variable")
    }

    # Create raster from slice
    env_raster <- raster::raster(t(env_slice),
                                 xmn = lon_range[1], xmx = lon_range[2],
                                 ymn = lat_range[1], ymx = lat_range[2])
    env_raster <- raster::flip(env_raster, direction = 2)
    raster::projection(env_raster) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'

    # Extract value at occurrence point (column 3 = long, column 2 = lat)
    df$env_value[i] <- raster::extract(env_raster, df[i, c(3, 2)])
  }

  ncdf4::nc_close(nc)

  # Remove NA values
  occ_values <- stats::na.omit(df$env_value)
  cat("Summary of extracted", var_name, "values:\n")
  print(summary(occ_values))

  # Kernel density estimation
  density_est <- stats::density(occ_values, kernel = "gaussian", bw = "nrd0")
  density_df <- data.frame(x = density_est$x, y = density_est$y)

  # Optimal value (peak of KDE)
  optimal_value <- density_df$x[which.max(density_df$y)]

  # Preference range
  if (pref_method == "FWHM") {
    # Full Width at Half Maximum
    peak_height <- max(density_df$y) / 2
    f <- stats::approxfun(density_df$x, density_df$y - peak_height)
    lower_pref <- tryCatch(
      stats::uniroot(f, interval = c(min(density_df$x), optimal_value))$root,
      error = function(e) NA
    )
    upper_pref <- tryCatch(
      stats::uniroot(f, interval = c(optimal_value, max(density_df$x)))$root,
      error = function(e) NA
    )
  } else { # quantile method
    pref_quant <- stats::quantile(occ_values, probs = pref_quantiles, na.rm = TRUE)
    lower_pref <- pref_quant[1]
    upper_pref <- pref_quant[2]
  }

  # Tolerance range (quantiles)
  tolerance_range <- stats::quantile(occ_values, probs = toler_quantiles, na.rm = TRUE)
  lower_tol <- tolerance_range[1]
  upper_tol <- tolerance_range[2]

  # Print results
  cat(sprintf("Preferred Range (%s): %.2f - %.2f (Optimal: %.2f)\n",
              ifelse(pref_method == "FWHM", "FWHM",
                     paste0(round(pref_quantiles[1]*100), "%-", round(pref_quantiles[2]*100), "%")),
              lower_pref, upper_pref, optimal_value))
  cat(sprintf("Tolerance Range (%d%%-%d%%): %.2f - %.2f\n",
              round(toler_quantiles[1]*100), round(toler_quantiles[2]*100),
              lower_tol, upper_tol))

  # Summary data frame
  summary_df <- data.frame(
    parameter = c("optimal_value", "lower_preferred", "upper_preferred",
                  "lower_tolerance", "upper_tolerance"),
    value = c(optimal_value, lower_pref, upper_pref, lower_tol, upper_tol),
    unit = rep(unit, 5)
  )

  # Build result list
  result <- list(
    occ_values = occ_values,
    summary = summary_df,
    density_df = density_df,
    pref_method = pref_method,
    pref_quantiles = if (pref_method == "quantile") pref_quantiles else NULL,
    env_label = env_label,
    unit = unit
  )
  class(result) <- "niche_occ"

  # Define S3 plot method
  plot.niche_occ <- function(x, ...) {
    # Extract values
    optimal <- x$summary$value[x$summary$parameter == "optimal_value"]
    lower_pref <- x$summary$value[x$summary$parameter == "lower_preferred"]
    upper_pref <- x$summary$value[x$summary$parameter == "upper_preferred"]
    lower_tol <- x$summary$value[x$summary$parameter == "lower_tolerance"]
    upper_tol <- x$summary$value[x$summary$parameter == "upper_tolerance"]

    # Create data frame for histogram
    occ_df <- data.frame(value = x$occ_values)

    # Axis label
    x_label <- ifelse(x$unit == "", x$env_label, paste(x$env_label, "(", x$unit, ")"))

    p <- ggplot2::ggplot() +
      ggplot2::geom_histogram(data = occ_df,
                              ggplot2::aes(x = value, y = ggplot2::after_stat(density)),
                              binwidth = diff(range(occ_df$value)) / 20,
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
      ggplot2::labs(title = paste("Species Preference and Tolerance Range Estimates for", x$env_label),
                    subtitle = paste("Based on", length(x$occ_values), "occurrence records"),
                    x = x_label,
                    y = "Density") +
      ggplot2::theme_minimal()
    return(p)
  }

  # Print and plot
  print(plot.niche_occ(result))

  # Save outputs if path provided
  if (!is.null(save_path)) {
    cat("Saving niche plot and data...\n")
    grDevices::tiff(paste0(save_path, "/Niche_fit.tif"),
                    width = width, height = height, units = "in", res = 300, compression = "lzw")
    print(plot.niche_occ(result))
    grDevices::dev.off()
    utils::write.csv(summary_df, paste0(save_path, "/summary_df.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  }

  return(result)
}
