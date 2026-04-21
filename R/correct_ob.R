#' Correct occurrence coordinates based on a standard raster layer of the predictor dataset
#'
#' @param occurrence_data Data frame. Occurrence data to be corrected. Must contain
#'   columns `decimalLongitude` and `decimalLatitude` (WGS84 coordinates).
#' @param predictor_raster_path Character. File path to the standard raster layer
#'   (e.g., an ASCII or GeoTIFF file) used as reference for correction.
#' @param max_distance_km Numeric. Maximum distance (in kilometres) that a point
#'   will be adjusted to the nearest grid cell with non‑`NA` value. Default is `10`.
#' @param out_dir Character. Directory path (including file name) where the corrected
#'   occurrence data will be saved as a CSV file.
#' @param verbose Logical. If `TRUE`, print progress messages. Default is `TRUE`.
#'
#' @return A data frame (invisibly) with the corrected occurrences. The returned
#'   data frame includes additional columns: `original_x`, `original_y`,
#'   `correction_status` (`"original"`, `"corrected"`, or `"removed"`),
#'   `correction_distance_km`, `correction_x`, `correction_y`, and `temperature`.
#'   The function also writes a CSV file to `out_dir`.
#'
#' @examples
#' \donttest{
#' corrected_ob <- correct_ob(
#'   occurrence_data = final_combined_data,
#'   predictor_raster_path = "C:/shark/bin/Q1_thetao_mean.asc",
#'   max_distance_km = 10,
#'   out_dir = "F:/shark/occ/R_typus_adj.csv",
#'   verbose = TRUE
#' )
#' }
#' @export
correct_ob <- function(occurrence_data, predictor_raster_path, max_distance_km = 10,
                       out_dir, verbose = TRUE) {

  # Validate input columns
  required_cols <- c("decimalLongitude", "decimalLatitude")
  if (!all(required_cols %in% colnames(occurrence_data))) {
    missing <- setdiff(required_cols, colnames(occurrence_data))
    stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  }

  # Read raster
  if (verbose) cat("Reading predictor raster...\n")
  temp_raster <- raster::raster(predictor_raster_path)

  original_n <- nrow(occurrence_data)

  # Add tracking columns using dplyr (or base R)
  occurrence_data <- occurrence_data |>
    dplyr::mutate(
      original_x = decimalLongitude,
      original_y = decimalLatitude,
      correction_status = "original",
      correction_distance_km = 0,
      correction_x = NA_real_,
      correction_y = NA_real_
    )

  # Convert points to SpatialPoints (using sp and raster)
  points_sp <- sp::SpatialPoints(
    occurrence_data[, c("decimalLongitude", "decimalLatitude")],
    proj4string = sp::CRS(raster::proj4string(temp_raster))
  )
  occurrence_data$temperature <- raster::extract(temp_raster, points_sp)

  degree_per_km <- 1 / 111.32
  max_distance_degrees <- max_distance_km * degree_per_km
  search_step <- raster::res(temp_raster)[1]

  invalid_indices <- which(is.na(occurrence_data$temperature))

  if (verbose) {
    cat(sprintf("Found %d points with NA temperature values\n", length(invalid_indices)))
  }

  for (i in invalid_indices) {
    original_x <- occurrence_data$decimalLongitude[i]
    original_y <- occurrence_data$decimalLatitude[i]
    found <- FALSE
    search_distance <- 0

    while (search_distance <= max_distance_degrees && !found) {
      search_distance <- search_distance + search_step
      n_points <- max(8, round(2 * pi * search_distance / search_step))
      angles <- seq(0, 2 * pi, length.out = n_points + 1)[1:n_points]
      search_points_x <- original_x + search_distance * cos(angles)
      search_points_y <- original_y + search_distance * sin(angles)

      for (j in seq_len(n_points)) {
        search_x <- search_points_x[j]
        search_y <- search_points_y[j]

        if (search_x >= raster::xmin(temp_raster) && search_x <= raster::xmax(temp_raster) &&
            search_y >= raster::ymin(temp_raster) && search_y <= raster::ymax(temp_raster)) {

          temp_value <- raster::extract(temp_raster, cbind(search_x, search_y))

          if (!is.na(temp_value)) {
            distance_km <- sqrt((search_x - original_x)^2 + (search_y - original_y)^2) * 111.32
            occurrence_data$correction_status[i] <- "corrected"
            occurrence_data$correction_distance_km[i] <- distance_km
            occurrence_data$correction_x[i] <- search_x
            occurrence_data$correction_y[i] <- search_y
            occurrence_data$temperature[i] <- temp_value
            occurrence_data$decimalLongitude[i] <- search_x
            occurrence_data$decimalLatitude[i] <- search_y
            found <- TRUE
            break
          }
        }
      }
    }
    if (!found) {
      occurrence_data$correction_status[i] <- "removed"
      if (verbose) {
        cat(sprintf("Point %d (lon=%f, lat=%f) could not be corrected – will be removed\n",
                    i, original_x, original_y))
      }
    }
  }

  # Filter out removed points
  occurrence_data <- occurrence_data |>
    dplyr::filter(correction_status != "removed")

  # Create final data frame with corrected coordinates (keep original columns)
  corrected_all_data <- occurrence_data |>
    dplyr::mutate(
      decimalLongitude = dplyr::if_else(
        correction_status == "corrected" & !is.na(correction_x),
        correction_x,
        original_x
      ),
      decimalLatitude = dplyr::if_else(
        correction_status == "corrected" & !is.na(correction_y),
        correction_y,
        original_y
      )
    )

  # Write to CSV
  utils::write.csv(corrected_all_data, out_dir, row.names = FALSE)

  if (verbose) {
    cat(sprintf("Original records: %d\n", original_n))
    cat(sprintf("Records removed: %d\n", sum(occurrence_data$correction_status == "removed")))
    cat(sprintf("Records corrected: %d\n", sum(occurrence_data$correction_status == "corrected")))
    cat(sprintf("Final records: %d\n", nrow(corrected_all_data)))
    cat("Output saved to:", out_dir, "\n")
  }

  invisible(corrected_all_data)
}
