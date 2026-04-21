#' Check if occurrences of a single species fall within a defined range map with a buffer
#'
#' @param species_name Character. Scientific name of the species.
#' @param occ_data Data frame. Occurrence data with columns `decimalLongitude`,
#'   `decimalLatitude`, and other optional columns (e.g., Species, year, month).
#' @param range_path Character. Path to the range map shapefile (must contain a column `sci_name`).
#' @param buffer_km Numeric. Buffer distance (km) around the range map to consider
#'   points as valid (default: 100).
#' @param distshore_km Numeric. Distance to shore (km) used to exclude points
#'   that are on land (default: 10).
#'
#' @return A list with elements:
#'   \item{valid_points}{Data frame of points that passed all checks.}
#'   \item{invalid_points}{Data frame of points that failed the checks.}
#'   \item{error}{Character string with error message if any, otherwise `NULL`.}
#'
#' @details The function performs three steps:
#'   1. Check if points are within the buffered range map.
#'   2. For points outside the buffer, try four sign combinations of coordinates.
#'   3. Remove points that fall on land (using `obistools::lookup_xy`).
#'
#'
#' @examples
#' \donttest{
#' # Assume `occ` is a data frame with decimalLongitude, decimalLatitude columns
#' occ<-read.csv("F:/whaleshark_sdm/occ/R_typus_adj.csv")
#' result <- range_check_single("Rhincodon typus",
#'                              occ,
#'                              "F:/whaleshark_sdm/range_map/R_typus_wip2.shp",
#'                               buffer_km = 100)
#' }
#' @export
range_check_single <- function(species_name, occ_data, range_path, buffer_km = 100, distshore_km = 10) {
  # ---- Input validation ----
  if (!file.exists(range_path)) stop("range_path file does not exist: ", range_path)
  if (!is.data.frame(occ_data)) stop("occ_data must be a data frame")
  required_cols <- c("decimalLongitude", "decimalLatitude")
  if (!all(required_cols %in% colnames(occ_data))) {
    stop("occ_data must contain columns: decimalLongitude, decimalLatitude")
  }

  cat("Processing species:", species_name, "- number of points:", nrow(occ_data), "\n")

  # ---- Read range map ----
  range.map <- sf::st_read(range_path, quiet = TRUE)

  # ---- Error handling wrapper ----
  result <- tryCatch({
    # If species not in range map, keep all points
    if (!species_name %in% range.map$sci_name) {
      cat("  Species not found in range map, keeping all points\n")
      return(list(
        valid_points = occ_data,
        invalid_points = data.frame(),
        error = NULL
      ))
    }

    # Extract species range polygon
    species_range <- range.map |> dplyr::filter(sci_name == species_name)

    if (nrow(species_range) == 0) {
      cat("  Unable to extract species range, keeping all points\n")
      return(list(
        valid_points = occ_data,
        invalid_points = data.frame(),
        error = NULL
      ))
    }

    # Make geometry valid
    species_range <- sf::st_make_valid(species_range)

    # Create buffered range (step 1 buffer)
    buffer_dist_m <- buffer_km * 1000
    species_range_buffered <- sf::st_buffer(species_range, dist = buffer_dist_m)

    # Convert points to sf object
    points_sf <- sf::st_as_sf(occ_data,
                              coords = c("decimalLongitude", "decimalLatitude"),
                              crs = sf::st_crs(range.map))

    # Step 1: Check points inside buffered range
    cat("  Step 1: Checking points inside buffered range map\n")
    step1_intersection <- sf::st_intersects(points_sf, species_range_buffered)
    step1_in_range <- sapply(step1_intersection, function(x) length(x) > 0)

    step1_valid_points <- occ_data[step1_in_range, ]
    if (nrow(step1_valid_points) > 0) {
      step1_valid_points <- as.data.frame(step1_valid_points)
      step1_valid_points$validation_step <- "step1_buffer"
    } else {
      step1_valid_points <- data.frame()
    }

    step1_invalid_points <- occ_data[!step1_in_range, ]
    cat("  Step 1 valid points:", nrow(step1_valid_points),
        "invalid points:", nrow(step1_invalid_points), "\n")

    # Step 2: Try sign combinations for invalid points
    step2_corrected_points <- data.frame()
    step2_final_invalid_points <- data.frame()

    if (nrow(step1_invalid_points) > 0) {
      cat("  Step 2: Trying sign combinations for coordinates\n")

      sign_combinations <- list(
        c(1, 1),   # original
        c(-1, 1),  # flip longitude
        c(1, -1),  # flip latitude
        c(-1, -1)  # flip both
      )

      for (i in seq_along(sign_combinations)) {
        combo <- sign_combinations[[i]]
        combo_name <- paste0("sign_combo_", i)

        corrected_coords <- step1_invalid_points
        corrected_coords$decimalLongitude <- corrected_coords$decimalLongitude * combo[1]
        corrected_coords$decimalLatitude <- corrected_coords$decimalLatitude * combo[2]

        corrected_coords_sf <- sf::st_as_sf(corrected_coords,
                                            coords = c("decimalLongitude", "decimalLatitude"),
                                            crs = sf::st_crs(range.map))

        step2_intersection <- sf::st_intersects(corrected_coords_sf, species_range_buffered)
        step2_in_range <- sapply(step2_intersection, function(x) length(x) > 0)

        if (any(step2_in_range)) {
          combo_corrected <- step1_invalid_points[step2_in_range, ]
          combo_corrected$decimalLongitude <- combo_corrected$decimalLongitude * combo[1]
          combo_corrected$decimalLatitude <- combo_corrected$decimalLatitude * combo[2]
          combo_corrected <- as.data.frame(combo_corrected)
          combo_corrected$correction_type <- combo_name
          combo_corrected$validation_step <- "step2_sign_corrected"

          step2_corrected_points <- dplyr::bind_rows(step2_corrected_points, combo_corrected)
          cat("    Combination", i, "corrected points:", nrow(combo_corrected), "\n")

          # Remove corrected points from invalid set
          step1_invalid_points <- step1_invalid_points[!step2_in_range, ]
        }
      }
    }

    # Step 3: Remove points on land using obistools
    occ_1_2 <- dplyr::bind_rows(step1_valid_points, step2_corrected_points)

    if (nrow(occ_1_2) > 0) {
      # Use obistools to compute shore distance
      xydata <- obistools::lookup_xy(occ_1_2, shoredistance = TRUE, grids = TRUE, areas = TRUE)
      dist <- as.data.frame(xydata$shoredistance)
      colnames(dist) <- "Distance"
      occ_1_2 <- cbind(occ_1_2, dist)

      # Keep points with distance > -1000 * distshore_km (i.e., not on land)
      all_valid_points <- occ_1_2[occ_1_2$Distance > (-1000 * distshore_km), ]
      step3_invalid <- occ_1_2[occ_1_2$Distance <= (-1000 * distshore_km), ]
      step3_invalid <- dplyr::select(step3_invalid, -Distance)
    } else {
      all_valid_points <- data.frame()
      step3_invalid <- data.frame()
    }

    final_invalid_points <- rbind(step1_invalid_points, step3_invalid)
    cat("  Final valid points:", nrow(all_valid_points),
        "final invalid points:", nrow(final_invalid_points), "\n")

    list(
      valid_points = all_valid_points,
      invalid_points = final_invalid_points,
      error = NULL
    )
  }, error = function(e) {
    cat("  Error processing species", species_name, ":", e$message, "\n")
    cat("  Keeping all original points as valid\n")
    list(
      valid_points = occ_data,   # fixed: was 'points_subset'
      invalid_points = data.frame(),
      error = e$message
    )
  })

  return(result)
}
