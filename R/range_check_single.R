#' Check if occurrences of a single species fall within a defined range map with a buffer
#'
#' @param species_name Character. Scientific name of the species.
#' @param occ_data Either a data frame with columns `decimalLongitude`, `decimalLatitude`,
#'   or a character string path to a CSV file containing those columns.
#' @param range_path Character. Path to the range map shapefile (must contain a column `sci_name`).
#' @param buffer_km Numeric. Buffer distance (km) around the range map to consider
#'   points as valid (default: 100).
#' @param plot Logical. If `TRUE`, generate a map showing point categories. Default `FALSE`.
#' @param plot_file Character. Optional file path to save the plot as TIFF (e.g., "plot.tif").
#'   If `NULL`, the plot is displayed in the graphics device.
#' @param xlim Numeric vector of length 2 (min, max) for longitude limits. If `NULL`,
#'   automatically determined from the buffered range.
#' @param ylim Numeric vector of length 2 (min, max) for latitude limits. If `NULL`,
#'   automatically determined from the buffered range.
#'
#' @return A list with elements:
#'   \item{valid_points}{Data frame of points that passed the checks (within buffer or corrected).}
#'   \item{invalid_points}{Data frame of points that failed all checks.}
#'   \item{correction_log}{Data frame with original and corrected coordinates for points that were fixed.}
#'   \item{error}{Character string with error message if any, otherwise `NULL`.}
#'
#' @details The function performs two steps:
#' 1. Check if points are within the buffered range map.
#' 2. For points outside the buffer, try four sign combinations of coordinates.
#' If `plot = TRUE`, a map is generated showing:
#'   - Black: valid points (no correction)
#'   - Red: original coordinates of corrected points
#'   - blue: corrected coordinates
#'   - Yellow: invalid points (still outside buffer after corrections)
#'   - Grey80 background: land polygons
#' The original species range polygon is shown as a red dashed border.
#'
#' @importFrom ggplot2 ggplot geom_sf geom_point scale_color_manual coord_sf theme_bw theme element_rect
#' @importFrom sf st_read st_as_sf st_crs st_buffer st_intersects st_make_valid st_bbox
#' @importFrom dplyr filter bind_rows
#' @importFrom rnaturalearth ne_countries
#' @export
range_check_single <- function(species_name, occ_data, range_path,
                               buffer_km = 100,
                               plot = FALSE,
                               plot_file = NULL,
                               xlim = NULL,
                               ylim = NULL) {

  # ---- Input validation ----
  if (!file.exists(range_path)) stop("range_path file does not exist: ", range_path)

  # ---- Determine if occ_data is a data frame or a file path ----
  if (is.data.frame(occ_data)) {
    occ_df <- occ_data
  } else if (is.character(occ_data) && length(occ_data) == 1 && file.exists(occ_data)) {
    occ_df <- utils::read.csv(occ_data, stringsAsFactors = FALSE)
  } else {
    stop("occ_data must be either a data frame or a valid file path to a CSV file")
  }

  required_cols <- c("decimalLongitude", "decimalLatitude")
  if (!all(required_cols %in% colnames(occ_df))) {
    stop("occ_data must contain columns: decimalLongitude, decimalLatitude")
  }

  cat("Processing species:", species_name, "- number of points:", nrow(occ_df), "\n")

  # ---- Read range map ----
  range.map <- sf::st_read(range_path, quiet = TRUE)

  # ---- Initialize correction log ----
  correction_log <- data.frame(
    original_long = numeric(),
    original_lat = numeric(),
    corrected_long = numeric(),
    corrected_lat = numeric(),
    combo = integer()
  )

  # ---- Error handling wrapper ----
  result <- tryCatch({
    if (!species_name %in% range.map$sci_name) {
      cat(" Species not found in range map, keeping all points\n")
      return(list(
        valid_points = occ_df,
        invalid_points = data.frame(),
        correction_log = correction_log,
        error = NULL
      ))
    }

    species_range <- range.map |> dplyr::filter(sci_name == species_name)
    if (nrow(species_range) == 0) {
      cat(" Unable to extract species range, keeping all points\n")
      return(list(
        valid_points = occ_df,
        invalid_points = data.frame(),
        correction_log = correction_log,
        error = NULL
      ))
    }

    species_range <- sf::st_make_valid(species_range)
    buffer_dist_m <- buffer_km * 1000
    species_range_buffered <- sf::st_buffer(species_range, dist = buffer_dist_m)

    points_sf <- sf::st_as_sf(occ_df,
                              coords = c("decimalLongitude", "decimalLatitude"),
                              crs = sf::st_crs(range.map))

    # Step 1: Check points inside buffered range
    cat(" Step 1: Checking points inside buffered range map\n")
    step1_intersection <- sf::st_intersects(points_sf, species_range_buffered)
    step1_in_range <- sapply(step1_intersection, function(x) length(x) > 0)
    step1_valid_points <- occ_df[step1_in_range, , drop = FALSE]
    step1_invalid_points <- occ_df[!step1_in_range, , drop = FALSE]
    cat(" Step 1 valid points:", nrow(step1_valid_points),
        "invalid points:", nrow(step1_invalid_points), "\n")

    # Step 2: Try sign combinations for invalid points
    step2_corrected_points <- data.frame()
    if (nrow(step1_invalid_points) > 0) {
      cat(" Step 2: Trying sign combinations for coordinates\n")
      sign_combinations <- list(c(1, 1), c(-1, 1), c(1, -1), c(-1, -1))
      for (i in seq_along(sign_combinations)) {
        combo <- sign_combinations[[i]]
        corrected_coords <- step1_invalid_points
        corrected_coords$decimalLongitude <- corrected_coords$decimalLongitude * combo[1]
        corrected_coords$decimalLatitude <- corrected_coords$decimalLatitude * combo[2]
        corrected_coords_sf <- sf::st_as_sf(corrected_coords,
                                            coords = c("decimalLongitude", "decimalLatitude"),
                                            crs = sf::st_crs(range.map))
        step2_intersection <- sf::st_intersects(corrected_coords_sf, species_range_buffered)
        step2_in_range <- sapply(step2_intersection, function(x) length(x) > 0)
        if (any(step2_in_range)) {
          combo_corrected <- step1_invalid_points[step2_in_range, , drop = FALSE]
          # Record original and corrected coordinates
          orig_long <- combo_corrected$decimalLongitude
          orig_lat  <- combo_corrected$decimalLatitude
          corr_long <- orig_long * combo[1]
          corr_lat  <- orig_lat * combo[2]
          # Update coordinates to corrected values
          combo_corrected$decimalLongitude <- corr_long
          combo_corrected$decimalLatitude <- corr_lat

          # Append to correction log
          log_entries <- data.frame(
            original_long = orig_long,
            original_lat = orig_lat,
            corrected_long = corr_long,
            corrected_lat = corr_lat,
            combo = i
          )
          correction_log <- dplyr::bind_rows(correction_log, log_entries)

          step2_corrected_points <- dplyr::bind_rows(step2_corrected_points, combo_corrected)
          cat(" Combination", i, "corrected points:", nrow(combo_corrected), "\n")
          step1_invalid_points <- step1_invalid_points[!step2_in_range, , drop = FALSE]
        }
      }
    }

    # Combine valid points: step1_valid_points + step2_corrected_points
    all_valid_points <- dplyr::bind_rows(step1_valid_points, step2_corrected_points)
    final_invalid_points <- step1_invalid_points

    cat(" Final valid points:", nrow(all_valid_points),
        "final invalid points:", nrow(final_invalid_points), "\n")

    # ---- Plot if requested ----
    if (plot) {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("Package 'ggplot2' is required for plotting. Please install it.")
      } else {
        # ---- Load land data ----
        land <- NULL
        if (requireNamespace("rnaturalearth", quietly = TRUE)) {
          tryCatch({
            land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
          }, error = function(e) {
            warning("Could not download land data. Proceeding without land background.")
          })
        } else {
          warning("Package 'rnaturalearth' not installed. Land background not available.")
        }

        # Build plot data frame
        # Valid original points (black)
        valid_orig <- step1_valid_points
        valid_orig$category <- "valid_original"

        # Corrected points (original coordinates in red, corrected in blue)
        if (nrow(correction_log) > 0) {
          corrected_orig <- data.frame(
            decimalLongitude = correction_log$original_long,
            decimalLatitude = correction_log$original_lat,
            category = "corrected_original"
          )
          corrected_new <- data.frame(
            decimalLongitude = correction_log$corrected_long,
            decimalLatitude = correction_log$corrected_lat,
            category = "corrected_new"
          )
        } else {
          corrected_orig <- data.frame()
          corrected_new <- data.frame()
        }

        # Invalid points (yellow)
        if (nrow(final_invalid_points) > 0) {
          invalid_pts <- final_invalid_points
          invalid_pts$category <- "invalid"
        } else {
          invalid_pts <- data.frame()
        }

        # Combine all points for plotting
        plot_points <- dplyr::bind_rows(
          valid_orig,
          corrected_orig,
          corrected_new,
          invalid_pts
        )

        # Convert to sf
        plot_points_sf <- sf::st_as_sf(plot_points,
                                       coords = c("decimalLongitude", "decimalLatitude"),
                                       crs = sf::st_crs(range.map))

        # Determine plot extent if not provided by user
        if (is.null(xlim) || is.null(ylim)) {
          buffered_bbox <- sf::st_bbox(species_range_buffered)
          if (is.null(xlim)) xlim <- c(buffered_bbox["xmin"], buffered_bbox["xmax"])
          if (is.null(ylim)) ylim <- c(buffered_bbox["ymin"], buffered_bbox["ymax"])
        }

        # Create plot
        p <- ggplot2::ggplot()

        # Add land background (grey80) if available
        if (!is.null(land)) {
          p <- p + ggplot2::geom_sf(data = land, fill = "grey80", color = NA)
        }

        # Add range polygon (red dashed border)
        p <- p + ggplot2::geom_sf(data = species_range,
                                  fill = NA, color = "red", linetype = "dashed", linewidth = 1)

        # Add points with category colors
        p <- p + ggplot2::geom_sf(data = plot_points_sf,
                                  ggplot2::aes(color = category),
                                  size = 1.5, alpha = 0.7) +
          ggplot2::scale_color_manual(
            values = c(
              "valid_original" = "black",
              "corrected_original" = "red",
              "corrected_new" = "blue",
              "invalid" = "yellow"
            ),
            labels = c(
              "valid_original" = "Valid (original)",
              "corrected_original" = "Corrected (original)",
              "corrected_new" = "Corrected (new)",
              "invalid" = "Invalid"
            )
          ) +
          ggplot2::coord_sf(xlim = xlim, ylim = ylim, crs = sf::st_crs(range.map)) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.position = "bottom"
          ) +
          ggplot2::labs(
            title = paste("Range Check for", species_name),
            subtitle = paste("Buffer:", buffer_km, "km")
          )

        # Save or display
        if (!is.null(plot_file)) {
          grDevices::tiff(plot_file, width = 8, height = 6, units = "in", res = 300, compression = "lzw")
          print(p)
          grDevices::dev.off()
          cat("Plot saved to:", plot_file, "\n")
        } else {
          print(p)
        }
      }
    }

    list(
      valid_points = all_valid_points,
      invalid_points = final_invalid_points,
      correction_log = correction_log,
      error = NULL
    )

  }, error = function(e) {
    cat(" Error processing species", species_name, ":", e$message, "\n")
    cat(" Keeping all original points as valid\n")
    list(
      valid_points = occ_df,
      invalid_points = data.frame(),
      correction_log = data.frame(),
      error = e$message
    )
  })

  return(result)
}
