#' Correct occurrence coordinates based on a standard raster layer of the predictor dataset
#'
#' @param occ_path Character. Path to occurrence CSV file. Must contain
#'   columns `decimalLongitude` and `decimalLatitude` (WGS84 coordinates).
#' @param predictor_raster_path Character. File path to the standard raster layer
#'   (e.g., an ASCII or GeoTIFF file) used as reference for correction.
#' @param max_distance_km Numeric. Maximum distance (in kilometres) that a point
#'   will be adjusted to the nearest grid cell with non-`NA` value. Default is `10`.
#' @param out_dir Character. **Full file path** (including `.csv` extension) where
#'   the corrected occurrence data will be saved. If the directory does not exist,
#'   it will be created. If a directory path is given instead, a default filename
#'   `corrected_occurrences.csv` will be appended.
#' @param verbose Logical. If `TRUE`, print progress messages. Default is `TRUE`.
#' @param plot Logical. If `TRUE`, generate a map showing correction status. Default `FALSE`.
#' @param plot_file Character. Optional file path to save the plot as TIFF (e.g., "plot.tif").
#'   If `NULL`, the plot is displayed in the graphics device.
#' @param xlim Numeric vector of length 2 (min, max) for longitude limits. If `NULL`,
#'   automatically determined from the raster extent.
#' @param ylim Numeric vector of length 2 (min, max) for latitude limits. If `NULL`,
#'   automatically determined from the raster extent.
#'
#' @return A data frame (invisibly) with the corrected occurrences.
#'
#' @importFrom terra rast extract res xmin xmax ymin ymax ext as.data.frame vect
#' @importFrom dplyr mutate filter if_else bind_rows
#' @importFrom sf st_as_sf
#' @importFrom ggplot2 ggplot geom_polygon geom_raster geom_sf scale_fill_gradient scale_color_manual coord_sf theme_bw theme labs element_blank element_text
#' @importFrom grDevices tiff dev.off
#' @export
correct_ob <- function(occ_path, predictor_raster_path,
                       max_distance_km = 10, out_dir, verbose = TRUE,
                       plot = FALSE, plot_file = NULL,
                       xlim = NULL, ylim = NULL) {

  # ---- 1. Validate input files ----
  if (missing(occ_path) || !file.exists(occ_path)) {
    stop("'occ_path' must be a valid file path to a CSV file.")
  }
  if (!file.exists(predictor_raster_path)) {
    stop("predictor_raster_path file does not exist: ", predictor_raster_path)
  }

  # ---- 2. Read occurrence data ----
  occurrence_data <- utils::read.csv(occ_path, stringsAsFactors = FALSE)
  required_cols <- c("decimalLongitude", "decimalLatitude")
  if (!all(required_cols %in% colnames(occurrence_data))) {
    missing <- setdiff(required_cols, colnames(occurrence_data))
    stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  }

  # ---- 3. Validate and prepare output file path ----
  if (missing(out_dir) || is.null(out_dir)) {
    stop("'out_dir' must be provided.")
  }
  if (dir.exists(out_dir) || grepl("[\\\\/]$", out_dir)) {
    out_dir <- file.path(out_dir, "corrected_occurrences.csv")
  }
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  parent_dir <- dirname(out_dir)
  if (!dir.exists(parent_dir)) {
    dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # ---- 4. Read raster ----
  if (verbose) cat("Reading predictor raster...\n")
  temp_raster <- terra::rast(predictor_raster_path)
  original_n <- nrow(occurrence_data)

  # ---- 5. Create copy for plotting (preserve original coordinates) ----
  plot_data <- occurrence_data
  plot_data$status <- "original"
  plot_data$corr_x <- NA_real_
  plot_data$corr_y <- NA_real_

  # ---- 6. Add tracking columns to main data ----
  occurrence_data <- occurrence_data |>
    dplyr::mutate(
      original_x = decimalLongitude,
      original_y = decimalLatitude,
      correction_status = "original",
      correction_distance_km = 0,
      correction_x = NA_real_,
      correction_y = NA_real_
    )

  # ---- 7. Extract initial values ----
  coords <- occurrence_data[, c("decimalLongitude", "decimalLatitude")]
  extracted <- terra::extract(temp_raster, coords)
  if (ncol(extracted) >= 2) {
    occurrence_data$env <- extracted[, -1]
  } else {
    occurrence_data$env <- NA_real_
  }

  # ---- 8. Find points with NA ----
  invalid_indices <- which(is.na(occurrence_data$env))
  if (verbose) {
    cat(sprintf("Found %d points with NA predictor values\n", length(invalid_indices)))
  }

  # ---- 9. Use two-SpatVector nearest to find closest non-NA cells ----
  if (length(invalid_indices) > 0) {
    # Subset points with NA
    pts_na_df <- occurrence_data[invalid_indices, c("decimalLongitude", "decimalLatitude")]

    # Get raster CRS (OGC:CRS84)
    raster_crs <- terra::crs(temp_raster)

    # Create points with the SAME CRS as the raster
    pts_na <- terra::vect(pts_na_df,
                          geom = c("decimalLongitude", "decimalLatitude"),
                          crs = raster_crs)

    # Extract all non-NA cells as points (inherits raster CRS)
    non_na_points <- terra::as.points(temp_raster, values = TRUE, na.rm = TRUE)

    # Now both objects have IDENTICAL CRS
    nearest_info <- terra::nearest(pts_na, non_na_points)
    nearest_df <- as.data.frame(nearest_info)

    # Column names: "from_id", "to_id", "distance"
    dist_col <- grep("distance", names(nearest_df), ignore.case = TRUE, value = TRUE)[1]
    to_col <- grep("to_id", names(nearest_df), ignore.case = TRUE, value = TRUE)[1]

    if (is.na(dist_col) || is.na(to_col)) {
      # Fallback: use first two numeric columns
      num_cols <- names(nearest_df)[sapply(nearest_df, is.numeric)]
      if (length(num_cols) >= 2) {
        to_col <- num_cols[1]
        dist_col <- num_cols[2]
      } else {
        stop("Cannot locate 'to_id' and 'distance' columns in nearest output.")
      }
    }

    distances_m <- nearest_df[[dist_col]]
    to_indices <- nearest_df[[to_col]]

    # Loop through each NA point and update
    for (k in seq_along(invalid_indices)) {
      i <- invalid_indices[k]
      to_idx <- to_indices[k]
      if (!is.na(to_idx)) {
        # Get coordinates of the nearest non-NA point
        coords_cell <- terra::geom(non_na_points)[to_idx, c("x", "y")]
        new_x <- coords_cell[1]
        new_y <- coords_cell[2]
        # Extract value from non_na_points
        raster_name <- names(temp_raster)[1]
        new_val <- as.data.frame(non_na_points)[to_idx, raster_name]
        # Distance in km
        dist_km <- distances_m[k] / 1000

        # ---- Check distance limit ----
        if (dist_km > max_distance_km) {
          # Point is too far -> remove it
          if (verbose) {
            cat(sprintf("Point %d removed: nearest cell at %.2f km (exceeds max %.1f km)\n",
                        i, dist_km, max_distance_km))
          }
          occurrence_data$correction_status[i] <- "removed"
          plot_data$status[i] <- "removed"
        } else {
          # Valid correction
          occurrence_data$original_x[i] <- occurrence_data$decimalLongitude[i]
          occurrence_data$original_y[i] <- occurrence_data$decimalLatitude[i]
          occurrence_data$decimalLongitude[i] <- new_x
          occurrence_data$decimalLatitude[i] <- new_y
          occurrence_data$env[i] <- new_val
          occurrence_data$correction_status[i] <- "corrected"
          occurrence_data$correction_distance_km[i] <- dist_km
          occurrence_data$correction_x[i] <- new_x
          occurrence_data$correction_y[i] <- new_y

          plot_data$status[i] <- "corrected"
          plot_data$corr_x[i] <- new_x
          plot_data$corr_y[i] <- new_y

          if (verbose) {
            cat(sprintf("Point %d corrected at distance %.2f km\n", i, dist_km))
          }
        }
      } else {
        # No non-NA cell found (should not happen)
        occurrence_data$correction_status[i] <- "removed"
        plot_data$status[i] <- "removed"
        if (verbose) {
          cat(sprintf("Point %d could not be corrected (no non-NA cells found)\n", i))
        }
      }
    }
  }

  # ---- 10. Plot if requested ----
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' is required for plotting. Please install it.")
    } else {
      # Prepare raster data frame
      raster_df <- as.data.frame(temp_raster, xy = TRUE, na.rm = FALSE)
      colnames(raster_df) <- c("x", "y", "value")
      raster_df <- raster_df[!is.na(raster_df$value), ]

      # Build point categories
      orig_pts <- plot_data[plot_data$status == "original", c("decimalLongitude", "decimalLatitude")]
      if (nrow(orig_pts) > 0) orig_pts$Category <- "valid_original"

      corr_idx <- which(plot_data$status == "corrected")
      if (length(corr_idx) > 0) {
        corr_orig <- data.frame(
          decimalLongitude = plot_data$decimalLongitude[corr_idx],
          decimalLatitude = plot_data$decimalLatitude[corr_idx],
          Category = "corrected_original"
        )
        corr_new <- data.frame(
          decimalLongitude = plot_data$corr_x[corr_idx],
          decimalLatitude = plot_data$corr_y[corr_idx],
          Category = "corrected_new"
        )
      } else {
        corr_orig <- data.frame()
        corr_new <- data.frame()
      }

      rem_idx <- which(plot_data$status == "removed")
      if (length(rem_idx) > 0) {
        rem_pts <- data.frame(
          decimalLongitude = plot_data$decimalLongitude[rem_idx],
          decimalLatitude = plot_data$decimalLatitude[rem_idx],
          Category = "removed"
        )
      } else {
        rem_pts <- data.frame()
      }

      plot_points <- dplyr::bind_rows(orig_pts, corr_orig, corr_new, rem_pts)
      plot_points_sf <- sf::st_as_sf(plot_points,
                                     coords = c("decimalLongitude", "decimalLatitude"),
                                     crs = 4326)

      # Determine plot extent
      if (is.null(xlim) || is.null(ylim)) {
        ext <- terra::ext(temp_raster)
        if (is.null(xlim)) xlim <- c(ext[1], ext[2])
        if (is.null(ylim)) ylim <- c(ext[3], ext[4])
      }

      # ---- Build plot ----
      p <- ggplot2::ggplot()

      # Land background (20% grey)
      if (requireNamespace("maps", quietly = TRUE)) {
        world <- ggplot2::map_data("world")
        p <- p + ggplot2::geom_polygon(data = world,
                                       ggplot2::aes(x = long, y = lat, group = group),
                                       fill = "grey80", color = NA)
      } else {
        message("Package 'maps' not available, land background skipped.")
      }

      # Raster layer (green to red gradient)
      raster_name <- names(temp_raster)[1]
      if (is.null(raster_name) || is.na(raster_name)) raster_name <- "Value"
      p <- p + ggplot2::geom_raster(data = raster_df,
                                    ggplot2::aes(x = x, y = y, fill = value),
                                    alpha = 0.8) +
        ggplot2::scale_fill_gradient(
          low = "green", high = "red",
          name = raster_name,
          na.value = "transparent",
          guide = ggplot2::guide_colorbar(direction = "vertical")
        )

      # Points with categories
      p <- p + ggplot2::geom_sf(data = plot_points_sf,
                                ggplot2::aes(color = Category),
                                size = 1.5, alpha = 0.7)
      p <- p + ggplot2::scale_color_manual(
        values = c(
          "valid_original" = "black",
          "corrected_original" = "red",
          "corrected_new" = "blue",
          "removed" = "yellow"
        ),
        labels = c(
          "valid_original" = "Valid (original)",
          "corrected_original" = "Corrected (original)",
          "corrected_new" = "Corrected (new)",
          "removed" = "Removed"
        )
      )

      p <- p + ggplot2::coord_sf(xlim = xlim, ylim = ylim, crs = 4326) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position = "right",
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::labs(
          title = "Occurrence Correction",
          subtitle = paste("Max search distance:", max_distance_km, "km"),
          x = NULL,
          y = NULL
        )

      # Save or display
      if (!is.null(plot_file)) {
        grDevices::tiff(plot_file, width = 8, height = 6, units = "in", res = 300, compression = "lzw")
        print(p)
        grDevices::dev.off()
        if (verbose) cat("Plot saved to:", plot_file, "\n")
      } else {
        print(p)
      }
    }
  }

  # ---- 11. Filter out removed points and prepare final data ----
  occurrence_data <- occurrence_data |>
    dplyr::filter(correction_status != "removed")

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

  # ---- 12. Write to CSV ----
  utils::write.csv(corrected_all_data, out_dir, row.names = FALSE)

  if (verbose) {
    cat(sprintf("Original records: %d\n", original_n))
    cat(sprintf("Records removed: %d\n", original_n - nrow(corrected_all_data)))
    cat(sprintf("Records corrected: %d\n", sum(occurrence_data$correction_status == "corrected")))
    cat(sprintf("Final records: %d\n", nrow(corrected_all_data)))
    cat("Output saved to:", out_dir, "\n")
  }

  invisible(corrected_all_data)
}
