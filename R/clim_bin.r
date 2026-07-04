#' Calculate mean water temperature (or any variable) for each grid cell and depth layer
#' within user-specified time periods. For seasonal aggregation, hemisphere is
#' automatically determined by latitude (north: >=0, south: <0).
#'
#' @param clim_nc Character. Path to the multi-depth layers over years
#'   and months of climate NetCDF file.
#' @param presence_asc List or SpatRaster. A named list where each element is a
#'   presence raster (values 1 = presence, NA/0 = absence) corresponding to a
#'   temporal group (folder name). For season, names must be one of
#'   "Winter", "Spring", "Summer", "Autumn".
#' @param time_type Character. Type of temporal aggregation. One of `"month"`,
#'   `"quarter"`, or `"season"`.
#' @param var_name Character. Variable name in NetCDF. If not found exactly,
#'   the function will attempt to match (case‑insensitive) and warn.
#'   e.g., `"thetao_mean"`.
#' @param save_path Character. Directory path where output files will be saved.
#' @param na_value Numeric. Value to replace `NA`. Default `-9999`.
#'
#' @return Invisibly returns a character vector of all generated file paths.
#' @export
clim_bin <- function(clim_nc, presence_asc, time_type = c("month", "quarter", "season"),
                     var_name = "thetao_mean",
                     save_path, na_value = -9999) {

  time_type <- match.arg(time_type)
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

  # ---- 1. Convert presence rasters to list of SpatRaster ----
  to_spatraster <- function(x) {
    if (inherits(x, "SpatRaster")) return(x)
    if (inherits(x, "RasterLayer")) return(terra::rast(x))
    if (is.character(x) && file.exists(x)) return(terra::rast(x))
    stop("Input must be a SpatRaster, RasterLayer, or file path to a raster.")
  }

  if (inherits(presence_asc, "RasterLayer") ||
      inherits(presence_asc, "SpatRaster") ||
      (is.character(presence_asc) && length(presence_asc) == 1 && file.exists(presence_asc))) {
    presence_asc <- to_spatraster(presence_asc)
    presence_asc <- list(presence_asc)
  } else if (is.list(presence_asc)) {
    presence_asc <- lapply(presence_asc, function(x) {
      if (inherits(x, "RasterLayer") || is.character(x)) to_spatraster(x) else x
    })
  } else {
    stop("presence_asc must be a SpatRaster, RasterLayer, file path, or a list of such objects.")
  }

  if (is.list(presence_asc) && is.null(names(presence_asc))) {
    stop("presence_asc must be a named list with folder names as names (e.g., Q1, Q2, ...)")
  }

  # ---- 2. Read NetCDF metadata ----
  cat("Reading NetCDF file...\n")
  nc <- ncdf4::nc_open(clim_nc)
  lon <- ncdf4::ncvar_get(nc, "longitude")
  lat <- ncdf4::ncvar_get(nc, "latitude")
  depth <- ncdf4::ncvar_get(nc, "depth")
  time <- ncdf4::ncvar_get(nc, "time")

  # Time conversion
  time_units <- ncdf4::ncatt_get(nc, "time", "units")$value
  time_origin <- stringr::str_extract(time_units, "\\d{4}-\\d{2}-\\d{2}")
  if (grepl("days since", time_units)) {
    dates <- as.Date(time, origin = time_origin)
  } else if (grepl("seconds since", time_units)) {
    dates <- as.Date(as.POSIXct(time, origin = time_origin, tz = "UTC"))
  } else {
    dates <- seq(as.Date("2015-01-01"), by = "month", length.out = length(time))
  }
  years <- lubridate::year(dates)
  months <- lubridate::month(dates)

  # ---- 3. Temporal grouping (only for month/quarter, season handled separately) ----
  if (time_type == "month") {
    folder_names <- paste0("M", sprintf("%02d", 1:12))
    # Actually, we only need the unique folder names that exist in presence_asc
    # We will filter later.
  } else if (time_type == "quarter") {
    folder_names <- paste0("Q", 1:4)
  } else { # season
    folder_names <- c("Winter", "Spring", "Summer", "Autumn")
  }

  # Only keep folder_names that have corresponding presence rasters
  valid_folders <- intersect(folder_names, names(presence_asc))
  if (length(valid_folders) == 0) {
    stop("None of the presence raster names match the expected folder names for time_type.")
  }
  folder_names <- valid_folders

  cat("Processing groups:", paste(folder_names, collapse = ", "), "\n")

  # ---- 4. Variable name resolution ----
  available_vars <- names(nc$var)
  if (var_name %in% available_vars) {
    cat("Using variable:", var_name, "\n")
  } else {
    matches <- grep(tolower(var_name), tolower(available_vars), value = TRUE, fixed = FALSE)
    if (length(matches) == 1) {
      var_name <- matches[1]
      cat("Variable name not found exactly. Using closest match:", var_name, "\n")
    } else if (length(matches) > 1) {
      ncdf4::nc_close(nc)
      stop("Multiple variables match '", var_name, "': ",
           paste(matches, collapse = ", "),
           ". Please specify one explicitly.")
    } else {
      ncdf4::nc_close(nc)
      stop("Variable '", var_name, "' not found in NetCDF. Available variables: ",
           paste(available_vars, collapse = ", "))
    }
  }

  # ---- 5. Define season month sets for both hemispheres ----
  north_season_months <- list(
    Winter = c(12, 1, 2),
    Spring = c(3, 4, 5),
    Summer = c(6, 7, 8),
    Autumn = c(9, 10, 11)
  )
  south_season_months <- list(
    Winter = c(6, 7, 8),
    Spring = c(9, 10, 11),
    Summer = c(12, 1, 2),
    Autumn = c(3, 4, 5)
  )

  # ---- 6. Create output folders ----
  for (folder_name in folder_names) {
    group_folder <- file.path(save_path, folder_name)
    if (!dir.exists(group_folder)) dir.create(group_folder, recursive = TRUE)
  }

  # ---- 7. Precompute presence cell indices and time indices for each group ----
  group_info <- list()
  for (folder_name in folder_names) {
    pres_raster <- presence_asc[[folder_name]]
    if (!inherits(pres_raster, "SpatRaster")) {
      pres_raster <- terra::rast(pres_raster)
    }
    pres_vals <- terra::values(pres_raster, mat = FALSE)
    presence_cells <- which(pres_vals == 1)
    if (length(presence_cells) == 0) {
      warning(sprintf("Group %s has no presence pixels, skipping", folder_name))
      next
    }
    presence_xy <- terra::xyFromCell(pres_raster, presence_cells)

    # Find nearest NC indices
    lon_idx <- apply(presence_xy[, 1, drop = FALSE], 1, function(x) which.min(abs(lon - x)))
    lat_idx <- apply(presence_xy[, 2, drop = FALSE], 1, function(x) which.min(abs(lat - x)))
    nc_lin_idx <- (lat_idx - 1) * length(lon) + lon_idx

    # Compute time indices for this group based on time_type
    if (time_type == "month") {
      # folder_name is like "M01", "M02", ...
      month_num <- as.numeric(substr(folder_name, 2, 3))
      time_idx <- which(months == month_num)
    } else if (time_type == "quarter") {
      # folder_name is like "Q1", "Q2", ...
      quarter_num <- as.numeric(substr(folder_name, 2, 2))
      time_idx <- which(lubridate::quarter(dates) == quarter_num)
    } else { # season
      time_idx <- NULL  # will be handled per cell based on latitude
    }

    group_info[[folder_name]] <- list(
      presence_cells = presence_cells,
      nc_lin_idx = nc_lin_idx,
      lat_vals = presence_xy[, 2],
      time_idx = time_idx
    )
  }

  if (length(group_info) == 0) stop("No valid presence groups found.")

  # ---- 8. Compute means for each depth layer and group ----
  result_files <- list()

  for (d in seq_along(depth)) {
    cat(sprintf("\nProcessing depth layer %.1f m...\n", depth[d]))

    # Read entire time series for this depth layer
    temp_data_full <- ncdf4::ncvar_get(nc, var_name, start = c(1,1,d,1), count = c(-1,-1,1,-1))
    # Reshape to (lon*lat) x time
    nlon <- length(lon)
    nlat <- length(lat)
    ntime <- length(time)
    temp_mat <- matrix(temp_data_full, nrow = nlon * nlat, ncol = ntime)

    for (folder_name in names(group_info)) {
      cat(sprintf("  Processing group %s...\n", folder_name))
      info <- group_info[[folder_name]]
      presence_cells <- info$presence_cells
      nc_lin_idx <- info$nc_lin_idx
      lat_vals <- info$lat_vals

      # Determine time indices and compute means
      if (time_type %in% c("month", "quarter")) {
        time_idx <- info$time_idx
        if (length(time_idx) == 0) {
          warning(sprintf("No time indices for group %s", folder_name))
          next
        }
        sub_mat <- temp_mat[nc_lin_idx, time_idx, drop = FALSE]
        group_means <- rowMeans(sub_mat, na.rm = TRUE)
      } else { # season
        # Split cells into north and south
        north_mask <- lat_vals >= 0
        south_mask <- !north_mask
        group_means <- numeric(length(presence_cells))

        if (any(north_mask)) {
          month_set <- north_season_months[[folder_name]]
          time_idx <- which(months %in% month_set)
          if (length(time_idx) > 0) {
            sub_mat <- temp_mat[nc_lin_idx[north_mask], time_idx, drop = FALSE]
            group_means[north_mask] <- rowMeans(sub_mat, na.rm = TRUE)
          } else {
            group_means[north_mask] <- NA
          }
        }
        if (any(south_mask)) {
          month_set <- south_season_months[[folder_name]]
          time_idx <- which(months %in% month_set)
          if (length(time_idx) > 0) {
            sub_mat <- temp_mat[nc_lin_idx[south_mask], time_idx, drop = FALSE]
            group_means[south_mask] <- rowMeans(sub_mat, na.rm = TRUE)
          } else {
            group_means[south_mask] <- NA
          }
        }
      }

      # Create output raster
      template_raster <- presence_asc[[folder_name]]
      if (!inherits(template_raster, "SpatRaster")) {
        template_raster <- terra::rast(template_raster)
      }
      result_raster <- terra::rast(template_raster)
      terra::values(result_raster) <- NA
      result_raster[presence_cells] <- group_means
      result_raster[is.na(result_raster)] <- na_value

      # Save
      group_folder <- file.path(save_path, folder_name)
      file_name <- sprintf("%s_temp_%.1f.asc", folder_name, depth[d])
      file_path <- file.path(group_folder, file_name)
      terra::writeRaster(result_raster, file_path, filetype = "AAIGrid",
                         overwrite = TRUE, NAflag = na_value)
      result_files <- c(result_files, file_path)
    }
  }

  ncdf4::nc_close(nc)
  cat("\nProcessing completed! Total files generated:", length(result_files), "\n")
  invisible(result_files)
}
