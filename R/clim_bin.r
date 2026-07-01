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
#'   the function will attempt to match (case‑insensitive) and warn. Default is
#'   `"thetao_mean"`.
#' @param save_path Character. Directory path where output files will be saved.
#' @param na_value Numeric. Value to replace `NA`. Default `-9999`.
#' @return Invisibly returns a character vector of all generated file paths.
#' @export
clim_bin <- function(clim_nc, presence_asc, time_type = c("month", "quarter", "season"),
                     var_name = "thetao_mean",
                     save_path, na_value = -9999) {

  time_type <- match.arg(time_type)
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

  # ---- 1. Convert presence rasters to SpatRaster (terra) ----
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

  # ---- 2. Read NetCDF information ----
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

  # ---- 3. Temporal grouping ----
  if (time_type == "month") {
    time_groups <- paste0(years, "_", sprintf("%02d", months))
    folder_names <- paste0("M", sprintf("%02d", months))
    unique_groups <- unique(time_groups)
  } else if (time_type == "quarter") {
    quarters <- lubridate::quarter(dates)
    time_groups <- paste0(years, "_Q", quarters)
    folder_names <- paste0("Q", quarters)
    unique_groups <- unique(time_groups)
  } else if (time_type == "season") {
    # Seasonal aggregation: use English season names
    season_names <- c("Winter", "Spring", "Summer", "Autumn")
    folder_names <- season_names
    unique_groups <- season_names
    time_groups <- rep(NA, length(time))  # placeholder
  }

  cat("Total temporal groups:", length(unique_groups), "\n")

  # ---- 4. Universal variable name resolution ----
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

  # ---- 5. Define season month sets for both hemispheres (using season names) ----
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

  # ---- 6. Create folders for each temporal group ----
  for (folder_name in unique(folder_names)) {
    group_folder <- file.path(save_path, folder_name)
    if (!dir.exists(group_folder)) dir.create(group_folder, recursive = TRUE)
  }

  # ---- 7. Compute means for each depth layer and temporal group ----
  result_files <- list()

  for (d in seq_along(depth)) {
    cat(sprintf("\nProcessing depth layer %.1f m...\n", depth[d]))

    # Read all time points for this depth layer
    temp_data_full <- ncdf4::ncvar_get(nc, var_name, start = c(1,1,d,1), count = c(-1,-1,1,-1))
    dimnames(temp_data_full) <- list(lon = lon, lat = lat, time = 1:length(time))

    for (g in seq_along(unique_groups)) {
      group <- unique_groups[g]
      folder_name <- folder_names[g]
      cat(sprintf(" Processing group %s...\n", folder_name))

      # Get presence raster for this group
      if (is.list(presence_asc) && !is.null(names(presence_asc))) {
        if (folder_name %in% names(presence_asc)) {
          pres_raster <- presence_asc[[folder_name]]
        } else {
          warning(sprintf("No presence raster for group %s, skipping", folder_name))
          next
        }
      } else {
        stop("Unable to find presence rasters")
      }

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
      find_nc_idx <- function(xy) {
        c(which.min(abs(lon - xy[1])), which.min(abs(lat - xy[2])))
      }
      presence_indices <- t(apply(presence_xy, 1, find_nc_idx))
      colnames(presence_indices) <- c("lon_idx", "lat_idx")
      unique_indices <- unique(presence_indices)

      # ---- Compute means for this group ----
      if (time_type == "season") {
        # For each unique grid point, determine hemisphere and corresponding month set
        group_means <- numeric(nrow(unique_indices))
        season_name <- folder_name  # e.g., "Winter"
        for (i in seq_len(nrow(unique_indices))) {
          lon_idx <- unique_indices[i, "lon_idx"]
          lat_idx <- unique_indices[i, "lat_idx"]
          lat_val <- lat[lat_idx]
          # Choose month set based on latitude (north: >=0, south: <0)
          if (lat_val >= 0) {
            month_set <- north_season_months[[season_name]]
          } else {
            month_set <- south_season_months[[season_name]]
          }
          time_idx <- which(months %in% month_set)
          if (length(time_idx) == 0) {
            group_means[i] <- NA
            next
          }
          temp_series <- temp_data_full[lon_idx, lat_idx, time_idx]
          group_means[i] <- mean(temp_series, na.rm = TRUE)
        }
      } else {
        # For month and quarter: use global time indices
        if (time_type == "month") {
          time_idx <- which(time_groups == group)
        } else if (time_type == "quarter") {
          time_idx <- which(time_groups == group)
        }
        if (length(time_idx) == 0) next
        group_means <- numeric(nrow(unique_indices))
        for (i in seq_len(nrow(unique_indices))) {
          lon_idx <- unique_indices[i, "lon_idx"]
          lat_idx <- unique_indices[i, "lat_idx"]
          temp_series <- temp_data_full[lon_idx, lat_idx, time_idx]
          group_means[i] <- mean(temp_series, na.rm = TRUE)
        }
      }

      # Create output raster and assign means
      result_raster <- terra::rast(pres_raster)
      terra::values(result_raster) <- NA
      for (j in seq_along(presence_cells)) {
        cell <- presence_cells[j]
        idx_match <- which(unique_indices[, "lon_idx"] == presence_indices[j, "lon_idx"] &
                             unique_indices[, "lat_idx"] == presence_indices[j, "lat_idx"])
        if (length(idx_match) == 1) {
          result_raster[cell] <- group_means[idx_match]
        }
      }
      result_raster[is.na(result_raster)] <- na_value

      # Save as ASCII grid
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
