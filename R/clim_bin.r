#' Calculate mean water temperature for each grid cell and depth layer within user-specified time periods
#'
#' @param clim_nc Character. Path to the multi-depth layers over years and months of climate NetCDF file (e.g., ocean temperature).
#' @param presence_tif List or RasterLayer. A named list where each element is a
#'   presence raster (values 1 = presence, NA/0 = absence) corresponding to a
#'   temporal group (folder name). Alternatively, a single raster (not recommended
#'   for multi‑group processing).
#' @param time_type Character. Type of temporal aggregation. One of `"month"`
#'   (monthly mean), `"quarter"` (quarterly mean), or `"season"` (seasonal mean).
#' @param save_path Character. Directory path where output files (organized in
#'   subfolders per temporal group) will be saved.
#' @param hemisphere Character. Hemisphere setting for seasonal aggregation.
#'   Either `"north"` (Northern Hemisphere) or `"south"` (Southern Hemisphere).
#'   Default is `"north"`.
#' @param na_value Numeric. Value used to replace `NA` in the output rasters.
#'   Default is `-9999`.
#'
#' @return Invisibly returns a character vector of all generated file paths.
#'   The function primarily writes ASCII raster files (`.asc`) into subfolders
#'   under `save_path`.
#'
#' @examples
#' \donttest{
#' library(raster)
#' quarter_files<-c(
#'   Q1="F:/whaleshark_sdm/result_sdm/Q1/Ensemble/SUP/MAX_TSS/Rhincodon_typus.tif",
#'   Q2="F:/whaleshark_sdm/result_sdm/Q2/Ensemble/SUP/MAX_TSS/Rhincodon_typus.tif",
#'   Q3="F:/whaleshark_sdm/result_sdm/Q3/Ensemble/SUP/MAX_TSS/Rhincodon_typus.tif",
#'   Q4="F:/whaleshark_sdm/result_sdm/Q4/Ensemble/SUP/MAX_TSS/Rhincodon_typus.tif")
#' presence_tif <- lapply(quarter_files, raster)
#' save_path<-"F:/whaleshark_sdm/4D"
#' clim_nc<-"F:/whaleshark_sdm/cmems/IndoWPac_temp_200m_y2015_y2024.nc"
#' clim_bin(clim_nc, presence_tif, time_type = "quarter",
#'          var_name="thetao_mean",
#'          save_path, hemisphere = "north", na_value = -9999)
#' }
#' @export
clim_bin <- function(clim_nc, presence_tif, time_type = c("month", "quarter", "season"),
                     var_name="thetao_mean",
                     save_path, hemisphere = "north", na_value = -9999) {

  time_type <- match.arg(time_type)
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

  # ---------- 1. Process species presence raster ----------
  # If presence_tif is a single raster, convert to a list (shared across all groups)

  # ---------- 2. Read NetCDF information ----------
  cat("Reading temperature NetCDF file...\n")
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

  # Temporal grouping
  if (time_type == "month") {
    time_groups <- paste0(years, "_", sprintf("%02d", months))
    group_names <- paste0("M", sprintf("%02d", months))
    folder_names <- paste0("M", sprintf("%02d", months))
  } else if (time_type == "quarter") {
    quarters <- lubridate::quarter(dates)
    time_groups <- paste0(years, "_Q", quarters)
    group_names <- paste0("Q", quarters)
    folder_names <- paste0("Q", quarters)
  } else if (time_type == "season") {
    if (hemisphere == "north") {
      seasons <- dplyr::case_when(
        months %in% c(3,4,5) ~ "MAM",
        months %in% c(6,7,8) ~ "JJA",
        months %in% c(9,10,11) ~ "SON",
        TRUE ~ "DJF"
      )
    } else {
      seasons <- dplyr::case_when(
        months %in% c(9,10,11) ~ "SON",
        months %in% c(12,1,2) ~ "DJF",
        months %in% c(3,4,5) ~ "MAM",
        TRUE ~ "JJA"
      )
    }
    time_groups <- paste0(years, "_", seasons)
    group_names <- seasons
    folder_names <- seasons
  }
  unique_groups <- unique(time_groups)
  cat("Total temporal groups:", length(unique_groups), "\n")

  # Temperature variable name
  if(var_name=="thetao_mean"){
    if (!(var_name %in% names(nc$var))) {
      var_names <- names(nc$var)
      var_name <- var_names[grep("thetao|temp|temperature", var_names, ignore.case = TRUE)][1]
      if (is.na(var_name)) stop("Unable to find temperature variable")
      cat("Using variable:", var_name, "\n")
    }
  }
  if(var_name=="o2"){
  if (!(var_name %in% names(nc$var))) {
    var_names <- names(nc$var)
    var_name <- var_names[grep("dissolved oxygen|oxygen", var_names, ignore.case = TRUE)][1]
    if (is.na(var_name)) stop("Unable to find dissolved oxygen variable")
    cat("Using variable:", var_name, "\n")
  }
  }

  # ---------- 3. Create folders for each temporal group ----------
  for (g in seq_along(unique_groups)) {
    folder_name <- folder_names[g]
    group_folder <- file.path(save_path, folder_name)
    if (!dir.exists(group_folder)) dir.create(group_folder, recursive = TRUE)
  }

  # ---------- 4. Compute means for each depth layer and temporal group ----------
  result_files <- list()
  # depth layer
  for (d in seq_along(depth)) {
    cat(sprintf("\nProcessing depth layer %.1f m...\n", depth[d]))

    # Read all time points for this depth layer (full grid)
    temp_data_full <- ncdf4::ncvar_get(nc, var_name, start = c(1,1,d,1), count = c(-1,-1,1,-1))
    # Dimensions: [lon, lat, time]
    dimnames(temp_data_full) <- list(lon = lon, lat = lat, time = 1:length(time))

    for (g in seq_along(unique_groups)) {
      group <- unique_groups[g]
      folder_name <- folder_names[g]
      cat(sprintf("  Processing group %s...\n", folder_name))

      # Get presence raster for this group
      if (is.list(presence_tif) && !is.null(names(presence_tif))) {
        # Named list, match by group name
        if (folder_name %in% names(presence_tif)) {
          pres_raster <- presence_tif[[folder_name]]
        } else {
          warning(sprintf("No presence raster for group %s, skipping", folder_name))
          next
        }
      } else {
        stop("Unable to find presence rasters")
      }

      # Extract valid presence pixels (value == 1)
      presence_cells <- which(raster::values(pres_raster) == 1)
      if (length(presence_cells) == 0) {
        warning(sprintf("Group %s has no presence pixels, skipping", folder_name))
        next
      }
      presence_xy <- raster::xyFromCell(pres_raster, presence_cells)

      # Find nearest NC grid indices for each presence pixel
      find_nc_idx <- function(xy) {
        lon_idx <- which.min(abs(lon - xy[1]))
        lat_idx <- which.min(abs(lat - xy[2]))
        return(c(lon_idx, lat_idx))
      }
      presence_indices <- t(apply(presence_xy, 1, find_nc_idx))
      colnames(presence_indices) <- c("lon_idx", "lat_idx")
      unique_indices <- unique(presence_indices)

      # Time indices belonging to this group
      time_idx <- which(time_groups == group)
      if (length(time_idx) == 0) next

      # Compute group mean for this depth layer (only for unique grid points)
      group_means <- numeric(nrow(unique_indices))
      for (i in seq_len(nrow(unique_indices))) {
        lon_idx <- unique_indices[i, "lon_idx"]
        lat_idx <- unique_indices[i, "lat_idx"]
        temp_series <- temp_data_full[lon_idx, lat_idx, time_idx]
        group_means[i] <- mean(temp_series, na.rm = TRUE)
      }

      # Create output raster (template from presence raster)
      result_raster <- raster::raster(pres_raster)
      raster::values(result_raster) <- NA

      # Assign means to corresponding presence pixels
      for (j in seq_along(presence_cells)) {
        cell <- presence_cells[j]
        idx_match <- which(unique_indices[, "lon_idx"] == presence_indices[j, "lon_idx"] &
                             unique_indices[, "lat_idx"] == presence_indices[j, "lat_idx"])
        if (length(idx_match) == 1) {
          result_raster[cell] <- group_means[idx_match]
        }
      }

      # Replace NA with specified value
      result_raster[is.na(result_raster)] <- na_value

      # Save
      group_folder <- file.path(save_path, folder_name)
      file_name <- sprintf("%s_temp_%.1f.asc", folder_name, depth[d])
      file_path <- file.path(group_folder, file_name)
      raster::writeRaster(result_raster, file_path, format = "ascii", overwrite = TRUE,
                          NAflag = na_value)
      result_files <- c(result_files, file_path)
    }
  }

  ncdf4::nc_close(nc)
  cat("\nProcessing completed! Total files generated:", length(result_files), "\n")
  invisible(result_files)
}
