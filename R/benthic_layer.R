#' Create the benthic layer of an environmental factor
#'
#' @param input_nc Character. File path to the input NetCDF file downloaded from
#'   CMEMS (Copernicus Marine Environment Monitoring Service). The file must
#'   contain multiple depth layers for the specified variable.
#' @param output_nc Character. File path where the output NetCDF file containing
#'   the benthic layer will be written.
#' @param var_name Character. Name of the environmental variable to process
#'   (e.g., `"thetao_mean"`), as defined in the CMEMS user manual.
#'
#' @return Invisibly returns the output file path (`output_nc`). The function
#'   primarily writes a NetCDF file to disk containing the benthic (deepest
#'   non‑NA) layer for each time step.
#' @examples
#' \donttest{
#' # Ensure you have downloaded a NetCDF file with multiple depth layers
#' input_nc <- "C:/cmems/Pacific_temp_05_50_2011.nc"
#' output_nc <- "C:/cmems/Pacific_bent_temp_05_50_2011.nc"
#' var_name <- "thetao_mean"
#' benthic_layer(input_nc, output_nc, var_name)
#' }
#' @export
benthic_layer <- function(input_nc, output_nc, var_name) {
  # Step 1: Read metadata and depth values
  nc <- ncdf4::nc_open(input_nc)
  lon <- ncdf4::ncvar_get(nc, "longitude")
  lat <- ncdf4::ncvar_get(nc, "latitude")
  time <- ncdf4::ncvar_get(nc, "time")
  depth_var <- ifelse("depth" %in% names(nc$dim), "depth",
                      ifelse("z" %in% names(nc$dim), "z", "lev"))  # Common depth names
  depth_vals <- ncdf4::ncvar_get(nc, depth_var)

  # Ensure depth is positive downward (oceanographic convention)
  # if (mean(depth_vals, na.rm = TRUE) < 0) depth_vals <- abs(depth_vals)

  # Step 2: Process benthic temperature by time slice
  benthic_stack <- terra::rast()  # Initialize empty SpatRaster for output

  for (t in seq_along(time)) {
    cat("Processing time", t, "/", length(time), "\n")

    # Read 3D temperature array [lon, lat, depth]
    temp_array <- ncdf4::ncvar_get(
      nc,
      var_name,  # Replace with your temperature variable name
      start = c(1, 1, 1, t),
      count = c(-1, -1, -1, 1)
    )

    # Reshape to 2D matrix (rows = grid cells, columns = depth levels)
    dims <- dim(temp_array)
    temp_mat <- matrix(temp_array, nrow = prod(dims[1:2]), ncol = dims[3])

    # Find deepest non-NA value per grid cell
    benthic_vec <- sapply(1:nrow(temp_mat), function(i) {
      vals <- temp_mat[i, ]
      non_na <- which(!is.na(vals))
      if (length(non_na) == 0) return(NA)
      deepest_idx <- non_na[which.max(depth_vals[non_na])]
      vals[deepest_idx]
    })

    # Convert to raster layer and add to stack
    benthic_layer <- terra::rast(
      matrix(benthic_vec, nrow = dims[1], ncol = dims[2]),
      extent = c(min(lon), max(lon), min(lat), max(lat)),
      crs = "+proj=longlat"
    )
    benthic_stack <- c(benthic_stack, benthic_layer)
  }
  ncdf4::nc_close(nc)

  # Step 3: Set time dimension and save
  terra::time(benthic_stack) <- time
  names(benthic_stack) <- paste0("time_", seq_along(time))

  inter <- strsplit(as.character(time[1]), split = "-")
  inter2 <- strsplit(as.character(time[length(time)]), split = "-")
  min_month <- as.numeric(sapply(inter, function(x) x[2]))
  min_year <- as.numeric(sapply(inter, function(x) x[1]))
  max_month <- as.numeric(sapply(inter2, function(x) x[2]))
  max_year <- as.numeric(sapply(inter2, function(x) x[1]))

  # Write to NetCDF
  terra::writeCDF(
    benthic_stack,
    filename = output_nc,
    varname = var_name,
    unit = "degC",
    longname = "Temperature at the deepest layer",
    overwrite = TRUE
  )
}
