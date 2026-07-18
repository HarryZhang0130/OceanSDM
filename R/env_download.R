#' Download environmental layers from CMEMS
#'
#' @param path_cmems_tool Character. Path to the copernicusmarine executable,
#'   e.g. `"/miniconda/envs/copernicusmarine/Scripts/copernicusmarine.exe"`.
#' @param out_dir Character. Directory where output data will be stored.
#' @param productId Character. ID of the CMEMS product (see dataset user manual
#'   in the CMEMS product webpage).
#'   e.g., `"cmems_mod_glo_phy-mnstd_my_0.25deg_P1M-m"`.
#' @param variable Character. Variable name to download. e.g., "thetao_mean".
#' @param date_min Date or character. Minimum date. e.g., "1993-01-01".
#' @param date_max Date or character. Maximum date. e.g., "2022-12-01".
#' @param lon Numeric vector of length 2. Longitude range (min, max).
#'   e.g., `c(-180, 179.75)`.
#' @param lat Numeric vector of length 2. Latitude range (min, max).
#'   e.g., `c(-80, 90)`.
#' @param depth Numeric vector of length 2. Depth range in meters (min, max).
#'   e.g., `c(0.51, 199.79)`.
#' @param out_name Character. Output file name ending with `.nc`.
#'
#' @return No return value. The function downloads a NetCDF file to `out_dir`.
#'
#' @examples
#' \donttest{
#' path_cms<-"C:/shark/copernicusmarine.exe" # replace to your own directory
#' out_dir <- "F:/whaleshark_sdm/cmems" # replace to your own directory
#' productId <- "cmems_mod_glo_phy-mnstd_my_0.25deg_P1M-m"
#' variable <- " --variable thetao_mean"
#' date_min <- lubridate::ymd(20110101)
#' date_max <- lubridate::ymd(20120101)
#' lon <- c(100, 130)
#' lat <- c(15, 30)
#' depth <- c(0.51, 50)
#' out_name <- "Pacific_temp_05_50_2011.nc"
#' env_download(path_cms,out_dir,productId,variable,
#'             date_min,date_max,lon,lat,
#'             depth,out_name)
#' }
#' @export
env_download<-function(path_cmems_tool,out_dir,productId,variable,
                       date_min,date_max,lon,lat,depth,out_name){

  # Ensure output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Build system command
  command <- paste(shQuote(path_cmems_tool), "subset -i", productId,
                   "-x", lon[1], "-X", lon[2],
                   "-y", lat[1], "-Y", lat[2],
                   "-t", date_min, "-T", date_max,
                   "-z", depth[1], "-Z", depth[2],
                   "-v", variable, "-o", shQuote(out_dir), "-f", out_name,
                    sep = " ")

  print(command)
  system(command)
}
