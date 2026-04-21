#' Download environmental layers from CMEMS
#'
#' @param path_cmems_tool Character. Path to the copernicusmarine executable,
#'   e.g. `"/miniconda/envs/copernicusmarine/Scripts/copernicusmarine.exe"`.
#' @param out_dir Character. Directory where output data will be stored.
#' @param productId Character. ID of the CMEMS product (see dataset user manual
#'   in the CMEMS product webpage).
#'   Default is `"cmems_mod_glo_phy-mnstd_my_0.25deg_P1M-m"`.
#' @param variable Character. Variable name to download, prefixed with "--variable ".
#'   Default is " --variable thetao_mean".
#' @param date_min Date or character. Minimum date. Default is `1993-01-01`.
#' @param date_max Date or character. Maximum date. Default is `2022-12-01`.
#' @param lon Numeric vector of length 2. Longitude range (min, max).
#'   Default is `c(-180, 179.75)`.
#' @param lat Numeric vector of length 2. Latitude range (min, max).
#'   Default is `c(-80, 90)`.
#' @param depth List of length 2. Depth range in meters (min, max).
#'   Default is `list(0.51, 199.79)`.
#' @param out_name Character. Output file name ending with `.nc`.
#'
#' @return No return value. The function downloads a NetCDF file to `out_dir`.
#'
#' @examples
#' \donttest{
#' path_cmems_tool <- "C:/shark/miniconda/envs/copernicusmarine/Scripts/copernicusmarine.exe"
#' out_dir <- "C:/shark/cmems"
#' productId <- "cmems_mod_glo_phy-mnstd_my_0.25deg_P1M-m"
#' variable <- " --variable thetao_mean"
#' date_min <- lubridate::ymd(20110101)
#' date_max <- lubridate::ymd(20120101)
#' lon <- c(100, 130)
#' lat <- c(15, 30)
#' depth <- list(0.51, 50)
#' out_name <- "Pacific_temp_05_50_2011.nc"
#' env_download(path_cmems_tool, out_dir, productId, variable,
#'              date_min, date_max, lon, lat, depth, out_name)
#' }
#' @export
env_download<-function(path_cmems_tool,out_dir,productId,variable,
                       date_min,date_max,lon,lat,depth,out_name){
  #### Step 1: install Anaconda ######
  # install.packages("reticulate")
  # library(reticulate)
  # install_miniconda(path = "C:/shark/miniconda", update = T) # give your own path where miniconda will be installed
  # conda_list(conda = "C:/shark/miniconda/_conda.exe")
  # use_condaenv(condaenv = "r-reticulate", conda = "C:/shark/miniconda/_conda.exe")
  #### Step 2: install copernicusmarine and login #####
  # go to the search pannel of your PC and find  Anaconda prompt terminal and follow this to install copernicusmarine:  https://www.bilibili.com/read/cv28415761/
  # python -m pip install copernicusmarine
  # conda create --name copernicusmarine conda-forge::copernicusmarine --yes
  # conda activate copernicusmarine
  # copernicusmarine --version && copernicusmarine --help
  # set COPERNICUSMARINE_SERVICE_USERNAME=your_username # set up username
  # set COPERNICUSMARINE_SERVICE_PASSWORD=your_password # set up password
  # copernicusmarine login
  # Default values
  if (missing(productId)) {
    productId <- "cmems_mod_glo_phy-mnstd_my_0.25deg_P1M-m"
  }
  if (missing(date_min)) {
    date_min <- lubridate::ymd(19930101)
  }
  if (missing(date_max)) {
    date_max <- lubridate::ymd(20231201)
  }
  if (missing(lon)) {
    lon <- c(-180, 179.75)
  }
  if (missing(lat)) {
    lat <- c(-80, 90)
  }
  if (missing(depth)) {
    depth <- list(0.51, 199.79)
  }
  if (missing(variable)) {
    variable <- " --variable thetao_mean"
  }

  # Ensure output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Build system command
  command <- paste(path_cmems_tool, "subset -i", productId,
                   "-x", lon[1], "-X", lon[2],
                   "-y", lat[1], "-Y", lat[2],
                   "-t", date_min, "-T", date_max,
                   "-z", depth[1], "-Z", depth[2],
                   variable, "-o", out_dir, "-f", out_name,
                   "--force-download", sep = " ")

  print(command)
  system(command)
}
