#' Create temporally binned predictor layer according to species traits
#'
#' @param sp_trait Character. Species vertical habitat use. One of `"pelagic"`,
#'   `"benthic"`, or `"benthopelagic"`.
#' @param env_path Character. Path to the downloaded CMEMS NetCDF file.
#' @param out_dir Character. Directory where the temporally‑binned raster layers
#'   will be saved.
#' @param var_name Character. Name of the variable (refer to CMEMS user manual).
#' @param lon_min Numeric. Minimum longitude (western boundary).
#' @param lon_max Numeric. Maximum longitude (eastern boundary).
#' @param lat_min Numeric. Minimum latitude (southern boundary).
#' @param lat_max Numeric. Maximum latitude (northern boundary).
#' @param date_min Date (or character convertible to Date). Start date.
#' @param date_max Date (or character convertible to Date). End date.
#' @param lag_days Numeric. Number of days ocean seasonality lags behind
#'   atmospheric seasonality. Default is `76` if missing.
#' @param temporal_bin Character. How a year will be split. One of `"monthly"`,
#'   `"quarterly"`, or `"seasonal"`.
#'
#' @return No return value. The function writes raster files (ASCII format) to
#'   the directory specified by `out_dir`. Files are named according to the
#'   temporal bin (e.g., `Jan_thetao_mean.asc`, `Q1_thetao_mean.asc`,
#'   `Spring_thetao_mean.asc`).
#'
#' @examples
#' \donttest{
#' env_bin(
#'   sp_trait = "pelagic",
#'   env_path = "F:/whaleshark_sdm/cmems/IndoWPac_temp_400m_y2015_y2024.nc",
#'   out_dir = "F:/whaleshark_sdm/bin",
#'   var_name = "thetao_mean",
#'   lon_min = 100, lon_max = 130,
#'   lat_min = 15, lat_max = 30,
#'   date_min = lubridate::ymd(19930101),
#'   date_max = lubridate::ymd(20221201),
#'   lag_days = 50,
#'   temporal_bin = "seasonal"
#' )
#' }
#' @export
env_bin <- function(sp_trait, env_path, out_dir, var_name,
                    lon_min, lon_max, lat_min, lat_max,
                    date_min, date_max, lag_days, temporal_bin) {

  # ---- parameter validation ----
  sp_trait <- match.arg(sp_trait, choices = c("pelagic", "benthic", "benthopelagic"))
  temporal_bin <- match.arg(temporal_bin, choices = c("monthly", "quarterly", "seasonal"))

  # default values (if missing)
  if (missing(lag_days)) lag_days <- 76
  if (missing(date_min)) date_min <- lubridate::ymd(19930101)
  if (missing(date_max)) date_max <- lubridate::ymd(20221201)
  if (missing(var_name)) var_name <- "thetao_mean"

  # ---- read NetCDF ----
  nc <- ncdf4::nc_open(env_path)

  # ---- helper: convert month indices to seasonal/quarterly groups ----
  # (the code uses month sequences, we keep the original approach)

  # ---- pelagic vs benthic/benthopelagic ----
  if (sp_trait == "pelagic") {
    # ----- pelagic: 3D or 4D data (depth dimension present) -----
    if (temporal_bin == "monthly") {
      inter <- strsplit(as.character(date_min), split = "-")
      inter2 <- strsplit(as.character(date_max), split = "-")
      min_month <- as.numeric(sapply(inter, function(x) x[2]))
      min_year  <- as.numeric(sapply(inter, function(x) x[1]))
      max_month <- as.numeric(sapply(inter2, function(x) x[2]))
      max_year  <- as.numeric(sapply(inter2, function(x) x[1]))
      length_month <- 12 * (max_year - min_year) + (max_month - min_month) + 1

      Jan <- seq(1, length_month - 11, by = 12)
      Feb <- seq(2, length_month - 10, by = 12)
      Mar <- seq(3, length_month - 9,  by = 12)
      Apr <- seq(4, length_month - 8,  by = 12)
      May <- seq(5, length_month - 7,  by = 12)
      Jun <- seq(6, length_month - 6,  by = 12)
      Jul <- seq(7, length_month - 5,  by = 12)
      Aug <- seq(8, length_month - 4,  by = 12)
      Sep <- seq(9, length_month - 3,  by = 12)
      Oct <- seq(10, length_month - 2, by = 12)
      Nov <- seq(11, length_month - 1, by = 12)
      Dec <- seq(12, length_month,      by = 12)

      for (i in 1:length_month) {
        # extract 2D slice (depth = 1) for pelagic
        if (length(dim(ncdf4::ncvar_get(nc, varid = var_name))) == 4) {
          predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , 1, i]
        } else {
          predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , i]
        }
        predictor_r1 <- raster::raster(t(predictor),
                                       xmn = lon_min, xmx = lon_max,
                                       ymn = lat_min, ymx = lat_max)
        predictor_r1 <- raster::flip(predictor_r1, direction = 2)
        raster::projection(predictor_r1) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0'

        # assign to month stacks
        if (i %in% Jan) {
          if (exists("Jan_predictor1")) Jan_predictor1 <- raster::stack(predictor_r1, Jan_predictor1)
          else Jan_predictor1 <- predictor_r1
        }
        if (i %in% Feb) {
          if (exists("Feb_predictor1")) Feb_predictor1 <- raster::stack(predictor_r1, Feb_predictor1)
          else Feb_predictor1 <- predictor_r1
        }
        if (i %in% Mar) {
          if (exists("Mar_predictor1")) Mar_predictor1 <- raster::stack(predictor_r1, Mar_predictor1)
          else Mar_predictor1 <- predictor_r1
        }
        if (i %in% Apr) {
          if (exists("Apr_predictor1")) Apr_predictor1 <- raster::stack(predictor_r1, Apr_predictor1)
          else Apr_predictor1 <- predictor_r1
        }
        if (i %in% May) {
          if (exists("May_predictor1")) May_predictor1 <- raster::stack(predictor_r1, May_predictor1)
          else May_predictor1 <- predictor_r1
        }
        if (i %in% Jun) {
          if (exists("Jun_predictor1")) Jun_predictor1 <- raster::stack(predictor_r1, Jun_predictor1)
          else Jun_predictor1 <- predictor_r1
        }
        if (i %in% Jul) {
          if (exists("Jul_predictor1")) Jul_predictor1 <- raster::stack(predictor_r1, Jul_predictor1)
          else Jul_predictor1 <- predictor_r1
        }
        if (i %in% Aug) {
          if (exists("Aug_predictor1")) Aug_predictor1 <- raster::stack(predictor_r1, Aug_predictor1)
          else Aug_predictor1 <- predictor_r1
        }
        if (i %in% Sep) {
          if (exists("Sep_predictor1")) Sep_predictor1 <- raster::stack(predictor_r1, Sep_predictor1)
          else Sep_predictor1 <- predictor_r1
        }
        if (i %in% Oct) {
          if (exists("Oct_predictor1")) Oct_predictor1 <- raster::stack(predictor_r1, Oct_predictor1)
          else Oct_predictor1 <- predictor_r1
        }
        if (i %in% Nov) {
          if (exists("Nov_predictor1")) Nov_predictor1 <- raster::stack(predictor_r1, Nov_predictor1)
          else Nov_predictor1 <- predictor_r1
        }
        if (i %in% Dec) {
          if (exists("Dec_predictor1")) Dec_predictor1 <- raster::stack(predictor_r1, Dec_predictor1)
          else Dec_predictor1 <- predictor_r1
        }
      }
      Jan_predictor <- raster::mean(Jan_predictor1)
      Feb_predictor <- raster::mean(Feb_predictor1)
      Mar_predictor <- raster::mean(Mar_predictor1)
      Apr_predictor <- raster::mean(Apr_predictor1)
      May_predictor <- raster::mean(May_predictor1)
      Jun_predictor <- raster::mean(Jun_predictor1)
      Jul_predictor <- raster::mean(Jul_predictor1)
      Aug_predictor <- raster::mean(Aug_predictor1)
      Sep_predictor <- raster::mean(Sep_predictor1)
      Oct_predictor <- raster::mean(Oct_predictor1)
      Nov_predictor <- raster::mean(Nov_predictor1)
      Dec_predictor <- raster::mean(Dec_predictor1)

      raster::writeRaster(Jan_predictor, paste0(out_dir,"/Jan_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Feb_predictor, paste0(out_dir,"/Feb_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Mar_predictor, paste0(out_dir,"/Mar_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Apr_predictor, paste0(out_dir,"/Apr_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(May_predictor, paste0(out_dir,"/May_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Jun_predictor, paste0(out_dir,"/Jun_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Jul_predictor, paste0(out_dir,"/Jul_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Aug_predictor, paste0(out_dir,"/Aug_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Sep_predictor, paste0(out_dir,"/Sep_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Oct_predictor, paste0(out_dir,"/Oct_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Nov_predictor, paste0(out_dir,"/Nov_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Dec_predictor, paste0(out_dir,"/Dec_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
    }  # end monthly

    if (temporal_bin == "quarterly") {
      inter <- strsplit(as.character(date_min), split = "-")
      inter2 <- strsplit(as.character(date_max), split = "-")
      min_month <- as.numeric(sapply(inter, function(x) x[2]))
      min_year  <- as.numeric(sapply(inter, function(x) x[1]))
      max_month <- as.numeric(sapply(inter2, function(x) x[2]))
      max_year  <- as.numeric(sapply(inter2, function(x) x[1]))
      length_month <- 12 * (max_year - min_year) + (max_month - min_month) + 1

      all_dates <- seq(from = as.Date(paste(min_year, min_month, "15", sep = "-")),
                       to   = as.Date(paste(max_year, max_month, "15", sep = "-")),
                       by = "month")
      months <- as.numeric(format(all_dates, "%m"))
      quarters <- ceiling(months / 3)
      row_numbers <- 1:length(months)
      Q1 <- row_numbers[quarters == 1]
      Q2 <- row_numbers[quarters == 2]
      Q3 <- row_numbers[quarters == 3]
      Q4 <- row_numbers[quarters == 4]

      for (i in 1:length_month) {
        if (length(dim(ncdf4::ncvar_get(nc, varid = var_name))) == 4) {
          predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , 1, i]
        } else {
          predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , i]
        }
        predictor_r1 <- raster::raster(t(predictor),
                                       xmn = lon_min, xmx = lon_max,
                                       ymn = lat_min, ymx = lat_max)
        predictor_r1 <- raster::flip(predictor_r1, direction = 2)
        raster::projection(predictor_r1) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0'

        if (i %in% Q1) {
          if (exists("Q1_predictor1")) Q1_predictor1 <- raster::stack(predictor_r1, Q1_predictor1)
          else Q1_predictor1 <- predictor_r1
        }
        if (i %in% Q2) {
          if (exists("Q2_predictor1")) Q2_predictor1 <- raster::stack(predictor_r1, Q2_predictor1)
          else Q2_predictor1 <- predictor_r1
        }
        if (i %in% Q3) {
          if (exists("Q3_predictor1")) Q3_predictor1 <- raster::stack(predictor_r1, Q3_predictor1)
          else Q3_predictor1 <- predictor_r1
        }
        if (i %in% Q4) {
          if (exists("Q4_predictor1")) Q4_predictor1 <- raster::stack(predictor_r1, Q4_predictor1)
          else Q4_predictor1 <- predictor_r1
        }
      }

      Q1_predictor <- raster::mean(Q1_predictor1)
      Q2_predictor <- raster::mean(Q2_predictor1)
      Q3_predictor <- raster::mean(Q3_predictor1)
      Q4_predictor <- raster::mean(Q4_predictor1)

      raster::writeRaster(Q1_predictor, paste0(out_dir,"/Q1_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Q2_predictor, paste0(out_dir,"/Q2_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Q3_predictor, paste0(out_dir,"/Q3_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Q4_predictor, paste0(out_dir,"/Q4_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
    }  # end quarterly

    if (temporal_bin == "seasonal") {
      inter <- strsplit(as.character(date_min), split = "-")
      inter2 <- strsplit(as.character(date_max), split = "-")
      min_month <- as.numeric(sapply(inter, function(x) x[2]))
      min_year  <- as.numeric(sapply(inter, function(x) x[1]))
      max_month <- as.numeric(sapply(inter2, function(x) x[2]))
      max_year  <- as.numeric(sapply(inter2, function(x) x[1]))
      length_month <- 12 * (max_year - min_year) + (max_month - min_month) + 1

      all_dates <- seq(from = as.Date(paste(min_year, min_month, "15", sep = "-")),
                       to   = as.Date(paste(max_year, max_month, "15", sep = "-")),
                       by = "month")
      adjusted_dates <- all_dates - lubridate::days(lag_days)
      adjusted_months <- lubridate::month(adjusted_dates)

      NH_spring <- which(adjusted_months %in% c(3,4,5))
      NH_summer <- which(adjusted_months %in% c(6,7,8))
      NH_autumn <- which(adjusted_months %in% c(9,10,11))
      NH_winter <- which(adjusted_months %in% c(12,1,2))

      SH_spring <- which(adjusted_months %in% c(9,10,11))
      SH_summer <- which(adjusted_months %in% c(12,1,2))
      SH_autumn <- which(adjusted_months %in% c(3,4,5))
      SH_winter <- which(adjusted_months %in% c(6,7,8))

      for (i in 1:length_month) {
        if (length(dim(ncdf4::ncvar_get(nc, varid = var_name))) == 4) {
          predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , 1, i]
        } else {
          predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , i]
        }
        predictor_r1 <- raster::raster(t(predictor),
                                       xmn = lon_min, xmx = lon_max,
                                       ymn = lat_min, ymx = lat_max)
        predictor_r1 <- raster::flip(predictor_r1, direction = 2)
        raster::projection(predictor_r1) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0'

        coords <- raster::coordinates(predictor_r1)
        lat_mask <- coords[,2]
        NH_raster <- predictor_r1
        NH_raster[lat_mask < 0] <- NA
        SH_raster <- predictor_r1
        SH_raster[lat_mask >= 0] <- NA

        if (i %in% NH_spring) {
          if (exists("NH_spring_stack")) NH_spring_stack <- raster::stack(NH_raster, NH_spring_stack)
          else NH_spring_stack <- NH_raster
        }
        if (i %in% NH_summer) {
          if (exists("NH_summer_stack")) NH_summer_stack <- raster::stack(NH_raster, NH_summer_stack)
          else NH_summer_stack <- NH_raster
        }
        if (i %in% NH_autumn) {
          if (exists("NH_autumn_stack")) NH_autumn_stack <- raster::stack(NH_raster, NH_autumn_stack)
          else NH_autumn_stack <- NH_raster
        }
        if (i %in% NH_winter) {
          if (exists("NH_winter_stack")) NH_winter_stack <- raster::stack(NH_raster, NH_winter_stack)
          else NH_winter_stack <- NH_raster
        }

        if (i %in% SH_spring) {
          if (exists("SH_spring_stack")) SH_spring_stack <- raster::stack(SH_raster, SH_spring_stack)
          else SH_spring_stack <- SH_raster
        }
        if (i %in% SH_summer) {
          if (exists("SH_summer_stack")) SH_summer_stack <- raster::stack(SH_raster, SH_summer_stack)
          else SH_summer_stack <- SH_raster
        }
        if (i %in% SH_autumn) {
          if (exists("SH_autumn_stack")) SH_autumn_stack <- raster::stack(SH_raster, SH_autumn_stack)
          else SH_autumn_stack <- SH_raster
        }
        if (i %in% SH_winter) {
          if (exists("SH_winter_stack")) SH_winter_stack <- raster::stack(SH_raster, SH_winter_stack)
          else SH_winter_stack <- SH_raster
        }
      }

      if (exists("NH_spring_stack") && exists("SH_spring_stack")) {
        spring_predictor <- raster::merge(raster::mean(NH_spring_stack, na.rm = TRUE),
                                          raster::mean(SH_spring_stack, na.rm = TRUE))
      } else if (exists("NH_spring_stack")) {
        spring_predictor <- raster::mean(NH_spring_stack, na.rm = TRUE)
      } else if (exists("SH_spring_stack")) {
        spring_predictor <- raster::mean(SH_spring_stack, na.rm = TRUE)
      }

      if (exists("NH_summer_stack") && exists("SH_summer_stack")) {
        summer_predictor <- raster::merge(raster::mean(NH_summer_stack, na.rm = TRUE),
                                          raster::mean(SH_summer_stack, na.rm = TRUE))
      } else if (exists("NH_summer_stack")) {
        summer_predictor <- raster::mean(NH_summer_stack, na.rm = TRUE)
      } else if (exists("SH_summer_stack")) {
        summer_predictor <- raster::mean(SH_summer_stack, na.rm = TRUE)
      }

      if (exists("NH_autumn_stack") && exists("SH_autumn_stack")) {
        autumn_predictor <- raster::merge(raster::mean(NH_autumn_stack, na.rm = TRUE),
                                          raster::mean(SH_autumn_stack, na.rm = TRUE))
      } else if (exists("NH_autumn_stack")) {
        autumn_predictor <- raster::mean(NH_autumn_stack, na.rm = TRUE)
      } else if (exists("SH_autumn_stack")) {
        autumn_predictor <- raster::mean(SH_autumn_stack, na.rm = TRUE)
      }

      if (exists("NH_winter_stack") && exists("SH_winter_stack")) {
        winter_predictor <- raster::merge(raster::mean(NH_winter_stack, na.rm = TRUE),
                                          raster::mean(SH_winter_stack, na.rm = TRUE))
      } else if (exists("NH_winter_stack")) {
        winter_predictor <- raster::mean(NH_winter_stack, na.rm = TRUE)
      } else if (exists("SH_winter_stack")) {
        winter_predictor <- raster::mean(SH_winter_stack, na.rm = TRUE)
      }

      raster::writeRaster(spring_predictor, paste0(out_dir,"/Spring_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(summer_predictor, paste0(out_dir,"/Summer_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(autumn_predictor, paste0(out_dir,"/Autumn_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(winter_predictor, paste0(out_dir,"/Winter_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
    }  # end seasonal
  }  # end pelagic

  # ----- benthic / benthopelagic (no depth dimension) -----
  if (sp_trait %in% c("benthic", "benthopelagic")) {
    if (temporal_bin == "monthly") {
      inter <- strsplit(as.character(date_min), split = "-")
      inter2 <- strsplit(as.character(date_max), split = "-")
      min_month <- as.numeric(sapply(inter, function(x) x[2]))
      min_year  <- as.numeric(sapply(inter, function(x) x[1]))
      max_month <- as.numeric(sapply(inter2, function(x) x[2]))
      max_year  <- as.numeric(sapply(inter2, function(x) x[1]))
      length_month <- 12 * (max_year - min_year) + (max_month - min_month) + 1

      Jan <- seq(1, length_month - 11, by = 12)
      Feb <- seq(2, length_month - 10, by = 12)
      Mar <- seq(3, length_month - 9,  by = 12)
      Apr <- seq(4, length_month - 8,  by = 12)
      May <- seq(5, length_month - 7,  by = 12)
      Jun <- seq(6, length_month - 6,  by = 12)
      Jul <- seq(7, length_month - 5,  by = 12)
      Aug <- seq(8, length_month - 4,  by = 12)
      Sep <- seq(9, length_month - 3,  by = 12)
      Oct <- seq(10, length_month - 2, by = 12)
      Nov <- seq(11, length_month - 1, by = 12)
      Dec <- seq(12, length_month,      by = 12)

      for (i in 1:length_month) {
        predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , i]
        predictor_r1 <- raster::raster(predictor,
                                       xmn = lon_min, xmx = lon_max,
                                       ymn = lat_min, ymx = lat_max)
        predictor_r1 <- raster::flip(predictor_r1, direction = 2)
        raster::projection(predictor_r1) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0'

        if (i %in% Jan) {
          if (exists("Jan_predictor1")) Jan_predictor1 <- raster::stack(predictor_r1, Jan_predictor1)
          else Jan_predictor1 <- predictor_r1
        }
        if (i %in% Feb) {
          if (exists("Feb_predictor1")) Feb_predictor1 <- raster::stack(predictor_r1, Feb_predictor1)
          else Feb_predictor1 <- predictor_r1
        }
        if (i %in% Mar) {
          if (exists("Mar_predictor1")) Mar_predictor1 <- raster::stack(predictor_r1, Mar_predictor1)
          else Mar_predictor1 <- predictor_r1
        }
        if (i %in% Apr) {
          if (exists("Apr_predictor1")) Apr_predictor1 <- raster::stack(predictor_r1, Apr_predictor1)
          else Apr_predictor1 <- predictor_r1
        }
        if (i %in% May) {
          if (exists("May_predictor1")) May_predictor1 <- raster::stack(predictor_r1, May_predictor1)
          else May_predictor1 <- predictor_r1
        }
        if (i %in% Jun) {
          if (exists("Jun_predictor1")) Jun_predictor1 <- raster::stack(predictor_r1, Jun_predictor1)
          else Jun_predictor1 <- predictor_r1
        }
        if (i %in% Jul) {
          if (exists("Jul_predictor1")) Jul_predictor1 <- raster::stack(predictor_r1, Jul_predictor1)
          else Jul_predictor1 <- predictor_r1
        }
        if (i %in% Aug) {
          if (exists("Aug_predictor1")) Aug_predictor1 <- raster::stack(predictor_r1, Aug_predictor1)
          else Aug_predictor1 <- predictor_r1
        }
        if (i %in% Sep) {
          if (exists("Sep_predictor1")) Sep_predictor1 <- raster::stack(predictor_r1, Sep_predictor1)
          else Sep_predictor1 <- predictor_r1
        }
        if (i %in% Oct) {
          if (exists("Oct_predictor1")) Oct_predictor1 <- raster::stack(predictor_r1, Oct_predictor1)
          else Oct_predictor1 <- predictor_r1
        }
        if (i %in% Nov) {
          if (exists("Nov_predictor1")) Nov_predictor1 <- raster::stack(predictor_r1, Nov_predictor1)
          else Nov_predictor1 <- predictor_r1
        }
        if (i %in% Dec) {
          if (exists("Dec_predictor1")) Dec_predictor1 <- raster::stack(predictor_r1, Dec_predictor1)
          else Dec_predictor1 <- predictor_r1
        }
      }
      Jan_predictor <- raster::mean(Jan_predictor1)
      Feb_predictor <- raster::mean(Feb_predictor1)
      Mar_predictor <- raster::mean(Mar_predictor1)
      Apr_predictor <- raster::mean(Apr_predictor1)
      May_predictor <- raster::mean(May_predictor1)
      Jun_predictor <- raster::mean(Jun_predictor1)
      Jul_predictor <- raster::mean(Jul_predictor1)
      Aug_predictor <- raster::mean(Aug_predictor1)
      Sep_predictor <- raster::mean(Sep_predictor1)
      Oct_predictor <- raster::mean(Oct_predictor1)
      Nov_predictor <- raster::mean(Nov_predictor1)
      Dec_predictor <- raster::mean(Dec_predictor1)

      raster::writeRaster(Jan_predictor, paste0(out_dir,"/Jan_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Feb_predictor, paste0(out_dir,"/Feb_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Mar_predictor, paste0(out_dir,"/Mar_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Apr_predictor, paste0(out_dir,"/Apr_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(May_predictor, paste0(out_dir,"/May_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Jun_predictor, paste0(out_dir,"/Jun_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Jul_predictor, paste0(out_dir,"/Jul_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Aug_predictor, paste0(out_dir,"/Aug_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Sep_predictor, paste0(out_dir,"/Sep_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Oct_predictor, paste0(out_dir,"/Oct_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Nov_predictor, paste0(out_dir,"/Nov_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Dec_predictor, paste0(out_dir,"/Dec_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)

    }  # end monthly (benthic)

    if (temporal_bin == "quarterly") {
      inter <- strsplit(as.character(date_min), split = "-")
      inter2 <- strsplit(as.character(date_max), split = "-")
      min_month <- as.numeric(sapply(inter, function(x) x[2]))
      min_year  <- as.numeric(sapply(inter, function(x) x[1]))
      max_month <- as.numeric(sapply(inter2, function(x) x[2]))
      max_year  <- as.numeric(sapply(inter2, function(x) x[1]))
      length_month <- 12 * (max_year - min_year) + (max_month - min_month) + 1

      all_dates <- seq(from = as.Date(paste(min_year, min_month, "15", sep = "-")),
                       to   = as.Date(paste(max_year, max_month, "15", sep = "-")),
                       by = "month")
      months <- as.numeric(format(all_dates, "%m"))
      quarters <- ceiling(months / 3)
      row_numbers <- 1:length(months)
      Q1 <- row_numbers[quarters == 1]
      Q2 <- row_numbers[quarters == 2]
      Q3 <- row_numbers[quarters == 3]
      Q4 <- row_numbers[quarters == 4]

      for (i in 1:length_month) {
        predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , i]
        predictor_r1 <- raster::raster(predictor,
                                       xmn = lon_min, xmx = lon_max,
                                       ymn = lat_min, ymx = lat_max)
        predictor_r1 <- raster::flip(predictor_r1, direction = 2)
        raster::projection(predictor_r1) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0'

        if (i %in% Q1) {
          if (exists("Q1_predictor1")) Q1_predictor1 <- raster::stack(predictor_r1, Q1_predictor1)
          else Q1_predictor1 <- predictor_r1
        }
        if (i %in% Q2) {
          if (exists("Q2_predictor1")) Q2_predictor1 <- raster::stack(predictor_r1, Q2_predictor1)
          else Q2_predictor1 <- predictor_r1
        }
        if (i %in% Q3) {
          if (exists("Q3_predictor1")) Q3_predictor1 <- raster::stack(predictor_r1, Q3_predictor1)
          else Q3_predictor1 <- predictor_r1
        }
        if (i %in% Q4) {
          if (exists("Q4_predictor1")) Q4_predictor1 <- raster::stack(predictor_r1, Q4_predictor1)
          else Q4_predictor1 <- predictor_r1
        }
      }
      Q1_predictor <- raster::mean(Q1_predictor1)
      Q2_predictor <- raster::mean(Q2_predictor1)
      Q3_predictor <- raster::mean(Q3_predictor1)
      Q4_predictor <- raster::mean(Q4_predictor1)

      raster::writeRaster(Q1_predictor, paste0(out_dir,"/Q1_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Q2_predictor, paste0(out_dir,"/Q2_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Q3_predictor, paste0(out_dir,"/Q3_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(Q4_predictor, paste0(out_dir,"/Q4_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
    }  # end quarterly (benthic)

    if (temporal_bin == "seasonal") {
      inter <- strsplit(as.character(date_min), split = "-")
      inter2 <- strsplit(as.character(date_max), split = "-")
      min_month <- as.numeric(sapply(inter, function(x) x[2]))
      min_year  <- as.numeric(sapply(inter, function(x) x[1]))
      max_month <- as.numeric(sapply(inter2, function(x) x[2]))
      max_year  <- as.numeric(sapply(inter2, function(x) x[1]))
      length_month <- 12 * (max_year - min_year) + (max_month - min_month) + 1

      all_dates <- seq(from = as.Date(paste(min_year, min_month, "15", sep = "-")),
                       to   = as.Date(paste(max_year, max_month, "15", sep = "-")),
                       by = "month")
      adjusted_dates <- all_dates - lubridate::days(lag_days)
      adjusted_months <- lubridate::month(adjusted_dates)

      NH_spring <- which(adjusted_months %in% c(3,4,5))
      NH_summer <- which(adjusted_months %in% c(6,7,8))
      NH_autumn <- which(adjusted_months %in% c(9,10,11))
      NH_winter <- which(adjusted_months %in% c(12,1,2))

      SH_spring <- which(adjusted_months %in% c(9,10,11))
      SH_summer <- which(adjusted_months %in% c(12,1,2))
      SH_autumn <- which(adjusted_months %in% c(3,4,5))
      SH_winter <- which(adjusted_months %in% c(6,7,8))

      for (i in 1:length_month) {
        predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , i]
        predictor_r1 <- raster::raster(predictor,
                                       xmn = lon_min, xmx = lon_max,
                                       ymn = lat_min, ymx = lat_max)
        predictor_r1 <- raster::flip(predictor_r1, direction = 2)
        raster::projection(predictor_r1) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0'

        coords <- raster::coordinates(predictor_r1)
        lat_mask <- coords[,2]
        NH_raster <- predictor_r1
        NH_raster[lat_mask < 0] <- NA
        SH_raster <- predictor_r1
        SH_raster[lat_mask >= 0] <- NA

        if (i %in% NH_spring) {
          if (exists("NH_spring_stack")) NH_spring_stack <- raster::stack(NH_raster, NH_spring_stack)
          else NH_spring_stack <- NH_raster
        }
        if (i %in% NH_summer) {
          if (exists("NH_summer_stack")) NH_summer_stack <- raster::stack(NH_raster, NH_summer_stack)
          else NH_summer_stack <- NH_raster
        }
        if (i %in% NH_autumn) {
          if (exists("NH_autumn_stack")) NH_autumn_stack <- raster::stack(NH_raster, NH_autumn_stack)
          else NH_autumn_stack <- NH_raster
        }
        if (i %in% NH_winter) {
          if (exists("NH_winter_stack")) NH_winter_stack <- raster::stack(NH_raster, NH_winter_stack)
          else NH_winter_stack <- NH_raster
        }

        if (i %in% SH_spring) {
          if (exists("SH_spring_stack")) SH_spring_stack <- raster::stack(SH_raster, SH_spring_stack)
          else SH_spring_stack <- SH_raster
        }
        if (i %in% SH_summer) {
          if (exists("SH_summer_stack")) SH_summer_stack <- raster::stack(SH_raster, SH_summer_stack)
          else SH_summer_stack <- SH_raster
        }
        if (i %in% SH_autumn) {
          if (exists("SH_autumn_stack")) SH_autumn_stack <- raster::stack(SH_raster, SH_autumn_stack)
          else SH_autumn_stack <- SH_raster
        }
        if (i %in% SH_winter) {
          if (exists("SH_winter_stack")) SH_winter_stack <- raster::stack(SH_raster, SH_winter_stack)
          else SH_winter_stack <- SH_raster
        }
      }



      if (exists("NH_spring_stack") && exists("SH_spring_stack")) {
        spring_predictor <- raster::merge(raster::mean(NH_spring_stack, na.rm = TRUE),
                                          raster::mean(SH_spring_stack, na.rm = TRUE))
      } else if (exists("NH_spring_stack")) {
        spring_predictor <- raster::mean(NH_spring_stack, na.rm = TRUE)
      } else if (exists("SH_spring_stack")) {
        spring_predictor <- raster::mean(SH_spring_stack, na.rm = TRUE)
      }

      if (exists("NH_summer_stack") && exists("SH_summer_stack")) {
        summer_predictor <- raster::merge(raster::mean(NH_summer_stack, na.rm = TRUE),
                                          raster::mean(SH_summer_stack, na.rm = TRUE))
      } else if (exists("NH_summer_stack")) {
        summer_predictor <- raster::mean(NH_summer_stack, na.rm = TRUE)
      } else if (exists("SH_summer_stack")) {
        summer_predictor <- raster::mean(SH_summer_stack, na.rm = TRUE)
      }

      if (exists("NH_autumn_stack") && exists("SH_autumn_stack")) {
        autumn_predictor <- raster::merge(raster::mean(NH_autumn_stack, na.rm = TRUE),
                                          raster::mean(SH_autumn_stack, na.rm = TRUE))
      } else if (exists("NH_autumn_stack")) {
        autumn_predictor <- raster::mean(NH_autumn_stack, na.rm = TRUE)
      } else if (exists("SH_autumn_stack")) {
        autumn_predictor <- raster::mean(SH_autumn_stack, na.rm = TRUE)
      }

      if (exists("NH_winter_stack") && exists("SH_winter_stack")) {
        winter_predictor <- raster::merge(raster::mean(NH_winter_stack, na.rm = TRUE),
                                          raster::mean(SH_winter_stack, na.rm = TRUE))
      } else if (exists("NH_winter_stack")) {
        winter_predictor <- raster::mean(NH_winter_stack, na.rm = TRUE)
      } else if (exists("SH_winter_stack")) {
        winter_predictor <- raster::mean(SH_winter_stack, na.rm = TRUE)
      }

      raster::writeRaster(spring_predictor, paste0(out_dir,"/Spring_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(summer_predictor, paste0(out_dir,"/Summer_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(autumn_predictor, paste0(out_dir,"/Autumn_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
      raster::writeRaster(winter_predictor, paste0(out_dir,"/Winter_", var_name, ".asc"),
                          format = "ascii", NAflag = -9999, overwrite = TRUE)
    }  # end seasonal (benthic)
  }  # end benthic/benthopelagic

  ncdf4::nc_close(nc)
}
