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
#'   atmospheric seasonality. e.g., `76` if missing.
#' @param temporal_bin Character. How a year will be split. One of `"monthly"`,
#'   `"quarterly"`, or `"seasonal"`.
#' @return No return value. The function writes raster files (ASCII format) to
#'   subdirectories of `out_dir` (e.g., `out_dir/Jan/thetao_mean.asc`,
#'   `out_dir/Q1/thetao_mean.asc`, `out_dir/Spring/thetao_mean.asc`).
#' @examples
#' \donttest{
#' ## sea surface height
# env_bin(sp_trait="pelagic",
#        env_path="F:/cmems/IndoWPac_slh_1m_y2015_y2024.nc",
#        out_dir="F:/whaleshark_sdm/bin",
#        var_name="zos_mean",
#        lon_min=90,
#        lon_max=180,
#        lat_min=-41,
#        lat_max=40,
#        date_min=ymd(20150101),
#        date_max = ymd(20241201),
#        temporal_bin="quarterly")
#' @export
env_bin <- function(sp_trait, env_path, out_dir, var_name,
                    lon_min, lon_max, lat_min, lat_max,
                    date_min, date_max, lag_days, temporal_bin) {

  cat("=== env_bin: Starting environmental binning ===\n")
  sp_trait <- match.arg(sp_trait, choices = c("pelagic", "benthic", "benthopelagic"))
  temporal_bin <- match.arg(temporal_bin, choices = c("monthly", "quarterly", "seasonal"))

  if (missing(lag_days)) lag_days <- 76
  if (missing(date_min)) date_min <- lubridate::ymd(19930101)
  if (missing(date_max)) date_max <- lubridate::ymd(20221201)
  if (missing(var_name)) var_name <- "thetao_mean"

  cat(sprintf(" Species trait: %s\n", sp_trait))
  cat(sprintf(" Temporal binning: %s\n", temporal_bin))
  cat(sprintf(" Variable: %s\n", var_name))
  cat(sprintf(" Date range: %s to %s\n", date_min, date_max))
  if (temporal_bin == "seasonal") cat(sprintf(" Lag days: %d\n", lag_days))
  cat(sprintf(" Output directory: %s\n", out_dir))

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Open NetCDF file
  cat("Opening NetCDF file...\n")
  nc <- ncdf4::nc_open(env_path)
  cat(sprintf(" Time steps: %d\n", nc$dim$time$len))
  if (!is.null(nc$dim$depth)) {
    cat(sprintf(" Depth levels: %d (using surface layer for pelagic)\n", nc$dim$depth$len))
  } else {
    cat(" No depth dimension (2D data)\n")
  }

  # Helper function to write raster (creates parent directory automatically)
  write_raster <- function(raster_obj, file_path) {
    dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
    terra::writeRaster(raster_obj, file_path,
                       filetype = "AAIGrid", NAflag = -9999,
                       overwrite = TRUE)
  }

  # ---- Main processing: pelagic vs benthic/benthopelagic ----
  if (sp_trait == "pelagic") {
    cat("Processing pelagic (using surface layer)\n")
    if (temporal_bin == "monthly") {
      cat("--- Monthly binning: computing monthly means ---\n")
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
      Dec <- seq(12, length_month,     by = 12)

      for (i in 1:length_month) {
        if (length(dim(ncdf4::ncvar_get(nc, varid = var_name))) == 4) {
          predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , 1, i]
        } else {
          predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , i]
        }
        predictor_r1 <- terra::rast(t(predictor))
        terra::ext(predictor_r1) <- c(lon_min, lon_max, lat_min, lat_max)
        terra::crs(predictor_r1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0"
        predictor_r1 <- terra::flip(predictor_r1, direction = "vertical")

        # Accumulate monthly stacks
        if (i %in% Jan) {
          if (exists("Jan_predictor1")) Jan_predictor1 <- c(Jan_predictor1, predictor_r1)
          else Jan_predictor1 <- predictor_r1
        }
        if (i %in% Feb) {
          if (exists("Feb_predictor1")) Feb_predictor1 <- c(Feb_predictor1, predictor_r1)
          else Feb_predictor1 <- predictor_r1
        }
        if (i %in% Mar) {
          if (exists("Mar_predictor1")) Mar_predictor1 <- c(Mar_predictor1, predictor_r1)
          else Mar_predictor1 <- predictor_r1
        }
        if (i %in% Apr) {
          if (exists("Apr_predictor1")) Apr_predictor1 <- c(Apr_predictor1, predictor_r1)
          else Apr_predictor1 <- predictor_r1
        }
        if (i %in% May) {
          if (exists("May_predictor1")) May_predictor1 <- c(May_predictor1, predictor_r1)
          else May_predictor1 <- predictor_r1
        }
        if (i %in% Jun) {
          if (exists("Jun_predictor1")) Jun_predictor1 <- c(Jun_predictor1, predictor_r1)
          else Jun_predictor1 <- predictor_r1
        }
        if (i %in% Jul) {
          if (exists("Jul_predictor1")) Jul_predictor1 <- c(Jul_predictor1, predictor_r1)
          else Jul_predictor1 <- predictor_r1
        }
        if (i %in% Aug) {
          if (exists("Aug_predictor1")) Aug_predictor1 <- c(Aug_predictor1, predictor_r1)
          else Aug_predictor1 <- predictor_r1
        }
        if (i %in% Sep) {
          if (exists("Sep_predictor1")) Sep_predictor1 <- c(Sep_predictor1, predictor_r1)
          else Sep_predictor1 <- predictor_r1
        }
        if (i %in% Oct) {
          if (exists("Oct_predictor1")) Oct_predictor1 <- c(Oct_predictor1, predictor_r1)
          else Oct_predictor1 <- predictor_r1
        }
        if (i %in% Nov) {
          if (exists("Nov_predictor1")) Nov_predictor1 <- c(Nov_predictor1, predictor_r1)
          else Nov_predictor1 <- predictor_r1
        }
        if (i %in% Dec) {
          if (exists("Dec_predictor1")) Dec_predictor1 <- c(Dec_predictor1, predictor_r1)
          else Dec_predictor1 <- predictor_r1
        }
      }

      # Compute monthly means
      Jan_predictor <- terra::mean(Jan_predictor1, na.rm = TRUE)
      Feb_predictor <- terra::mean(Feb_predictor1, na.rm = TRUE)
      Mar_predictor <- terra::mean(Mar_predictor1, na.rm = TRUE)
      Apr_predictor <- terra::mean(Apr_predictor1, na.rm = TRUE)
      May_predictor <- terra::mean(May_predictor1, na.rm = TRUE)
      Jun_predictor <- terra::mean(Jun_predictor1, na.rm = TRUE)
      Jul_predictor <- terra::mean(Jul_predictor1, na.rm = TRUE)
      Aug_predictor <- terra::mean(Aug_predictor1, na.rm = TRUE)
      Sep_predictor <- terra::mean(Sep_predictor1, na.rm = TRUE)
      Oct_predictor <- terra::mean(Oct_predictor1, na.rm = TRUE)
      Nov_predictor <- terra::mean(Nov_predictor1, na.rm = TRUE)
      Dec_predictor <- terra::mean(Dec_predictor1, na.rm = TRUE)

      # Plot (optional)
      terra::plot(Jan_predictor); terra::plot(Feb_predictor); terra::plot(Mar_predictor)
      terra::plot(Apr_predictor); terra::plot(May_predictor); terra::plot(Jun_predictor)
      terra::plot(Jul_predictor); terra::plot(Aug_predictor); terra::plot(Sep_predictor)
      terra::plot(Oct_predictor); terra::plot(Nov_predictor); terra::plot(Dec_predictor)

      # ===== MODIFIED: save each month into its own subfolder =====
      write_raster(Jan_predictor, file.path(out_dir, "Jan", paste0(var_name, ".asc")))
      write_raster(Feb_predictor, file.path(out_dir, "Feb", paste0(var_name, ".asc")))
      write_raster(Mar_predictor, file.path(out_dir, "Mar", paste0(var_name, ".asc")))
      write_raster(Apr_predictor, file.path(out_dir, "Apr", paste0(var_name, ".asc")))
      write_raster(May_predictor, file.path(out_dir, "May", paste0(var_name, ".asc")))
      write_raster(Jun_predictor, file.path(out_dir, "Jun", paste0(var_name, ".asc")))
      write_raster(Jul_predictor, file.path(out_dir, "Jul", paste0(var_name, ".asc")))
      write_raster(Aug_predictor, file.path(out_dir, "Aug", paste0(var_name, ".asc")))
      write_raster(Sep_predictor, file.path(out_dir, "Sep", paste0(var_name, ".asc")))
      write_raster(Oct_predictor, file.path(out_dir, "Oct", paste0(var_name, ".asc")))
      write_raster(Nov_predictor, file.path(out_dir, "Nov", paste0(var_name, ".asc")))
      write_raster(Dec_predictor, file.path(out_dir, "Dec", paste0(var_name, ".asc")))
      # ============================================================

      cat("--- Monthly binning completed. 12 files generated. ---\n")

    } else if (temporal_bin == "quarterly") {
      cat("--- Quarterly binning: computing quarterly means ---\n")
      inter <- strsplit(as.character(date_min), split = "-")
      inter2 <- strsplit(as.character(date_max), split = "-")
      min_month <- as.numeric(sapply(inter, function(x) x[2]))
      min_year  <- as.numeric(sapply(inter, function(x) x[1]))
      max_month <- as.numeric(sapply(inter2, function(x) x[2]))
      max_year  <- as.numeric(sapply(inter2, function(x) x[1]))
      length_month <- 12 * (max_year - min_year) + (max_month - min_month) + 1

      all_dates <- seq(from = as.Date(paste(min_year, min_month, "15", sep = "-")),
                       to   = as.Date(paste(max_year, max_month, "15", sep = "-")),
                       by   = "month")
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
        predictor_r1 <- terra::rast(t(predictor))
        terra::ext(predictor_r1) <- c(lon_min, lon_max, lat_min, lat_max)
        terra::crs(predictor_r1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0"
        predictor_r1 <- terra::flip(predictor_r1, direction = "vertical")

        if (i %in% Q1) {
          if (exists("Q1_predictor1")) Q1_predictor1 <- c(Q1_predictor1, predictor_r1)
          else Q1_predictor1 <- predictor_r1
        }
        if (i %in% Q2) {
          if (exists("Q2_predictor1")) Q2_predictor1 <- c(Q2_predictor1, predictor_r1)
          else Q2_predictor1 <- predictor_r1
        }
        if (i %in% Q3) {
          if (exists("Q3_predictor1")) Q3_predictor1 <- c(Q3_predictor1, predictor_r1)
          else Q3_predictor1 <- predictor_r1
        }
        if (i %in% Q4) {
          if (exists("Q4_predictor1")) Q4_predictor1 <- c(Q4_predictor1, predictor_r1)
          else Q4_predictor1 <- predictor_r1
        }
      }

      Q1_predictor <- terra::mean(Q1_predictor1, na.rm = TRUE)
      Q2_predictor <- terra::mean(Q2_predictor1, na.rm = TRUE)
      Q3_predictor <- terra::mean(Q3_predictor1, na.rm = TRUE)
      Q4_predictor <- terra::mean(Q4_predictor1, na.rm = TRUE)

      terra::plot(Q1_predictor); terra::plot(Q2_predictor)
      terra::plot(Q3_predictor); terra::plot(Q4_predictor)

      # ===== MODIFIED: save each quarter into its own subfolder =====
      write_raster(Q1_predictor, file.path(out_dir, "Q1", paste0(var_name, ".asc")))
      write_raster(Q2_predictor, file.path(out_dir, "Q2", paste0(var_name, ".asc")))
      write_raster(Q3_predictor, file.path(out_dir, "Q3", paste0(var_name, ".asc")))
      write_raster(Q4_predictor, file.path(out_dir, "Q4", paste0(var_name, ".asc")))
      # =============================================================

      cat("--- Quarterly binning completed. 4 files generated. ---\n")

    } else if (temporal_bin == "seasonal") {
      cat("--- Seasonal binning: computing seasonal means (Spring, Summer, Autumn, Winter) ---\n")
      cat(" Ocean seasons lag atmospheric seasons by", lag_days, "days\n")

      inter <- strsplit(as.character(date_min), split = "-")
      inter2 <- strsplit(as.character(date_max), split = "-")
      min_month <- as.numeric(sapply(inter, function(x) x[2]))
      min_year  <- as.numeric(sapply(inter, function(x) x[1]))
      max_month <- as.numeric(sapply(inter2, function(x) x[2]))
      max_year  <- as.numeric(sapply(inter2, function(x) x[1]))
      length_month <- 12 * (max_year - min_year) + (max_month - min_month) + 1

      all_dates <- seq(from = as.Date(paste(min_year, min_month, "15", sep = "-")),
                       to   = as.Date(paste(max_year, max_month, "15", sep = "-")),
                       by   = "month")
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
        predictor_r1 <- terra::rast(t(predictor))
        terra::ext(predictor_r1) <- c(lon_min, lon_max, lat_min, lat_max)
        terra::crs(predictor_r1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0"
        predictor_r1 <- terra::flip(predictor_r1, direction = "vertical")

        # Split hemispheres
        lat_rast <- terra::init(predictor_r1, fun = "y")
        NH_raster <- predictor_r1; NH_raster[lat_rast < 0] <- NA
        SH_raster <- predictor_r1; SH_raster[lat_rast >= 0] <- NA

        if (i %in% NH_spring) {
          if (exists("NH_spring_stack")) NH_spring_stack <- c(NH_spring_stack, NH_raster)
          else NH_spring_stack <- NH_raster
        }
        if (i %in% NH_summer) {
          if (exists("NH_summer_stack")) NH_summer_stack <- c(NH_summer_stack, NH_raster)
          else NH_summer_stack <- NH_raster
        }
        if (i %in% NH_autumn) {
          if (exists("NH_autumn_stack")) NH_autumn_stack <- c(NH_autumn_stack, NH_raster)
          else NH_autumn_stack <- NH_raster
        }
        if (i %in% NH_winter) {
          if (exists("NH_winter_stack")) NH_winter_stack <- c(NH_winter_stack, NH_raster)
          else NH_winter_stack <- NH_raster
        }
        if (i %in% SH_spring) {
          if (exists("SH_spring_stack")) SH_spring_stack <- c(SH_spring_stack, SH_raster)
          else SH_spring_stack <- SH_raster
        }
        if (i %in% SH_summer) {
          if (exists("SH_summer_stack")) SH_summer_stack <- c(SH_summer_stack, SH_raster)
          else SH_summer_stack <- SH_raster
        }
        if (i %in% SH_autumn) {
          if (exists("SH_autumn_stack")) SH_autumn_stack <- c(SH_autumn_stack, SH_raster)
          else SH_autumn_stack <- SH_raster
        }
        if (i %in% SH_winter) {
          if (exists("SH_winter_stack")) SH_winter_stack <- c(SH_winter_stack, SH_raster)
          else SH_winter_stack <- SH_raster
        }
      }

      # Compute seasonal means (merge hemispheres)
      if (exists("NH_spring_stack") && exists("SH_spring_stack")) {
        spring_predictor <- terra::merge(
          terra::mean(NH_spring_stack, na.rm = TRUE),
          terra::mean(SH_spring_stack, na.rm = TRUE)
        )
      } else if (exists("NH_spring_stack")) {
        spring_predictor <- terra::mean(NH_spring_stack, na.rm = TRUE)
      } else if (exists("SH_spring_stack")) {
        spring_predictor <- terra::mean(SH_spring_stack, na.rm = TRUE)
      }

      if (exists("NH_summer_stack") && exists("SH_summer_stack")) {
        summer_predictor <- terra::merge(
          terra::mean(NH_summer_stack, na.rm = TRUE),
          terra::mean(SH_summer_stack, na.rm = TRUE)
        )
      } else if (exists("NH_summer_stack")) {
        summer_predictor <- terra::mean(NH_summer_stack, na.rm = TRUE)
      } else if (exists("SH_summer_stack")) {
        summer_predictor <- terra::mean(SH_summer_stack, na.rm = TRUE)
      }

      if (exists("NH_autumn_stack") && exists("SH_autumn_stack")) {
        autumn_predictor <- terra::merge(
          terra::mean(NH_autumn_stack, na.rm = TRUE),
          terra::mean(SH_autumn_stack, na.rm = TRUE)
        )
      } else if (exists("NH_autumn_stack")) {
        autumn_predictor <- terra::mean(NH_autumn_stack, na.rm = TRUE)
      } else if (exists("SH_autumn_stack")) {
        autumn_predictor <- terra::mean(SH_autumn_stack, na.rm = TRUE)
      }

      if (exists("NH_winter_stack") && exists("SH_winter_stack")) {
        winter_predictor <- terra::merge(
          terra::mean(NH_winter_stack, na.rm = TRUE),
          terra::mean(SH_winter_stack, na.rm = TRUE)
        )
      } else if (exists("NH_winter_stack")) {
        winter_predictor <- terra::mean(NH_winter_stack, na.rm = TRUE)
      } else if (exists("SH_winter_stack")) {
        winter_predictor <- terra::mean(SH_winter_stack, na.rm = TRUE)
      }

      terra::plot(spring_predictor); terra::plot(summer_predictor)
      terra::plot(autumn_predictor); terra::plot(winter_predictor)

      # ===== MODIFIED: save each season into its own subfolder =====
      write_raster(spring_predictor, file.path(out_dir, "Spring", paste0(var_name, ".asc")))
      write_raster(summer_predictor, file.path(out_dir, "Summer", paste0(var_name, ".asc")))
      write_raster(autumn_predictor, file.path(out_dir, "Autumn", paste0(var_name, ".asc")))
      write_raster(winter_predictor, file.path(out_dir, "Winter", paste0(var_name, ".asc")))
      # =============================================================

      cat("--- Seasonal binning completed. 4 files generated. ---\n")
    }

  } else {
    # ---- Benthic / Benthopelagic: no depth dimension ----
    cat("Processing benthic/benthopelagic (2D data)\n")
    if (temporal_bin == "monthly") {
      cat("--- Monthly binning: computing monthly means ---\n")
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
      Dec <- seq(12, length_month,     by = 12)

      for (i in 1:length_month) {
        predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , i]
        predictor_r1 <- terra::rast(t(predictor))
        terra::ext(predictor_r1) <- c(lon_min, lon_max, lat_min, lat_max)
        terra::crs(predictor_r1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0"
        predictor_r1 <- terra::flip(predictor_r1, direction = "vertical")

        if (i %in% Jan) {
          if (exists("Jan_predictor1")) Jan_predictor1 <- c(Jan_predictor1, predictor_r1)
          else Jan_predictor1 <- predictor_r1
        }
        if (i %in% Feb) {
          if (exists("Feb_predictor1")) Feb_predictor1 <- c(Feb_predictor1, predictor_r1)
          else Feb_predictor1 <- predictor_r1
        }
        if (i %in% Mar) {
          if (exists("Mar_predictor1")) Mar_predictor1 <- c(Mar_predictor1, predictor_r1)
          else Mar_predictor1 <- predictor_r1
        }
        if (i %in% Apr) {
          if (exists("Apr_predictor1")) Apr_predictor1 <- c(Apr_predictor1, predictor_r1)
          else Apr_predictor1 <- predictor_r1
        }
        if (i %in% May) {
          if (exists("May_predictor1")) May_predictor1 <- c(May_predictor1, predictor_r1)
          else May_predictor1 <- predictor_r1
        }
        if (i %in% Jun) {
          if (exists("Jun_predictor1")) Jun_predictor1 <- c(Jun_predictor1, predictor_r1)
          else Jun_predictor1 <- predictor_r1
        }
        if (i %in% Jul) {
          if (exists("Jul_predictor1")) Jul_predictor1 <- c(Jul_predictor1, predictor_r1)
          else Jul_predictor1 <- predictor_r1
        }
        if (i %in% Aug) {
          if (exists("Aug_predictor1")) Aug_predictor1 <- c(Aug_predictor1, predictor_r1)
          else Aug_predictor1 <- predictor_r1
        }
        if (i %in% Sep) {
          if (exists("Sep_predictor1")) Sep_predictor1 <- c(Sep_predictor1, predictor_r1)
          else Sep_predictor1 <- predictor_r1
        }
        if (i %in% Oct) {
          if (exists("Oct_predictor1")) Oct_predictor1 <- c(Oct_predictor1, predictor_r1)
          else Oct_predictor1 <- predictor_r1
        }
        if (i %in% Nov) {
          if (exists("Nov_predictor1")) Nov_predictor1 <- c(Nov_predictor1, predictor_r1)
          else Nov_predictor1 <- predictor_r1
        }
        if (i %in% Dec) {
          if (exists("Dec_predictor1")) Dec_predictor1 <- c(Dec_predictor1, predictor_r1)
          else Dec_predictor1 <- predictor_r1
        }
      }

      Jan_predictor <- terra::mean(Jan_predictor1, na.rm = TRUE)
      Feb_predictor <- terra::mean(Feb_predictor1, na.rm = TRUE)
      Mar_predictor <- terra::mean(Mar_predictor1, na.rm = TRUE)
      Apr_predictor <- terra::mean(Apr_predictor1, na.rm = TRUE)
      May_predictor <- terra::mean(May_predictor1, na.rm = TRUE)
      Jun_predictor <- terra::mean(Jun_predictor1, na.rm = TRUE)
      Jul_predictor <- terra::mean(Jul_predictor1, na.rm = TRUE)
      Aug_predictor <- terra::mean(Aug_predictor1, na.rm = TRUE)
      Sep_predictor <- terra::mean(Sep_predictor1, na.rm = TRUE)
      Oct_predictor <- terra::mean(Oct_predictor1, na.rm = TRUE)
      Nov_predictor <- terra::mean(Nov_predictor1, na.rm = TRUE)
      Dec_predictor <- terra::mean(Dec_predictor1, na.rm = TRUE)

      terra::plot(Jan_predictor); terra::plot(Feb_predictor); terra::plot(Mar_predictor)
      terra::plot(Apr_predictor); terra::plot(May_predictor); terra::plot(Jun_predictor)
      terra::plot(Jul_predictor); terra::plot(Aug_predictor); terra::plot(Sep_predictor)
      terra::plot(Oct_predictor); terra::plot(Nov_predictor); terra::plot(Dec_predictor)

      # ===== MODIFIED: save each month into its own subfolder =====
      write_raster(Jan_predictor, file.path(out_dir, "Jan", paste0(var_name, ".asc")))
      write_raster(Feb_predictor, file.path(out_dir, "Feb", paste0(var_name, ".asc")))
      write_raster(Mar_predictor, file.path(out_dir, "Mar", paste0(var_name, ".asc")))
      write_raster(Apr_predictor, file.path(out_dir, "Apr", paste0(var_name, ".asc")))
      write_raster(May_predictor, file.path(out_dir, "May", paste0(var_name, ".asc")))
      write_raster(Jun_predictor, file.path(out_dir, "Jun", paste0(var_name, ".asc")))
      write_raster(Jul_predictor, file.path(out_dir, "Jul", paste0(var_name, ".asc")))
      write_raster(Aug_predictor, file.path(out_dir, "Aug", paste0(var_name, ".asc")))
      write_raster(Sep_predictor, file.path(out_dir, "Sep", paste0(var_name, ".asc")))
      write_raster(Oct_predictor, file.path(out_dir, "Oct", paste0(var_name, ".asc")))
      write_raster(Nov_predictor, file.path(out_dir, "Nov", paste0(var_name, ".asc")))
      write_raster(Dec_predictor, file.path(out_dir, "Dec", paste0(var_name, ".asc")))
      # ============================================================

      cat("--- Monthly binning completed. 12 files generated. ---\n")

    } else if (temporal_bin == "quarterly") {
      cat("--- Quarterly binning: computing quarterly means ---\n")
      inter <- strsplit(as.character(date_min), split = "-")
      inter2 <- strsplit(as.character(date_max), split = "-")
      min_month <- as.numeric(sapply(inter, function(x) x[2]))
      min_year  <- as.numeric(sapply(inter, function(x) x[1]))
      max_month <- as.numeric(sapply(inter2, function(x) x[2]))
      max_year  <- as.numeric(sapply(inter2, function(x) x[1]))
      length_month <- 12 * (max_year - min_year) + (max_month - min_month) + 1

      all_dates <- seq(from = as.Date(paste(min_year, min_month, "15", sep = "-")),
                       to   = as.Date(paste(max_year, max_month, "15", sep = "-")),
                       by   = "month")
      months <- as.numeric(format(all_dates, "%m"))
      quarters <- ceiling(months / 3)
      row_numbers <- 1:length(months)
      Q1 <- row_numbers[quarters == 1]
      Q2 <- row_numbers[quarters == 2]
      Q3 <- row_numbers[quarters == 3]
      Q4 <- row_numbers[quarters == 4]

      for (i in 1:length_month) {
        predictor <- ncdf4::ncvar_get(nc, varid = var_name)[, , i]
        predictor_r1 <- terra::rast(t(predictor))
        terra::ext(predictor_r1) <- c(lon_min, lon_max, lat_min, lat_max)
        terra::crs(predictor_r1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0"
        predictor_r1 <- terra::flip(predictor_r1, direction = "vertical")

        if (i %in% Q1) {
          if (exists("Q1_predictor1")) Q1_predictor1 <- c(Q1_predictor1, predictor_r1)
          else Q1_predictor1 <- predictor_r1
        }
        if (i %in% Q2) {
          if (exists("Q2_predictor1")) Q2_predictor1 <- c(Q2_predictor1, predictor_r1)
          else Q2_predictor1 <- predictor_r1
        }
        if (i %in% Q3) {
          if (exists("Q3_predictor1")) Q3_predictor1 <- c(Q3_predictor1, predictor_r1)
          else Q3_predictor1 <- predictor_r1
        }
        if (i %in% Q4) {
          if (exists("Q4_predictor1")) Q4_predictor1 <- c(Q4_predictor1, predictor_r1)
          else Q4_predictor1 <- predictor_r1
        }
      }

      Q1_predictor <- terra::mean(Q1_predictor1, na.rm = TRUE)
      Q2_predictor <- terra::mean(Q2_predictor1, na.rm = TRUE)
      Q3_predictor <- terra::mean(Q3_predictor1, na.rm = TRUE)
      Q4_predictor <- terra::mean(Q4_predictor1, na.rm = TRUE)

      terra::plot(Q1_predictor); terra::plot(Q2_predictor)
      terra::plot(Q3_predictor); terra::plot(Q4_predictor)

      # ===== MODIFIED: save each quarter into its own subfolder =====
      write_raster(Q1_predictor, file.path(out_dir, "Q1", paste0(var_name, ".asc")))
      write_raster(Q2_predictor, file.path(out_dir, "Q2", paste0(var_name, ".asc")))
      write_raster(Q3_predictor, file.path(out_dir, "Q3", paste0(var_name, ".asc")))
      write_raster(Q4_predictor, file.path(out_dir, "Q4", paste0(var_name, ".asc")))
      # =============================================================

      cat("--- Quarterly binning completed. 4 files generated. ---\n")

    } else if (temporal_bin == "seasonal") {
      cat("--- Seasonal binning: computing seasonal means (Spring, Summer, Autumn, Winter) ---\n")
      cat(" Ocean seasons lag atmospheric seasons by", lag_days, "days\n")

      inter <- strsplit(as.character(date_min), split = "-")
      inter2 <- strsplit(as.character(date_max), split = "-")
      min_month <- as.numeric(sapply(inter, function(x) x[2]))
      min_year  <- as.numeric(sapply(inter, function(x) x[1]))
      max_month <- as.numeric(sapply(inter2, function(x) x[2]))
      max_year  <- as.numeric(sapply(inter2, function(x) x[1]))
      length_month <- 12 * (max_year - min_year) + (max_month - min_month) + 1

      all_dates <- seq(from = as.Date(paste(min_year, min_month, "15", sep = "-")),
                       to   = as.Date(paste(max_year, max_month, "15", sep = "-")),
                       by   = "month")
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
        predictor_r1 <- terra::rast(t(predictor))
        terra::ext(predictor_r1) <- c(lon_min, lon_max, lat_min, lat_max)
        terra::crs(predictor_r1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0"
        predictor_r1 <- terra::flip(predictor_r1, direction = "vertical")

        lat_rast <- terra::init(predictor_r1, fun = "y")
        NH_raster <- predictor_r1; NH_raster[lat_rast < 0] <- NA
        SH_raster <- predictor_r1; SH_raster[lat_rast >= 0] <- NA

        if (i %in% NH_spring) {
          if (exists("NH_spring_stack")) NH_spring_stack <- c(NH_spring_stack, NH_raster)
          else NH_spring_stack <- NH_raster
        }
        if (i %in% NH_summer) {
          if (exists("NH_summer_stack")) NH_summer_stack <- c(NH_summer_stack, NH_raster)
          else NH_summer_stack <- NH_raster
        }
        if (i %in% NH_autumn) {
          if (exists("NH_autumn_stack")) NH_autumn_stack <- c(NH_autumn_stack, NH_raster)
          else NH_autumn_stack <- NH_raster
        }
        if (i %in% NH_winter) {
          if (exists("NH_winter_stack")) NH_winter_stack <- c(NH_winter_stack, NH_raster)
          else NH_winter_stack <- NH_raster
        }
        if (i %in% SH_spring) {
          if (exists("SH_spring_stack")) SH_spring_stack <- c(SH_spring_stack, SH_raster)
          else SH_spring_stack <- SH_raster
        }
        if (i %in% SH_summer) {
          if (exists("SH_summer_stack")) SH_summer_stack <- c(SH_summer_stack, SH_raster)
          else SH_summer_stack <- SH_raster
        }
        if (i %in% SH_autumn) {
          if (exists("SH_autumn_stack")) SH_autumn_stack <- c(SH_autumn_stack, SH_raster)
          else SH_autumn_stack <- SH_raster
        }
        if (i %in% SH_winter) {
          if (exists("SH_winter_stack")) SH_winter_stack <- c(SH_winter_stack, SH_raster)
          else SH_winter_stack <- SH_raster
        }
      }

      if (exists("NH_spring_stack") && exists("SH_spring_stack")) {
        spring_predictor <- terra::merge(
          terra::mean(NH_spring_stack, na.rm = TRUE),
          terra::mean(SH_spring_stack, na.rm = TRUE)
        )
      } else if (exists("NH_spring_stack")) {
        spring_predictor <- terra::mean(NH_spring_stack, na.rm = TRUE)
      } else if (exists("SH_spring_stack")) {
        spring_predictor <- terra::mean(SH_spring_stack, na.rm = TRUE)
      }

      if (exists("NH_summer_stack") && exists("SH_summer_stack")) {
        summer_predictor <- terra::merge(
          terra::mean(NH_summer_stack, na.rm = TRUE),
          terra::mean(SH_summer_stack, na.rm = TRUE)
        )
      } else if (exists("NH_summer_stack")) {
        summer_predictor <- terra::mean(NH_summer_stack, na.rm = TRUE)
      } else if (exists("SH_summer_stack")) {
        summer_predictor <- terra::mean(SH_summer_stack, na.rm = TRUE)
      }

      if (exists("NH_autumn_stack") && exists("SH_autumn_stack")) {
        autumn_predictor <- terra::merge(
          terra::mean(NH_autumn_stack, na.rm = TRUE),
          terra::mean(SH_autumn_stack, na.rm = TRUE)
        )
      } else if (exists("NH_autumn_stack")) {
        autumn_predictor <- terra::mean(NH_autumn_stack, na.rm = TRUE)
      } else if (exists("SH_autumn_stack")) {
        autumn_predictor <- terra::mean(SH_autumn_stack, na.rm = TRUE)
      }

      if (exists("NH_winter_stack") && exists("SH_winter_stack")) {
        winter_predictor <- terra::merge(
          terra::mean(NH_winter_stack, na.rm = TRUE),
          terra::mean(SH_winter_stack, na.rm = TRUE)
        )
      } else if (exists("NH_winter_stack")) {
        winter_predictor <- terra::mean(NH_winter_stack, na.rm = TRUE)
      } else if (exists("SH_winter_stack")) {
        winter_predictor <- terra::mean(SH_winter_stack, na.rm = TRUE)
      }

      terra::plot(spring_predictor); terra::plot(summer_predictor)
      terra::plot(autumn_predictor); terra::plot(winter_predictor)

      # ===== MODIFIED: save each season into its own subfolder =====
      write_raster(spring_predictor, file.path(out_dir, "Spring", paste0(var_name, ".asc")))
      write_raster(summer_predictor, file.path(out_dir, "Summer", paste0(var_name, ".asc")))
      write_raster(autumn_predictor, file.path(out_dir, "Autumn", paste0(var_name, ".asc")))
      write_raster(winter_predictor, file.path(out_dir, "Winter", paste0(var_name, ".asc")))
      # =============================================================

      cat("--- Seasonal binning completed. 4 files generated. ---\n")
    }
  }

  ncdf4::nc_close(nc)
  cat("=== env_bin: Processing completed successfully ===\n")
}
