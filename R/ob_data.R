#' Download & clean OBIS occurrences
#'
#' @param key Character. The species scientific name (e.g., `"Rhincodon typus"`).
#' @param sp_geometry Character. The geographic range of the study area as a WKT polygon
#'   (e.g., `"POLYGON ((20 -41, 20 40, 180 40, 180 -41, 20 -41))"`).
#' @param start_date Character. Start date used in `robis::occurrence()`, format `"YYYY-MM-DD"`.
#' @param end_date Character. End date used in `robis::occurrence()`, format `"YYYY-MM-DD"`.
#' @param event_date Character. Year range used in `rgbif::occ_search()`, format
#'   `"YYYY-MM-DD,YYYY-MM-DD"`.
#'
#' @return A list containing:
#'   \item{depth_summary}{Summary statistics of depth values.}
#'   \item{sp_occ}{Data frame of cleaned occurrence records with columns:
#'     `Species`, `decimalLatitude`, `decimalLongitude`, `depth`, `date`, `year`, `month`, `day`.}
#'   \item{plot}{A map plot of occurrences (returned by `obistools::plot_map`).}
#'
#'
#' @examples
#' \donttest{
#' ws_occ <- ob_data(
#'   key = "Rhincodon typus",
#'   sp_geometry =  "POLYGON ((90 -41, 90 40, 180 40, 180 -41, 90 -41))",
#'   start_date = "2010-01-01",
#'   end_date = "2025-12-30",
#'   event_date = "2010-01-01,2025-12-30"
#' )
#' }
#' @export
ob_data <- function(key, sp_geometry, start_date, end_date, event_date) {
  # ---- Download from OBIS via robis ----
  df <- robis::occurrence(
    scientificname = key,
    geometry = sp_geometry,
    startdate = start_date,
    enddate = end_date
  )
  ob0 <- as.data.frame(df)
  ob1 <- subset(ob0, !(is.na(decimalLatitude) | is.na(decimalLongitude)))
  ob2 <- subset(ob1, !((ob1$decimalLatitude == 0) & (ob1$decimalLongitude == 0)))

  if (nrow(ob2) > 0) {
    sp_ob <- ob2[c("scientificName", "decimalLatitude", "decimalLongitude", "depth", "eventDate")]
  } else {
    warning("No data were selected from OBIS")
    sp_ob <- data.frame()
  }

  # ---- Download from GBIF via rgbif ----
  df_gbif <- rgbif::occ_search(
    scientificName = key,
    geometry = sp_geometry,
    hasCoordinate = TRUE,
    eventDate = event_date,
    limit = 100000
  )
  gb0 <- as.data.frame(df_gbif$data)
  gb <- gb0[c("scientificName", "decimalLatitude", "decimalLongitude", "depth", "eventDate")]

  # Combine OBIS and GBIF data
  ob_gb <- rbind(gb, sp_ob)

  # Replace scientific names with the provided key
  ob_gb$scientificName <- key

  # ---- Format dates ----
  ob_gb$date <- NA
  ob_gb$year <- NA
  ob_gb$month <- NA
  ob_gb$day <- NA

  for (i in 1:nrow(ob_gb)) {
    ev <- as.character(ob_gb$eventDate[i])
    nch <- nchar(ev)

    if (nch > 10) {
      # Format like "2015-01-01T12:00:00Z"
      date0 <- strsplit(ev, split = "T") |> sapply(`[`, 1)
      inter <- strsplit(date0, split = "-") |> unlist() |> as.numeric()
      ob_gb$date[i] <- as.Date(date0, tryFormats = "%Y-%m-%d")
      ob_gb$year[i] <- inter[1]
      ob_gb$month[i] <- inter[2]
      ob_gb$day[i] <- inter[3]
    } else if (nch == 10) {
      # Format "2015-01-01"
      inter <- strsplit(ev, split = "-") |> unlist() |> as.numeric()
      ob_gb$date[i] <- as.Date(ev, tryFormats = "%Y-%m-%d")
      ob_gb$year[i] <- inter[1]
      ob_gb$month[i] <- inter[2]
      ob_gb$day[i] <- inter[3]
    } else if (nch == 4) {
      # Only year provided
      ob_gb$date[i] <- NA
      ob_gb$year[i] <- as.numeric(ev)
      ob_gb$month[i] <- NA
      ob_gb$day[i] <- NA
    } else {
      ob_gb$date[i] <- NA
      ob_gb$year[i] <- NA
      ob_gb$month[i] <- NA
      ob_gb$day[i] <- NA
    }
  }

  # Remove rows without a valid date
  ob_gb <- ob_gb[!is.na(ob_gb$date), ]

  # ---- Clean and rename columns ----
  depth <- stats::na.omit(ob_gb$depth)
  depth_summary <- summary(as.numeric(depth))

  ob_gb$decimalLatitude <- round(ob_gb$decimalLatitude, 5)
  ob_gb$decimalLongitude <- round(ob_gb$decimalLongitude, 5)

  ob_gb <- ob_gb |>
    dplyr::select(-eventDate)

  names(ob_gb) <- c("Species", "decimalLatitude", "decimalLongitude", "depth",
                    "date", "year", "month", "day")

  ob_gb <- unique(ob_gb)

  # ---- Plot occurrences on map ----
  pf <- obistools::plot_map(ob_gb, zoom = TRUE)

  # ---- Print summary information ----
  cat("Summary of depth:\n")
  print(depth_summary)
  cat("Original occurrences:\n")
  print(str(ob_gb))
  cat("Plot occurrences on worldmap:\n")
  print(pf)

  # ---- Return results ----
  return(list(
    depth_summary = depth_summary,
    sp_occ = ob_gb,
    plot = pf
  ))
}
