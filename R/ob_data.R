#' Download & clean OBIS and GBIF occurrence data
#'
#' @param key Character. The species scientific name (e.g., "Rhincodon typus").
#' @param sp_geometry Character. The geographic range of the study area as a WKT polygon
#'   (e.g., "POLYGON ((20 -41, 20 40, 180 40, 180 -41, 20 -41))").
#' @param start_date Character. Start date used in `robis::occurrence()`, format "YYYY-MM-DD".
#' @param end_date Character. End date used in `robis::occurrence()`, format "YYYY-MM-DD".
#' @param event_date Character. Year range used in `rgbif::occ_search()`, format "YYYY-MM-DD,YYYY-MM-DD".
#'
#' @return A list containing:
#'   \item{sp_occ}{Data frame of cleaned occurrence records with columns:
#'     `Species`, `decimalLatitude`, `decimalLongitude`, `date`, `year`, `month`, `day`.}
#'   \item{plot}{A ggplot map plot of occurrences.}
#'
#' @importFrom robis occurrence
#' @importFrom rgbif occ_search
#' @importFrom ggplot2 ggplot geom_polygon geom_point map_data coord_quickmap theme_minimal labs
#' @importFrom dplyr select bind_rows
#'
#' @examples
#' \donttest{
#' ws_occ <- ob_data(
#'   key = "Rhincodon typus",
#'   sp_geometry = "POLYGON ((90 -41, 90 40, 180 40, 180 -41, 90 -41))",
#'   start_date = "2010-01-01",
#'   end_date = "2025-12-30",
#'   event_date = "2010-01-01,2025-12-30"
#' )
#' }
#'
#' @export
ob_data <- function(key, sp_geometry, start_date, end_date, event_date) {

  # ---- Internal field standardization function ----
  options(warn = 1)

  map_fields <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)

    field_synonyms <- list(
      species = c("species","scientificName", "scientificname", "scientific_name", "name"),
      decimalLatitude = c("decimalLatitude", "decimallatitude", "lat", "latitude", "y", "decimal_latitude"),
      decimalLongitude = c("decimalLongitude", "decimallongitude", "lon", "long", "longitude", "x", "decimal_longitude"),
      eventDate = c("eventDate", "eventdate", "date", "samplingDate", "occurrenceDate")
    )

    for (target in names(field_synonyms)) {
      if (target %in% names(df) && !all(is.na(df[[target]]))) {
        message(sprintf("\nInfo: Standard column '%s' already present and used.", target))
        next
      }
      synonym_columns <- field_synonyms[[target]]
      found_col <- synonym_columns[synonym_columns %in% names(df)]
      if (length(found_col) > 0) {
        source_col <- found_col[1]
        df[[target]] <- df[[source_col]]
        message(sprintf("\nInfo: Standard column '%s' mapped to standard name '%s'", source_col, target))
      } else {
        df[[target]] <- NA
        warning(sprintf("\nWarning: No column found for required field '%s'. Created with NAs.", target))
      }
    }

    final_cols <- c("species", "decimalLatitude", "decimalLongitude", "eventDate")
    df <- df[, intersect(final_cols, names(df)), drop = FALSE]
    return(df)
  }

  # ---- Internal plotting function ----
  plot_map_alt <- function(data, zoom = TRUE,
                           point_size = 1, point_alpha = 0.6, point_color = "steelblue",
                           land_fill = "grey90", land_color = "grey70", ...) {
    world_df <- ggplot2::map_data("world")
    p <- ggplot2::ggplot() +
      ggplot2::geom_polygon(data = world_df,
                            ggplot2::aes(x = long, y = lat, group = group),
                            fill = land_fill, color = land_color, linewidth = 0.2)
    p <- p +
      ggplot2::geom_point(data = data,
                          ggplot2::aes(x = decimalLongitude, y = decimalLatitude),
                          size = point_size, alpha = point_alpha,
                          color = point_color, ...)
    if (zoom) {
      x_range <- range(data$decimalLongitude, na.rm = TRUE)
      y_range <- range(data$decimalLatitude, na.rm = TRUE)
      x_margin <- diff(x_range) * 0.05
      y_margin <- diff(y_range) * 0.05
      p <- p +
        ggplot2::coord_quickmap(xlim = c(x_range[1] - x_margin, x_range[2] + x_margin),
                                ylim = c(y_range[1] - y_margin, y_range[2] + y_margin))
    } else {
      p <- p + ggplot2::coord_quickmap()
    }
    p <- p +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Longitude", y = "Latitude")
    return(p)
  }

  # ---- Main Data Processing ----

  # Define empty template for no-data cases
  empty_df <- data.frame(
    species = character(),
    decimalLatitude = numeric(),
    decimalLongitude = numeric(),
    eventDate = character(),
    stringsAsFactors = FALSE
  )

  # ---- Download from OBIS ----
  df <- robis::occurrence(
    scientificname = key,
    geometry = sp_geometry,
    startdate = start_date,
    enddate = end_date
  )
  cat("\n")
  ob0 <- as.data.frame(df)

  sp_ob <- empty_df  # initialize empty

  if (nrow(ob0) > 0) {
    ob1 <- map_fields(ob0)
    ob1 <- subset(ob1, !(is.na(ob1$decimalLatitude) | is.na(ob1$decimalLongitude)))
    sp_ob <- subset(ob1, !((ob1$decimalLatitude == 0) & (ob1$decimalLongitude == 0)))
  } else {
    warning("No data were selected from OBIS")
  }

  # ---- Download from GBIF ----
  df_gbif <- rgbif::occ_search(
    scientificName = key,
    geometry = sp_geometry,
    hasCoordinate = TRUE,
    eventDate = event_date,
    limit = 100000
  )
  gb0 <- as.data.frame(df_gbif$data)

  sp_gb <- empty_df  # initialize empty

  if (nrow(gb0) > 0) {
    gb1 <- map_fields(gb0)
    gb1 <- subset(gb1, !(is.na(gb1$decimalLatitude) | is.na(gb1$decimalLongitude)))
    sp_gb <- subset(gb1, !((gb1$decimalLatitude == 0) & (gb1$decimalLongitude == 0)))
  } else {
    warning("No data were selected from GBIF")
  }

  # ---- Combine and Clean Data ----
  ob_gb <- dplyr::bind_rows(sp_gb, sp_ob)

  # Replace scientific names with the provided key
  ob_gb$scientificName <- key

  # ---- Format dates ----
  ob_gb$date <- NA
  ob_gb$year <- NA
  ob_gb$month <- NA
  ob_gb$day <- NA

  for (i in seq_len(nrow(ob_gb))) {
    ev <- as.character(ob_gb$eventDate[i])
    nch <- nchar(ev)
    if (nch > 10) {
      # Format like "2015-01-01T12:00:00Z"
      date0 <- strsplit(ev, split = "T")[[1]][1]
      inter <- suppressWarnings(as.numeric(unlist(strsplit(date0, split = "-"))))
      ob_gb$date[i] <- as.Date(date0, tryFormats = "%Y-%m-%d")
      ob_gb$year[i] <- inter[1]
      ob_gb$month[i] <- inter[2]
      ob_gb$day[i] <- inter[3]
    } else if (nch == 10) {
      # Format "2015-01-01"
      inter <- suppressWarnings(as.numeric(unlist(strsplit(ev, split = "-"))))
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

  # ---- Final cleaning ----
  ob_gb$decimalLatitude <- round(ob_gb$decimalLatitude, 5)
  ob_gb$decimalLongitude <- round(ob_gb$decimalLongitude, 5)

  # Select and rename final output columns
  final_cols <- c("species", "decimalLatitude", "decimalLongitude", "date", "year", "month", "day")
  ob_gb <- ob_gb[, intersect(final_cols, names(ob_gb))]
  ob_gb <- unique(ob_gb)

  # ---- Plot occurrences on map ----
  pf <- plot_map_alt(ob_gb, zoom = TRUE)

  # ---- Print summary information ----
  cat("Original occurrences:\n")
  print(str(ob_gb))
  cat("Plot occurrences on worldmap:\n")
  print(pf)

  # ---- Return results ----
  return(list(
    sp_occ = ob_gb,
    plot = pf
  ))
}
