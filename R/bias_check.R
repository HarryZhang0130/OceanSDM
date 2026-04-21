#' Check spatial and temporal sampling bias
#'
#' @param occ_data Data frame. Species occurrence records. Must contain columns:
#'   `Species`, `x`, `y`, `year`, `month`, `day`.
#' @param temporal_level Character. Temporal level at which temporal sampling
#'   bias is to be tested. One of `month`, `quarter`, or `season`.
#' @param spatial_method Character. Method used to generate species habitat area
#'   for spatial random sampling. One of `convex_hull` or `range_map`.
#' @param world_path Character. File path to a world map shapefile (required when
#'   `spatial_method = "convex_hull"`).
#' @param range_path Character. File path to a species range map shapefile
#'   (required when `spatial_method = "range_map"`). Must contain a column
#'   `sci_name`.
#' @param prj Character. Coordinate reference system to be used for spatial maps,
#'   given as a PROJ string (e.g., "+proj=longlat +datum=WGS84").
#'
#' @return A list of class `BiasCheckResult` containing the following components:
#'   \item{Temporal_level}{The temporal level used for testing.}
#'   \item{Spatial_method}{The spatial method used.}
#'   \item{Statistical_tests}{List of statistical test results for temporal
#'     bias and for spatial bias per temporal group.}
#'   \item{Plots}{List of `ggplot` objects showing temporal and spatial bias.}
#'   \item{Summary}{List with summary statistics (e.g., total records,
#'     number of temporal groups).}
#' @examples
#' \donttest{
#' colnames(corrected_ob)
#' occ_ob<-corrected_ob[,c(1:8)]
#' colnames(occ_ob)[2:3]<-c("y","x")
#' ## Add the quarter column
#' library(dplyr)
#' occ_ob2 <- occ_ob %>%
#'  mutate(
#'      quarter = case_when(
#'            month %in% 1:3   ~ "Q1",
#'               month %in% 4:6   ~ "Q2",
#'                 month %in% 7:9  ~ "Q3",
#'                  month %in% 10:12 ~ "Q4",
#'                   TRUE ~ NA_character_
#'                   ))%>% filter(!is.na(quarter))
#' range_path<-"F:/whaleshark_sdm/range_map/R_typus_wip2.shp"
#' bias_test <- bias_check(occ_data=occ_ob,
#'                         temporal_level="quarter",
#'                         spatial_method="range_map",
#'                         range_path=range_path,
#'                         prj = "+proj=longlat +datum=WGS84")
#' }
#' @export
bias_check <- function(occ_data,
                       temporal_level,
                       spatial_method,
                       world_path,
                       range_path,
                       prj = "+proj=longlat +datum=WGS84") {

  # Check required columns
  required_cols <- c("Species", "x", "y", "year", "month", "day")
  if (!all(required_cols %in% colnames(occ_data))) {
    missing_cols <- setdiff(required_cols, colnames(occ_data))
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Add quarter and season variables if not present
  if (!"quarter" %in% colnames(occ_data)) {
    occ_data <- occ_data |>
      dplyr::mutate(
        quarter = dplyr::case_when(
          month %in% 1:3 ~ "Q1",
          month %in% 4:6 ~ "Q2",
          month %in% 7:9 ~ "Q3",
          month %in% 10:12 ~ "Q4",
          TRUE ~ NA_character_
        )) |>
      dplyr::filter(!is.na(quarter))
  }

  if (!"season" %in% colnames(occ_data)) {
    occ_data <- occ_data |>
      dplyr::mutate(date = lubridate::make_date(year, month, day),
                    adjusted_date = date - lubridate::days(76),  # Adjust for sea temperature lag (73 days)
                    adjusted_month = lubridate::month(adjusted_date),
                    hemisphere = ifelse(y >= 0, "Northern", "Southern"),
                    season = dplyr::case_when(
                      hemisphere == "Northern" & adjusted_month %in% 3:5 ~ "Spring",
                      hemisphere == "Northern" & adjusted_month %in% 6:8 ~ "Summer",
                      hemisphere == "Northern" & adjusted_month %in% 9:11 ~ "Autumn",
                      hemisphere == "Northern" & adjusted_month %in% c(12,1,2) ~ "Winter",
                      hemisphere == "Southern" & adjusted_month %in% 3:5 ~ "Autumn",
                      hemisphere == "Southern" & adjusted_month %in% 6:8 ~ "Winter",
                      hemisphere == "Southern" & adjusted_month %in% 9:11 ~ "Spring",
                      hemisphere == "Southern" & adjusted_month %in% c(12,1,2) ~ "Summer",
                      TRUE ~ NA_character_
                    )) |>
      dplyr::filter(!is.na(season))
  }

  ## 2. Parameter validation
  temporal_level <- match.arg(arg = temporal_level,
                              choices = c("month", "quarter", "season"),
                              several.ok = FALSE)
  spatial_method <- match.arg(spatial_method,
                              choices = c("convex_hull", "range_map"),
                              several.ok = FALSE)

  # Validate file paths if provided
  if (spatial_method == "convex_hull" && !is.null(world_path) && !file.exists(world_path)) {
    stop(paste("world_path file does not exist:", world_path))
  }
  if (spatial_method == "range_map" && !is.null(range_path) && !file.exists(range_path)) {
    stop(paste("range_path file does not exist:", range_path))
  }

  ## 3. Create lists to store results
  res.list <- list()
  plot.list <- list()

  ## 4. Temporal bias test
  # Generate frequency table
  freq_table <- table(occ_data[[temporal_level]])
  # Expected frequencies under uniform distribution
  expected_freq <- rep(1 / length(freq_table), length(freq_table))
  # Chi-squared test
  chi_test <- stats::chisq.test(x = freq_table, p = expected_freq)
  res.list[["Temporal_bias"]] <- list(
    test_type = "Chi-squared test",
    statistic = chi_test$statistic,
    p_value = chi_test$p.value,
    df = chi_test$parameter,
    observed = as.numeric(freq_table),
    expected = expected_freq * sum(freq_table)
  )

  # Temporal bias plot
  time_plot <- ggplot2::ggplot(occ_data, ggplot2::aes(x = .data[[temporal_level]])) +
    ggplot2::geom_bar(fill = "steelblue", color = "black", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = mean(freq_table),
                        linetype = "dashed", color = "red", linewidth = 1) +
    ggplot2::labs(title = paste("Temporal bias -", temporal_level),
                  subtitle = paste("Chi-squared =", round(chi_test$statistic, 3),
                                   ", p =", format.pval(chi_test$p.value, digits = 3)),
                  x = temporal_level,
                  y = "Number of occurrences") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  plot.list[["Temporal_bias_plot"]] <- time_plot

  ## 5. Spatial bias test (grouped by temporal level)
  unique_temporal_values <- unique(occ_data[[temporal_level]])

  for (temp_val in unique_temporal_values) {
    # Subset occurrences for current time period
    occ_temp <- occ_data[occ_data[[temporal_level]] == temp_val, ]
    if (nrow(occ_temp) < 3) {
      warning(paste("Time period", temp_val,
                    "has fewer than 3 occurrence records. Skipping spatial bias test."))
      next
    }

    # Define spatial area
    if (spatial_method == "convex_hull") {
      # Create convex hull
      points_sf <- sf::st_as_sf(occ_temp, coords = c("x", "y"), crs = prj)
      if (nrow(points_sf) >= 3) {
        convex_hull <- sf::st_convex_hull(sf::st_union(points_sf))
      } else {
        stop(paste("Time period", temp_val, "has fewer than 3 valid points."))
      }

      # Check if we need to intersect with land (only if world map provided)
      if (!is.null(world_path) && file.exists(world_path)) {
        worldmap <- sf::st_read(world_path, quiet = TRUE)
        worldmap <- sf::st_transform(worldmap, sf::st_crs(convex_hull))
        # Check if convex hull intersects land
        intersection_test <- sf::st_intersects(convex_hull, worldmap, sparse = FALSE)
        if (any(intersection_test)) {
          # Convex hull intersects land, remove land part
          message(paste("Time period", temp_val,
                        ": convex hull intersects land. Removing land portion."))
          tryCatch({
            ocean_area <- sf::st_difference(convex_hull, sf::st_union(worldmap))
            # Check if resulting ocean area is valid
            if (sf::st_is_empty(ocean_area) || sf::st_area(ocean_area) == 0) {
              warning(paste("Time period", temp_val,
                            ": unable to obtain valid ocean area. Using original convex hull."))
              area <- convex_hull
            } else {
              area <- ocean_area
            }
          }, error = function(e) {
            warning(paste("Time period", temp_val,
                          ": error during land removal. Using original convex hull. Error:", e$message))
            area <- convex_hull
          })
        } else {
          # No intersection with land, use convex hull directly
          message(paste("Time period", temp_val,
                        ": convex hull does not intersect land. Using convex hull directly."))
          area <- convex_hull
        }
      } else {
        # No world map provided, use convex hull directly
        message(paste("Time period", temp_val,
                      ": no world map provided. Using convex hull directly."))
        area <- convex_hull
      }

      # Ensure area is valid
      area <- sf::st_make_valid(area)
      message(paste("Time period", temp_val, "area created successfully."))

    } else if (spatial_method == "range_map") {
      if (is.null(range_path)) {
        stop("When using 'range_map' method, 'range_path' parameter is required.")
      }
      range_map <- sf::st_read(range_path, quiet = TRUE)
      # Assume range map has a column 'sci_name'
      if (!"sci_name" %in% colnames(range_map)) {
        stop("The range_map must contain a column named 'sci_name'.")
      }
      species_name <- unique(occ_temp$Species)
      if (length(species_name) > 1) {
        warning("Multiple species names found. Using the first one.")
        species_name <- species_name[1]
      }
      species_range <- range_map |> dplyr::filter(sci_name == species_name)
      if (nrow(species_range) == 0) {
        warning(paste("Range not found for species", species_name,
                      ". Using convex hull instead."))
        points_sf <- sf::st_as_sf(occ_temp, coords = c("x", "y"), crs = prj)
        area <- sf::st_convex_hull(sf::st_union(points_sf))
      } else {
        # Create buffer if needed (not used here, just making valid)
        area <- sf::st_make_valid(species_range)
      }
    }

    # Ensure points are within the area
    points_sf <- sf::st_as_sf(occ_temp, coords = c("x", "y"), crs = prj)
    points_in_area <- sf::st_intersection(points_sf, area)
    if (nrow(points_in_area) < nrow(occ_temp)) {
      warning(paste("Time period", temp_val, ": some points fall outside the defined area."))
    }

    # Compute nearest neighbor distances
    if (nrow(points_in_area) >= 2) {
      coords_matrix <- sf::st_coordinates(points_in_area)
      dist_matrix <- as.matrix(stats::dist(coords_matrix))
      diag(dist_matrix) <- NA
      min_nndist_actual <- apply(dist_matrix, 1, min, na.rm = TRUE)

      # Generate random points
      area_valid <- sf::st_make_valid(area)
      n_random_points <- nrow(points_in_area)
      # Ensure enough space for random points
      if (as.numeric(sf::st_area(area_valid)) > 0) {
        random_points <- sf::st_sample(area_valid, size = n_random_points, type = "random")
        random_coords <- sf::st_coordinates(random_points)
        # Compute nearest neighbor distances for random points
        random_dist_matrix <- as.matrix(stats::dist(random_coords))
        diag(random_dist_matrix) <- NA
        min_nndist_random <- apply(random_dist_matrix, 1, min, na.rm = TRUE)

        # t-test comparison
        if (length(min_nndist_actual) > 1 && length(min_nndist_random) > 1) {
          t_test <- stats::t.test(min_nndist_actual, min_nndist_random)
          res.list[[paste0("Spatial_bias_", temporal_level, "_", temp_val)]] <- list(
            test_type = "Two-sample t-test",
            statistic = t_test$statistic,
            p_value = t_test$p.value,
            df = t_test$parameter,
            mean_actual = mean(min_nndist_actual),
            mean_random = mean(min_nndist_random),
            sd_actual = stats::sd(min_nndist_actual),
            sd_random = stats::sd(min_nndist_random),
            n_actual = length(min_nndist_actual),
            n_random = length(min_nndist_random)
          )

          # Create spatial plot
          range_bbox <- sf::st_bbox(area_valid)
          occ_bbox <- sf::st_bbox(points_in_area)
          plot_bbox <- c(
            xmin = min(occ_bbox["xmin"], range_bbox["xmin"], na.rm = TRUE),
            ymin = min(occ_bbox["ymin"], range_bbox["ymin"], na.rm = TRUE),
            xmax = max(occ_bbox["xmax"], range_bbox["xmax"], na.rm = TRUE),
            ymax = max(occ_bbox["ymax"], range_bbox["ymax"], na.rm = TRUE)
          )

          spatial_plot <- ggplot2::ggplot() +
            ggplot2::geom_sf(data = area, fill = "lightgray", color = "darkgray", alpha = 0.5) +
            ggplot2::geom_sf(data = points_in_area, color = "black", size = 2, alpha = 0.7) +
            ggplot2::geom_point(data = as.data.frame(random_coords),
                                ggplot2::aes(x = X, y = Y), color = "red", size = 1, alpha = 0.3) +
            ggplot2::labs(title = paste("Spatial bias -", temporal_level, ":", temp_val),
                          subtitle = paste("t =", round(t_test$statistic, 3),
                                           ", p =", format.pval(t_test$p.value, digits = 3)),
                          x = "Longitude",
                          y = "Latitude") +
            ggplot2::theme_minimal() +
            ggplot2::coord_sf(
              xlim = c(plot_bbox["xmin"], plot_bbox["xmax"]),
              ylim = c(plot_bbox["ymin"], plot_bbox["ymax"]),
              expand = FALSE
            )

          plot.list[[paste0("Spatial_bias_plot_", temporal_level, "_", temp_val)]] <- spatial_plot

        } else {
          warning(paste("Time period", temp_val, ": insufficient points for t-test."))
        }
      } else {
        warning(paste("Time period", temp_val, ": valid area is too small."))
      }
    } else {
      warning(paste("Time period", temp_val, ": fewer than 2 valid points in area."))
    }
  }

  ## 6. Return results
  result <- list(
    Temporal_level = temporal_level,
    Spatial_method = spatial_method,
    Statistical_tests = res.list,
    Plots = plot.list,
    Summary = list(
      n_total_records = nrow(occ_data),
      n_temporal_groups = length(unique(occ_data[[temporal_level]])),
      temporal_groups = unique_temporal_values
    )
  )
  class(result) <- "BiasCheckResult"
  return(result)
}
