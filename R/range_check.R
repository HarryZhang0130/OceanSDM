#' Check if occurrences of multiple species fall within a defined range map with a buffer
#' and remove species with limited data
#'
#' @param occ_path Character. Path to the occurrence CSV file. Must contain columns:
#'   `species`, `decimalLongitude`, `decimalLatitude` (or `x`, `y`). Adjust extraction accordingly.
#' @param range_path Character. Path to the species range map (shapefile) containing
#'   a column `sci_name` with species scientific names.
#' @param buffer_km Numeric. Buffer distance (in kilometres) around the range map
#'   to consider occurrences as valid.
#'
#' @return A list containing:
#'   \item{successful.sp.checklist}{Character vector of species that passed the
#'     range check.}
#'   \item{valid.occ}{Data frame of valid occurrence points (those falling within
#'     the buffered range).}
#'   \item{problematic.sp.checklist}{Character vector of species that either
#'     lacked a range map, had processing errors, or had insufficient valid points.}
#'   \item{problematic.occ}{Data frame of occurrences belonging to problematic species.}
#'
#' @details This function requires a helper function `range_check_single` (not provided here)
#'   to process individual species. The user must define `range_check_single` with
#'   arguments `(sci_name, species_data, range_path, buffer_km)` returning
#'   a list with elements `valid_points`, `invalid_points`, and `error` (or `NULL`).
#'
#'
#' @examples
#' \donttest{
#' occ_path <-system.file("extdata", "R_typus.csv", package = "OceanSDM")
#' shp_path <- system.file("extdata", "R_typus_wip2.shp", package = "OceanSDM")
#' result <- range_check(
#'   occ_path = occ_path,
#'   range_path = shp_path
#' )
#' }
#' @export
range_check <- function(occ_path, range_path, buffer_km) {
  # ---- Validate inputs ----
  if (!file.exists(occ_path)) stop("occ_path file does not exist: ", occ_path)
  if (!file.exists(range_path)) stop("range_path file does not exist: ", range_path)

  # ---- Read occurrence data ----
  cat("Start to check which species do not have range maps....\n")
  occ_data <- utils::read.csv(occ_path)

  # Check required columns (adjust column names as needed)
  required_cols <- c("species", "decimalLongitude", "decimalLatitude")
  if (!all(required_cols %in% colnames(occ_data))) {
    # Try alternative names
    if ("x" %in% colnames(occ_data) && "y" %in% colnames(occ_data)) {
      occ_data$decimalLongitude <- occ_data$x
      occ_data$decimalLatitude <- occ_data$y
    } else {
      stop("Occurrence data must contain columns 'species', and either ('decimalLongitude','decimalLatitude') or ('x','y')")
    }
  }

  cat("Species checklist:", paste(unique(occ_data$species), collapse = ", "), "\n")

  # ---- Read range map ----
  range.map <- sf::st_read(range_path, quiet = TRUE)
  cat("Species that have range map:", paste(unique(range.map$sci_name), collapse = ", "), "\n")

  # ---- Identify species with/without range maps ----
  species_in_data <- unique(occ_data$species)
  species_in_range_map <- unique(range.map$sci_name)

  species_with_maps <- intersect(species_in_data, species_in_range_map)
  species_without_maps <- setdiff(species_in_data, species_in_range_map)

  cat("Species without range maps:\n", paste(species_without_maps, collapse = ", "), "\n")
  cat("No. of species without range maps:\n", length(species_without_maps), "\n")

  # Subset occurrence data for species that have maps
  df2 <- occ_data[occ_data$species %in% species_with_maps, ]
  df3 <- occ_data[occ_data$species %in% species_without_maps, ]

  # ---- Process each species with a range map ----
  cat("Start processing occurrence points for all species...\n")

  species_groups <- split(df2, df2$species)
  species_names <- names(species_groups)

  pb <- utils::txtProgressBar(min = 0, max = length(species_names), style = 3)

  # Assuming range_check_single is defined elsewhere
  results <- purrr::map(seq_along(species_names), function(i) {
    sci.name <- species_names[i]
    utils::setTxtProgressBar(pb, i)
    cat("\nProcessing species (", i, "/", length(species_names), "): ", sci.name, "\n")

    # Call the helper function (must be provided by user)
    result <- range_check_single(sci.name, species_groups[[sci.name]],
                                 range_path, buffer_km)
    result$species_name <- sci.name
    return(result)
  })

  close(pb)
  names(results) <- species_names

  # ---- Extract valid and invalid points across all species ----
  all_valid_points <- purrr::map_dfr(results, function(x) {
    if (!is.null(x$error)) {
      df <- x$valid_points
      if (nrow(df) > 0) {
        df$processing_error <- x$error
        df$species <- x$species_name
        return(df)
      } else {
        return(data.frame(Species = x$species_name, processing_error = x$error))
      }
    } else {
      df <- x$valid_points
      if (nrow(df) > 0) {
        df$species <- x$species_name
        return(df)
      } else {
        return(data.frame(Species = x$species_name))
      }
    }
  })

  invalid_points <- purrr::map_dfr(results, function(x) {
    if (is.null(x$error) && nrow(x$invalid_points) > 0) {
      df <- x$invalid_points
      df$species <- x$species_name
      return(df)
    } else {
      return(data.frame())
    }
  })

  # ---- Collect errors ----
  errors <- purrr::map(results, ~ .x$error) |>
    purrr::compact() |>
    tibble::enframe(name = "Species", value = "error") |>
    tidyr::unnest(error)

  # ---- Report summary ----
  if (nrow(errors) > 0) {
    cat("\nWarning: Errors occurred while processing the following species:\n")
    print(errors)
  } else {
    cat("\nAll species processed successfully, no errors.\n")
  }

  valid_points <- all_valid_points  # rename for clarity

  cat("Processing complete!\n")
  cat("Original total points:", nrow(df2), "\n")
  cat("Final valid points:", nrow(valid_points), "\n")
  cat("Final invalid points:", nrow(invalid_points), "\n")
  cat("Overall retention rate:", round(nrow(valid_points) / nrow(df2) * 100, 2), "%\n")

  # ---- Identify problematic species (insufficient valid points) ----
  all_species <- unique(df2$species)
  valid_counts <- table(valid_points$species)

  species_stats <- data.frame(
    Species = all_species,
    valid_points_count = 0,
    stringsAsFactors = FALSE
  )
  for (sp in names(valid_counts)) {
    species_stats$valid_points_count[species_stats$species == sp] <- valid_counts[sp]
  }

  error_species <- if (nrow(errors) > 0) errors$Species else character(0)
  species_stats$had_error <- FALSE
  if (length(error_species) > 0) {
    species_stats$had_error <- species_stats$Species %in% error_species
  }

  problematic_species <- species_stats[species_stats$had_error, ]

  cat("Number of problematic species:", nrow(problematic_species), "\n")
  cat("Species with processing errors:", sum(problematic_species$had_error), "\n")

  # ---- Extract problematic species data ----
  if (nrow(problematic_species) > 0) {
    problematic_species_data <- df2[df2$species %in% problematic_species$species, ]
    problematic_species_data <- merge(problematic_species_data, problematic_species, by = "Species")
    problematic_species_data$problem_type <- ifelse(
      problematic_species_data$had_error,
      "processing_error",
      "insufficient_points"
    )
    problematic_summary <- aggregate(
      list(total_points = problematic_species_data$species),
      by = list(
        Species = problematic_species_data$species,
        problem_type = problematic_species_data$problem_type,
        valid_points_count = problematic_species_data$valid_points_count
      ),
      FUN = length
    )

    cat("\nProblematic species details:\n")
    for (i in 1:nrow(problematic_summary)) {
      sp <- problematic_summary$species[i]
      ptype <- problematic_summary$problem_type[i]
      vpts <- problematic_summary$valid_points_count[i]
      tpts <- problematic_summary$total_points[i]
      if (ptype == "processing_error") {
        cat(sprintf("  %s: processing error (original points: %d)\n", sp, tpts))
      } else {
        cat(sprintf("  %s: insufficient valid points (valid/total: %d/%d)\n", sp, vpts, tpts))
      }
    }
  } else {
    cat("No problematic species found.\n")
    problematic_species_data <- data.frame()
    problematic_summary <- data.frame()
  }

  df_dd <- problematic_species_data
  list_check <- unique(problematic_summary$species)
  df_dc <- df2[df2$species %in% list_check, ]
  df4 <- rbind(df3, df_dc)
  list_e <- unique(df4$species)
  if(nrow(problematic_species)>0){
    cat("Problematic species identified:\n")
    print(list_e)
  }
  problematic_points <- occ_data[occ_data$species %in% list_e, ]
  # ---- Frequency table of valid points ----
  freq_table <- table(valid_points$species)
  print(freq_table)
  report_table <- as.data.frame(freq_table)
  colnames(report_table) <- c("Species", "Count")
  print(report_table)

  # ---- Return results ----
  list(
    successful.sp.checklist = list_check,
    valid.occ = valid_points,
    problematic.sp.checklist = list_e,
    problematic.occ = problematic_points
  )
}
