#' Calculate depth range within specified temperature range for each grid cell
#' using adaptive interpolation based on temp_bin results
#'
#' @param temp_bin_path Character. Root path of temp_bin generated files.
#' @param PTR Numeric vector of length 2. Temperature range, e.g., `c(23.5, 30.2)`.
#' @param time_type Character. Time type, should match that used in temp_bin:
#'   `"month"`, `"quarter"`, or `"season"`.
#' @param save_path Character. Path to save output files.
#' @param na_value Numeric. Value to replace NA (default: -9999).
#' @param use_makima Logical. Whether to use makima interpolation (requires
#'   `boostmath` package). If `FALSE` or package not available, falls back to spline.
#' @param N_cells Integer. Number of random cells to plot (5–30). If `0`, no plot is generated.
#' @param layer_option Character. Which depth interval to extract: `"top"` (shallowest)
#'   or `"bottom"` (deepest). Default is `"top"`.
#'
#' @return Invisibly returns a character vector of all generated file paths.
#' @export
#'
#' @examples
#' \donttest{
#' # Assuming temp_bin output exists in "path/to/temp_bin"
#' TB_depth(temp_bin_path = "path/to/temp_bin",
#'          PTR = c(23.5, 30.2), # thermal ranges, unit ℃
#'          time_type = "quarter",
#'          save_path = "./output",
#'          N_cells = 20)
#' }
TB_depth <- function(temp_bin_path, PTR, time_type = c("month", "quarter", "season"),
                     layer_option = "top", save_path, na_value = -9999,
                     use_makima = TRUE, N_cells = 100) {

  # ---- Argument validation ----
  time_type <- match.arg(time_type)
  if (N_cells < 0 || N_cells > 100) stop("N_cells must be between 0 and 100")
  if (N_cells > 0 && N_cells < 5) {
    warning("N_cells is less than 5, setting to 5")
    N_cells <- 5
  }
  if (length(PTR) != 2 || PTR[1] >= PTR[2]) {
    stop("PTR must be a vector of length 2 in increasing order, e.g., c(23.5, 30.2)")
  }

  # ---- Check optional package for makima ----
  makima_available <- FALSE
  if (use_makima && requireNamespace("boostmath", quietly = TRUE)) {
    makima_available <- TRUE
    cat("Using boostmath::makima for interpolation (when data points >= 4)\n")
  } else if (use_makima) {
    cat("Warning: boostmath package not available, using base R spline interpolation\n")
  }

  # ---- Create save directory ----
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

  # ---- Helper: find continuous intervals where temperature lies within PTR ----
  find_continuous_intervals <- function(depths, temps, PTR) {
    if (length(depths) < 2) {
      if (length(depths) == 1 && temps[1] >= PTR[1] && temps[1] <= PTR[2]) {
        return(list(list(min = depths[1], max = depths[1])))
      } else {
        return(list())
      }
    }
    # Ensure depths are sorted
    ord <- order(depths)
    depths <- depths[ord]
    temps <- temps[ord]

    intersect_depth <- function(d1, d2, t1, t2, thresh) {
      if (t1 == t2) return(NULL)
      d1 + (thresh - t1) * (d2 - d1) / (t2 - t1)
    }

    points <- numeric(0)
    n <- length(depths)
    for (i in 1:(n-1)) {
      d1 <- depths[i]; d2 <- depths[i+1]
      t1 <- temps[i]; t2 <- temps[i+1]
      if ((t1 - PTR[1]) * (t2 - PTR[1]) <= 0 && t1 != t2) {
        cross <- intersect_depth(d1, d2, t1, t2, PTR[1])
        if (!is.null(cross)) points <- c(points, cross)
      }
      if ((t1 - PTR[2]) * (t2 - PTR[2]) <= 0 && t1 != t2) {
        cross <- intersect_depth(d1, d2, t1, t2, PTR[2])
        if (!is.null(cross)) points <- c(points, cross)
      }
      if (t1 >= PTR[1] && t1 <= PTR[2]) points <- c(points, d1)
      if (t2 >= PTR[1] && t2 <= PTR[2]) points <- c(points, d2)
    }

    points <- unique(points)
    if (length(points) == 0) return(list())
    points <- sort(points)
    # Merge very close points
    merged_points <- points[1]
    for (p in points[-1]) {
      if (p - tail(merged_points, 1) > 0.001) merged_points <- c(merged_points, p)
    }
    points <- merged_points

    intervals <- list()
    start <- points[1]
    for (i in 2:length(points)) {
      d_mid <- (points[i-1] + points[i]) / 2
      seg <- which(depths[-1] >= d_mid)[1]
      if (is.na(seg)) seg <- n-1
      d_left <- depths[seg]; d_right <- depths[seg+1]
      t_left <- temps[seg]; t_right <- temps[seg+1]
      if (d_right != d_left) {
        t_mid <- t_left + (t_right - t_left) * (d_mid - d_left) / (d_right - d_left)
      } else {
        t_mid <- (t_left + t_right) / 2
      }
      if (t_mid >= PTR[1] && t_mid <= PTR[2]) {
        if (i == length(points)) {
          intervals <- c(intervals, list(list(min = start, max = points[i])))
        }
      } else {
        intervals <- c(intervals, list(list(min = start, max = points[i-1])))
        start <- points[i]
        if (i == length(points)) {
          intervals <- c(intervals, list(list(min = start, max = points[i])))
        }
      }
    }
    # Merge adjacent intervals
    if (length(intervals) > 1) {
      merged <- list()
      current <- intervals[[1]]
      for (i in 2:length(intervals)) {
        if (intervals[[i]]$min - current$max <= 0.001) {
          current$max <- intervals[[i]]$max
        } else {
          merged <- c(merged, list(current))
          current <- intervals[[i]]
        }
      }
      merged <- c(merged, list(current))
      intervals <- merged
    }
    return(intervals)
  }

  # ---- Determine time group folders ----
  if (time_type == "month") {
    folders <- paste0("M", sprintf("%02d", 1:12))
    group_names <- paste0("M", sprintf("%02d", 1:12))
  } else if (time_type == "quarter") {
    folders <- paste0("Q", 1:4)
    group_names <- paste0("Q", 1:4)
  } else { # season
    folders <- c("DJF", "MAM", "JJA", "SON")
    group_names <- folders
  }

  time_dirs <- list.dirs(temp_bin_path, recursive = FALSE, full.names = TRUE)
  valid_dirs <- c()
  valid_group_names <- c()
  for (f in seq_along(folders)) {
    folder_name <- folders[f]
    dir_path <- file.path(temp_bin_path, folder_name)
    if (dir.exists(dir_path)) {
      valid_dirs <- c(valid_dirs, dir_path)
      valid_group_names <- c(valid_group_names, group_names[f])
    }
  }
  if (length(valid_dirs) == 0) stop("No time group folders found in specified path")
  cat("Found", length(valid_dirs), "time groups:", paste(valid_group_names, collapse = ", "), "\n")

  # ---- Get depth layer information ----
  first_dir <- valid_dirs[1]
  first_files <- list.files(first_dir, pattern = "temp_.*\\.asc$", full.names = TRUE)
  if (length(first_files) == 0) stop("No temperature raster files found")
  depth_values <- as.numeric(stringr::str_extract(basename(first_files), "\\d+\\.?\\d*(?=\\.asc)"))
  depth_order <- order(depth_values)
  depth_values <- depth_values[depth_order]
  cat("Found", length(depth_values), "depth layers:", paste(depth_values, collapse = ", "), "\n")

  # ---- Plotting parameters ----
  method_linetypes <- c(
    "single_point" = NA,
    "linear" = "dashed",
    "makima" = "dotdash",
    "spline" = "dotted",
    "linear_fallback" = "longdash"
  )
  method_colors <- c(
    "linear" = "#E69F00",
    "makima" = "#56B4E9",
    "spline" = "#009E73",
    "linear_fallback" = "#CC79A7",
    "single_point" = "#999999"
  )
  method_linewidths <- c(
    "linear" = 2,
    "makima" = 2.5,
    "spline" = 2,
    "linear_fallback" = 1.5,
    "single_point" = 1
  )
  method_display_names <- c(
    "single_point" = "Single Point (no fit)",
    "linear" = "Linear (2-3 pts)",
    "makima" = "Makima (≥ 4 pts)",
    "spline" = "Spline (≥ 4 pts)",
    "linear_fallback" = "Linear Fallback"
  )

  # ---- Process each time group ----
  result_files <- list()

  for (d_idx in seq_along(valid_dirs)) {
    time_dir <- valid_dirs[d_idx]
    group_name <- valid_group_names[d_idx]
    cat(sprintf("\nProcessing time group %s...\n", group_name))

    group_save_path <- file.path(save_path, group_name)
    if (!dir.exists(group_save_path)) dir.create(group_save_path, recursive = TRUE)

    # Read all depth layer rasters for this group
    layer_files <- list.files(time_dir, pattern = "temp_.*\\.asc$", full.names = TRUE)
    if (length(layer_files) == 0) {
      cat(sprintf(" Group %s has no valid data, skipping\n", group_name))
      next
    }
    file_depths <- as.numeric(stringr::str_extract(basename(layer_files), "\\d+\\.?\\d*(?=\\.asc)"))
    depth_order <- order(file_depths)
    layer_files <- layer_files[depth_order]
    file_depths <- file_depths[depth_order]

    layer_rasters <- list()
    for (i in seq_along(layer_files)) {
      r <- raster::raster(layer_files[i])
      r[r == na_value] <- NA
      layer_rasters[[as.character(file_depths[i])]] <- r
    }
    template_raster <- layer_rasters[[1]]

    # Initialize result rasters
    min_depth_raster <- max_depth_raster <- range_depth_raster <- template_raster
    raster::values(min_depth_raster) <- NA
    raster::values(max_depth_raster) <- NA
    raster::values(range_depth_raster) <- NA

    n_cells <- raster::ncell(template_raster)
    cat(sprintf(" Interpolating and calculating depth range for %d cells...\n", n_cells))

    plot_data <- list()
    pb <- utils::txtProgressBar(min = 0, max = n_cells, style = 3)

    for (cell in 1:n_cells) {
      utils::setTxtProgressBar(pb, cell)

      temp_profile <- sapply(layer_rasters, function(r) r[cell])
      valid_idx <- which(!is.na(temp_profile))
      n_valid <- length(valid_idx)
      if (n_valid == 0) next

      valid_depths <- file_depths[valid_idx]
      valid_temps <- temp_profile[valid_idx]

      if (max(valid_temps) < PTR[1] || min(valid_temps) > PTR[2]) next

      interp_method <- "unknown"

      # Case 1: Single point
      if (n_valid == 1) {
        if (valid_temps[1] >= PTR[1] && valid_temps[1] <= PTR[2]) {
          min_depth_raster[cell] <- max_depth_raster[cell] <- range_depth_raster[cell] <- valid_depths[1]
          interp_method <- "single_point"
        }
        if (N_cells > 0) {
          plot_data[[length(plot_data) + 1]] <- list(
            cell = cell, depths = valid_depths, temps = valid_temps,
            method = interp_method, interp_depths = NULL, interp_temps = NULL
          )
        }
        next
      }

      # Case 2: 2-3 points -> linear interpolation
      if (n_valid <= 3) {
        depth_min <- min(valid_depths)
        depth_max <- max(valid_depths)
        fine_depths <- seq(depth_min, depth_max, length.out = 100)
        interp_func <- stats::approxfun(valid_depths, valid_temps, method = "linear", rule = 1)
        pred_temps <- interp_func(fine_depths)
        interp_method <- "linear"

        intervals <- find_continuous_intervals(fine_depths, pred_temps, PTR)
        if (length(intervals) > 0) {
          if (layer_option == "top") {
            selected <- intervals[[which.min(sapply(intervals, function(x) x$min))]]
          } else {
            selected <- intervals[[which.max(sapply(intervals, function(x) x$max))]]
          }
          min_depth_raster[cell] <- selected$min
          max_depth_raster[cell] <- selected$max
          range_depth_raster[cell] <- selected$max - selected$min
        }

        if (N_cells > 0) {
          plot_data[[length(plot_data) + 1]] <- list(
            cell = cell, depths = valid_depths, temps = valid_temps,
            method = interp_method, interp_depths = fine_depths, interp_temps = pred_temps
          )
        }
        next
      }

      # Case 3: >=4 points -> try makima, fallback to spline, then linear
      depth_min <- min(valid_depths)
      depth_max <- max(valid_depths)
      fine_depths <- seq(depth_min, depth_max, length.out = 200)
      pred_temps <- NULL

      if (makima_available) {
        tryCatch({
          interpolator <- boostmath::makima(valid_depths, valid_temps)
          pred_temps <- sapply(fine_depths, function(xi) interpolator$interpolate(xi))
          interp_method <- "makima"
        }, error = function(e) {
          pred_temps <- tryCatch({
            stats::spline(valid_depths, valid_temps, xout = fine_depths, method = "fmm")$y
          }, error = function(e2) {
            stats::approx(valid_depths, valid_temps, xout = fine_depths, rule = 1)$y
          })
          interp_method <<- if (length(pred_temps) == length(fine_depths)) "spline" else "linear_fallback"
        })
      } else {
        pred_temps <- tryCatch({
          stats::spline(valid_depths, valid_temps, xout = fine_depths, method = "fmm")$y
        }, error = function(e) {
          stats::approx(valid_depths, valid_temps, xout = fine_depths, rule = 1)$y
        })
        interp_method <- if (length(pred_temps) == length(fine_depths)) "spline" else "linear_fallback"
      }

      intervals <- find_continuous_intervals(fine_depths, pred_temps, PTR)
      if (length(intervals) > 0) {
        if (layer_option == "top") {
          selected <- intervals[[which.min(sapply(intervals, function(x) x$min))]]
        } else {
          selected <- intervals[[which.max(sapply(intervals, function(x) x$max))]]
        }
        min_depth_raster[cell] <- selected$min
        max_depth_raster[cell] <- selected$max
        range_depth_raster[cell] <- selected$max - selected$min
      }

      if (N_cells > 0) {
        plot_data[[length(plot_data) + 1]] <- list(
          cell = cell, depths = valid_depths, temps = valid_temps,
          method = interp_method, interp_depths = fine_depths, interp_temps = pred_temps
        )
      }
    }

    close(pb)

    n_valid_cells <- length(which(!is.na(raster::values(min_depth_raster))))
    cat(sprintf(" Group %s: %d cells have valid depth ranges\n", group_name, n_valid_cells))

    # Save rasters
    min_depth_raster[is.na(min_depth_raster)] <- na_value
    max_depth_raster[is.na(max_depth_raster)] <- na_value
    range_depth_raster[is.na(range_depth_raster)] <- na_value

    min_file <- file.path(group_save_path, paste0(group_name, "_min_depth.asc"))
    max_file <- file.path(group_save_path, paste0(group_name, "_max_depth.asc"))
    range_file <- file.path(group_save_path, paste0(group_name, "_depth_range.asc"))

    raster::writeRaster(min_depth_raster, min_file, format = "ascii", overwrite = TRUE, NAflag = na_value)
    raster::writeRaster(max_depth_raster, max_file, format = "ascii", overwrite = TRUE, NAflag = na_value)
    raster::writeRaster(range_depth_raster, range_file, format = "ascii", overwrite = TRUE, NAflag = na_value)

    result_files <- c(result_files, min_file, max_file, range_file)
    cat(sprintf(" Group %s processing complete, 3 files saved\n", group_name))

    # ---- Plotting (if requested) ----
    if (N_cells > 0 && length(plot_data) > 0) {
      tryCatch({
        n_available <- length(plot_data)
        n_sample <- min(N_cells, n_available)
        sampled_idx <- sample(1:n_available, n_sample)

        all_depths <- unlist(lapply(plot_data[sampled_idx], function(x) x$depths))
        depth_range <- range(all_depths, na.rm = TRUE)

        plot_file <- file.path(save_path, paste0("Temperature_profile_", group_name, ".tif"))
        grDevices::tiff(plot_file, width = 14, height = 8, units = "in", res = 300, compression = "lzw")

        graphics::layout(matrix(c(1,2), nrow = 1, ncol = 2), widths = c(3, 1))
        graphics::par(mar = c(5, 5, 4, 1), xpd = FALSE)

        graphics::plot(NULL,
                       xlim = range(unlist(lapply(plot_data[sampled_idx], function(x) x$temps)), na.rm = TRUE),
                       ylim = rev(range(all_depths)),
                       xlab = "Temperature (°C)", ylab = "Depth (m)",
                       main = paste(group_name, "- Temperature Profiles (Random Sample of", n_sample, "Cells)"),
                       cex.lab = 1.2, cex.main = 1.2)
        graphics::grid(col = "gray90", lty = "dotted")

        depth_colors <- grDevices::colorRampPalette(c("lightblue", "darkblue"))(100)
        methods_used <- c()

        for (i in sampled_idx) {
          data <- plot_data[[i]]
          if (!(data$method %in% methods_used)) methods_used <- c(methods_used, data$method)

          depths <- data$depths
          temps <- data$temps
          depth_norm <- (depths - depth_range[1]) / (depth_range[2] - depth_range[1] + 0.001)
          point_colors <- depth_colors[round(depth_norm * 99) + 1]
          point_sizes <- 0.5 + 2 * depth_norm
          graphics::points(temps, depths, pch = 16, col = point_colors, cex = point_sizes, lwd = 1.5)

          if (!is.null(data$interp_depths) && !is.null(data$interp_temps)) {
            graphics::lines(data$interp_temps, data$interp_depths,
                            col = method_colors[data$method], lty = method_linetypes[data$method],
                            lwd = method_linewidths[data$method])
          }
        }
        graphics::abline(v = PTR[1], col = "red", lty = "dashed", lwd = 1)
        graphics::abline(v = PTR[2], col = "red", lty = "dashed", lwd = 1)

        # Legend
        graphics::par(mar = c(5, 1, 4, 2), xpd = NA)
        graphics::plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "", main = "")

        used_methods <- unique(methods_used)
        used_methods <- used_methods[!is.na(used_methods) & used_methods != "unknown"]
        y_top <- 0.9
        y_spacing <- 0.07

        if (length(used_methods) > 0) {
          line_methods <- used_methods[used_methods != "single_point"]
          point_methods <- used_methods[used_methods == "single_point"]
          if (length(line_methods) > 0) {
            graphics::legend(0, y_top, title = "Interpolation Method",
                             legend = method_display_names[line_methods],
                             lty = method_linetypes[line_methods],
                             col = method_colors[line_methods],
                             lwd = method_linewidths[line_methods],
                             bty = "n", cex = 0.9, xpd = NA)
            y_top <- y_top - length(line_methods) * y_spacing - 0.1
          }
          if (length(point_methods) > 0) {
            graphics::legend(0, y_top, legend = method_display_names["single_point"],
                             pch = 16, col = method_colors["single_point"], pt.cex = 1.2,
                             bty = "n", cex = 0.9, xpd = NA)
            y_top <- y_top - 2 * y_spacing
          }
        }

        legend_depths <- seq(depth_range[1], depth_range[2], length.out = 5)
        legend_norm <- (legend_depths - depth_range[1]) / (depth_range[2] - depth_range[1] + 0.001)
        legend_colors <- depth_colors[round(legend_norm * 99) + 1]
        legend_sizes <- 0.5 + 2 * legend_norm
        graphics::legend(0, y_top, title = "Depth (m)", legend = round(legend_depths, 1),
                         pch = 16, col = legend_colors, pt.cex = legend_sizes, bty = "n", cex = 0.9, xpd = NA)

        grDevices::dev.off()
        cat(sprintf(" Saved temperature profile plot: %s\n", basename(plot_file)))
      }, error = function(e) {
        warning(sprintf("Plotting failed: %s", e$message))
      })
    }
  }

  cat("\nTB_depth processing complete! Generated", length(result_files), "files\n")
  invisible(result_files)
}
