#' Calculate depth range within specified temperature range for each grid cell
#' using adaptive interpolation based on temp_bin results
#'
#' @param temp_bin_path Character. Root path of temp_bin generated files.
#' @param PTR Numeric vector of length 2. Temperature range, e.g., `c(23.5, 30.2)`.
#' @param time_type Character. Time type, should match that used in temp_bin:
#'   `"month"`, `"quarter"`, or `"season"`.
#' @param layer_option Character. Which depth interval to extract: `"top"` (shallowest)
#'   or `"bottom"` (deepest). Default is `"top"`.
#' @param save_path Character. Path to save output files.
#' @param na_value Numeric. Value to replace NA (default: -9999).
#' @param use_makima Logical. Whether to use makima interpolation. If `FALSE` , falls back to spline.
#' @param N_cells Integer. Number of random cells to plot (5--30). If `0`, no plot is generated.
#' @param ncores Integer. Number of cores for parallel processing (uses `pbapply`). Default: 1.
#'
#' @return Invisibly returns a character vector of all generated file paths.
#' @examples
#' \donttest{
#' df<-read.csv(system.file("extdata","summary_niche_fit.csv",package = "OceanSDM"))
#' PTR<-df[2:3,2] # preferred temperature range
#' TTR<-df[4:5,2] # tolerable temperature range
#' temp_bin_path<-"F:/whaleshark_sdm/4D" # the same path as the save_path of the clim_bin()
#' save_path_TPD<-"F:/whaleshark_sdm/4D/TPD" ## new path for storing data of Preferred Depth (PD)
#' save_path_TTD<-"F:/whaleshark_sdm/4D/TTD" ## new path for storing data of Tolerable Depth (TD)
#' TTD_quarter<-TB_depth(temp_bin_path, TTR, time_type = "quarter",
#'                      layer_option="top",
#'                      save_path_TTD, na_value = -9999,
#'                         use_makima = TRUE,
#'                          N_cells = 50)
#' }
#' @export
TB_depth <- function(temp_bin_path, PTR, time_type = c("month", "quarter", "season"),
                     layer_option = "top", save_path, na_value = -9999,
                     use_makima = TRUE, N_cells = 50, ncores = 1) {

  # ----Define makima_interp function----
  makima_interp <- function(x, y, xout) {
    # ---- Input validation ----
    if (length(x) < 2) stop("At least two points are required.")
    if (length(x) != length(y)) stop("x and y must have the same length.")

    # Sort data by x (needed for interval searching)
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    n <- length(x)

    # ---- Compute slopes at each data point (Akima's geometric weighting) ----
    slopes <- numeric(n)

    # Endpoint slopes: use the slope of the first/last segment
    slopes[1] <- (y[2] - y[1]) / (x[2] - x[1])
    slopes[n] <- (y[n] - y[n-1]) / (x[n] - x[n-1])

    # Interior slopes: weighted average of adjacent secant slopes
    if (n >= 3) {
      for (i in 2:(n-1)) {
        m_left  <- (y[i] - y[i-1]) / (x[i] - x[i-1])
        m_right <- (y[i+1] - y[i]) / (x[i+1] - x[i])

        # Weights based on local curvature
        w_left  <- abs(m_right - m_left)
        w_right <- abs(m_left - m_right)

        if (w_left + w_right > 0) {
          slopes[i] <- (w_left * m_left + w_right * m_right) / (w_left + w_right)
        } else {
          slopes[i] <- (m_left + m_right) / 2
        }
      }
    }

    # ---- Interpolate at each query point ----
    result <- numeric(length(xout))

    for (k in seq_along(xout)) {
      xi <- xout[k]

      # Linear extrapolation outside the data range
      if (xi <= x[1]) {
        result[k] <- y[1] + slopes[1] * (xi - x[1])
      } else if (xi >= x[n]) {
        result[k] <- y[n] + slopes[n] * (xi - x[n])
      } else {
        # Find the interval containing xi
        idx <- findInterval(xi, x)
        if (idx < 1) idx <- 1
        if (idx >= n) idx <- n - 1

        # Cubic Hermite interpolation (standard formula)
        h <- x[idx+1] - x[idx]          # interval width
        t <- (xi - x[idx]) / h          # normalized position (0..1)

        # Hermite basis functions
        t2 <- t * t
        t3 <- t2 * t
        h00 <-  2*t3 - 3*t2 + 1
        h10 <-   t3 - 2*t2 + t
        h01 <- -2*t3 + 3*t2
        h11 <-   t3 - t2

        # Evaluate the cubic Hermite interpolant
        result[k] <- h00 * y[idx] + h10 * slopes[idx] * h +
          h01 * y[idx+1] + h11 * slopes[idx+1] * h
      }
    }

    return(result)
  }
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
  if (use_makima) {
    makima_available <- TRUE
    cat("Using makima for interpolation (when data points >= 4)\n")
  } else if (use_makima) {
    cat("Warning: using base R spline interpolation\n")
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

  # ---- Process each time group ----
  result_files <- list()

  # Plotting parameters (unchanged)
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

  for (d_idx in seq_along(valid_dirs)) {
    time_dir <- valid_dirs[d_idx]
    group_name <- valid_group_names[d_idx]
    cat(sprintf("\nProcessing time group %s...\n", group_name))
    group_save_path <- file.path(save_path, group_name)
    if (!dir.exists(group_save_path)) dir.create(group_save_path, recursive = TRUE)

    # ---- Read all depth layers as a raster stack ----
    layer_files <- list.files(time_dir, pattern = "temp_.*\\.asc$", full.names = TRUE)
    if (length(layer_files) == 0) {
      cat(sprintf(" Group %s has no valid data, skipping\n", group_name))
      next
    }

    file_depths <- as.numeric(stringr::str_extract(basename(layer_files), "\\d+\\.?\\d*(?=\\.asc)"))
    depth_order <- order(file_depths)
    layer_files <- layer_files[depth_order]
    file_depths <- file_depths[depth_order]

    # Read all layers into a SpatRaster stack
    cat(" Reading depth layers into stack...\n")
    raster_stack <- terra::rast(layer_files)
    names(raster_stack) <- as.character(file_depths)

    # Replace NA values
    raster_stack[raster_stack == na_value] <- NA

    # ---- Extract values matrix (cells x layers) ----
    cat(" Extracting values matrix...\n")
    values_matrix <- terra::values(raster_stack)  # n_cells x n_layers
    n_cells <- nrow(values_matrix)
    n_layers <- ncol(values_matrix)
    cat(sprintf(" Matrix dimensions: %d cells x %d depth layers\n", n_cells, n_layers))

    # ---- Pre-allocate result vectors ----
    min_depth <- rep(NA_real_, n_cells)
    max_depth <- rep(NA_real_, n_cells)
    range_depth <- rep(NA_real_, n_cells)
    method_used <- rep(NA_character_, n_cells)
    plot_data <- list()

    # ---- Find cells with any valid data ----
    valid_cells <- which(rowSums(!is.na(values_matrix)) > 0)
    n_valid_cells_total <- length(valid_cells)
    cat(sprintf(" Cells with at least one valid depth layer: %d\n", n_valid_cells_total))

    if (n_valid_cells_total == 0) {
      cat(sprintf(" Group %s: no valid cells, skipping\n", group_name))
      next
    }

    # ---- Pre-compute fine depths for interpolation ----
    # Use global depth range for all cells (can be refined per cell if needed)
    global_depth_min <- min(file_depths, na.rm = TRUE)
    global_depth_max <- max(file_depths, na.rm = TRUE)
    fine_depths_global <- seq(global_depth_min, global_depth_max, length.out = 200)

    # ---- Process cells with progress bar ----
    cat(" Interpolating and calculating depth ranges...\n")

    # Use pbapply for progress bar, or fallback to lapply with txtProgressBar
    if (requireNamespace("pbapply", quietly = TRUE) && ncores > 1) {
      # Parallel processing with pbapply
      cl <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)

      results <- pbapply::pblapply(valid_cells, function(cell) {
        temp_vals <- values_matrix[cell, ]
        valid_idx <- which(!is.na(temp_vals))
        n_valid <- length(valid_idx)

        if (n_valid == 0) {
          return(list(min = NA, max = NA, range = NA, method = NA, plot_info = NULL))
        }

        valid_depths <- file_depths[valid_idx]
        valid_temps <- temp_vals[valid_idx]

        # Check if any temps are within PTR
        if (max(valid_temps) < PTR[1] || min(valid_temps) > PTR[2]) {
          return(list(min = NA, max = NA, range = NA, method = NA, plot_info = NULL))
        }

        interp_method <- "unknown"
        min_val <- max_val <- range_val <- NA_real_
        interp_depths <- NULL
        interp_temps <- NULL

        # Case 1: Single point
        if (n_valid == 1) {
          if (valid_temps[1] >= PTR[1] && valid_temps[1] <= PTR[2]) {
            min_val <- max_val <- range_val <- valid_depths[1]
            interp_method <- "single_point"
          }
          return(list(min = min_val, max = max_val, range = range_val,
                      method = interp_method, plot_info = list(
                        cell = cell, depths = valid_depths, temps = valid_temps,
                        method = interp_method, interp_depths = NULL, interp_temps = NULL
                      )))
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
            min_val <- selected$min
            max_val <- selected$max
            range_val <- selected$max - selected$min
          }

          return(list(min = min_val, max = max_val, range = range_val,
                      method = interp_method, plot_info = list(
                        cell = cell, depths = valid_depths, temps = valid_temps,
                        method = interp_method, interp_depths = fine_depths, interp_temps = pred_temps
                      )))
        }

        # Case 3: >=4 points -> try makima, fallback to spline, then linear
        depth_min <- min(valid_depths)
        depth_max <- max(valid_depths)
        fine_depths <- seq(depth_min, depth_max, length.out = 200)
        pred_temps <- NULL

        if (makima_available) {
          tryCatch({
            pred_temps <- makima_interp(valid_depths, valid_temps, fine_depths)

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
          min_val <- selected$min
          max_val <- selected$max
          range_val <- selected$max - selected$min
        }

        return(list(min = min_val, max = max_val, range = range_val,
                    method = interp_method, plot_info = list(
                      cell = cell, depths = valid_depths, temps = valid_temps,
                      method = interp_method, interp_depths = fine_depths, interp_temps = pred_temps
                    )))
      })

      # Extract results
      for (idx in seq_along(valid_cells)) {
        cell <- valid_cells[idx]
        res <- results[[idx]]
        min_depth[cell] <- res$min
        max_depth[cell] <- res$max
        range_depth[cell] <- res$range
        method_used[cell] <- res$method
        if (N_cells > 0 && !is.null(res$plot_info) && !is.na(res$method)) {
          plot_data[[length(plot_data) + 1]] <- res$plot_info
        }
      }

    } else {
      # Sequential processing with txtProgressBar
      pb <- utils::txtProgressBar(min = 0, max = length(valid_cells), style = 3)

      for (idx in seq_along(valid_cells)) {
        utils::setTxtProgressBar(pb, idx)
        cell <- valid_cells[idx]
        temp_vals <- values_matrix[cell, ]
        valid_idx <- which(!is.na(temp_vals))
        n_valid <- length(valid_idx)

        if (n_valid == 0) next

        valid_depths <- file_depths[valid_idx]
        valid_temps <- temp_vals[valid_idx]

        if (max(valid_temps) < PTR[1] || min(valid_temps) > PTR[2]) next

        interp_method <- "unknown"
        min_val <- max_val <- range_val <- NA_real_
        interp_depths <- NULL
        interp_temps <- NULL

        # Case 1: Single point
        if (n_valid == 1) {
          if (valid_temps[1] >= PTR[1] && valid_temps[1] <= PTR[2]) {
            min_val <- max_val <- range_val <- valid_depths[1]
            interp_method <- "single_point"
          }
          min_depth[cell] <- min_val
          max_depth[cell] <- max_val
          range_depth[cell] <- range_val
          method_used[cell] <- interp_method
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
            min_val <- selected$min
            max_val <- selected$max
            range_val <- selected$max - selected$min
          }

          min_depth[cell] <- min_val
          max_depth[cell] <- max_val
          range_depth[cell] <- range_val
          method_used[cell] <- interp_method
          if (N_cells > 0) {
            plot_data[[length(plot_data) + 1]] <- list(
              cell = cell, depths = valid_depths, temps = valid_temps,
              method = interp_method, interp_depths = fine_depths, interp_temps = pred_temps
            )
          }
          next
        }

        # Case 3: >=4 points
        depth_min <- min(valid_depths)
        depth_max <- max(valid_depths)
        fine_depths <- seq(depth_min, depth_max, length.out = 200)
        pred_temps <- NULL

        if (makima_available) {
          tryCatch({
            pred_temps <- makima_interp(valid_depths, valid_temps, fine_depths)

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
          min_val <- selected$min
          max_val <- selected$max
          range_val <- selected$max - selected$min
        }

        min_depth[cell] <- min_val
        max_depth[cell] <- max_val
        range_depth[cell] <- range_val
        method_used[cell] <- interp_method
        if (N_cells > 0) {
          plot_data[[length(plot_data) + 1]] <- list(
            cell = cell, depths = valid_depths, temps = valid_temps,
            method = interp_method, interp_depths = fine_depths, interp_temps = pred_temps
          )
        }
      }
      close(pb)
    }

    # ---- Create rasters from vectors ----
    min_depth_raster <- raster_stack[[1]]
    terra::values(min_depth_raster) <- min_depth
    max_depth_raster <- raster_stack[[1]]
    terra::values(max_depth_raster) <- max_depth
    range_depth_raster <- raster_stack[[1]]
    terra::values(range_depth_raster) <- range_depth

    n_valid_cells <- length(which(!is.na(min_depth)))
    cat(sprintf(" Group %s: %d cells have valid depth ranges\n", group_name, n_valid_cells))

    # ---- Save rasters ----
    min_depth_raster[is.na(min_depth_raster)] <- na_value
    max_depth_raster[is.na(max_depth_raster)] <- na_value
    range_depth_raster[is.na(range_depth_raster)] <- na_value

    min_file <- file.path(group_save_path, paste0(group_name, "_min_depth.asc"))
    max_file <- file.path(group_save_path, paste0(group_name, "_max_depth.asc"))
    range_file <- file.path(group_save_path, paste0(group_name, "_depth_range.asc"))

    terra::writeRaster(min_depth_raster, min_file, filetype = "AAIGrid",
                       overwrite = TRUE, NAflag = na_value)
    terra::writeRaster(max_depth_raster, max_file, filetype = "AAIGrid",
                       overwrite = TRUE, NAflag = na_value)
    terra::writeRaster(range_depth_raster, range_file, filetype = "AAIGrid",
                       overwrite = TRUE, NAflag = na_value)

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
