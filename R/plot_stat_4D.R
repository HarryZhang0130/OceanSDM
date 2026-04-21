#' Plot depth statistics and compute habitat volume (quarterly, monthly, or seasonal)
#'
#' @param file_list Character vector. File paths to raster files. For quarterly/seasonal:
#'   12 files (Q1_min, Q1_max, Q1_range, Q2_min, ...; or Spring_min, ...).
#'   For monthly: 36 files (Jan_min, Jan_max, Jan_range, Feb_min, ...).
#' @param save_path Character. Directory path to save plots and report.
#' @param time_unit Character. One of `"quarter"`, `"month"`, or `"season"` (default: `"quarter"`).
#'
#' @return A list with the following components:
#'   \item{volum}{Data frame containing temporal bins and computed volume (m³).}
#'   \item{depth_range}{Combined plot of range maps (returned by `cowplot::plot_grid`).}
#'   \item{depth_min}{Combined plot of minimum depth maps.}
#'   \item{depth_max}{Combined plot of maximum depth maps.}
#'
#'
#' @examples
#' \donttest{
#' # Define file lists for thermally preferred depth (TPD) and tolerable depth (TTD)
#' TPD_quarter <- list(
#'   "F:/whaleshark_sdm/4D/TPD/Q1/Q1_min_depth.asc",
#'   "F:/whaleshark_sdm/4D/TPD/Q1/Q1_max_depth.asc",
#'   "F:/whaleshark_sdm/4D/TPD/Q1/Q1_depth_range.asc",
#'   "F:/whaleshark_sdm/4D/TPD/Q2/Q2_min_depth.asc",
#'   "F:/whaleshark_sdm/4D/TPD/Q2/Q2_max_depth.asc",
#'   "F:/whaleshark_sdm/4D/TPD/Q2/Q2_depth_range.asc",
#'   "F:/whaleshark_sdm/4D/TPD/Q3/Q3_min_depth.asc",
#'   "F:/whaleshark_sdm/4D/TPD/Q3/Q3_max_depth.asc",
#'   "F:/whaleshark_sdm/4D/TPD/Q3/Q3_depth_range.asc",
#'   "F:/whaleshark_sdm/4D/TPD/Q4/Q4_min_depth.asc",
#'   "F:/whaleshark_sdm/4D/TPD/Q4/Q4_max_depth.asc",
#'   "F:/whaleshark_sdm/4D/TPD/Q4/Q4_depth_range.asc"
#' )
#'
#' TTD_quarter <- list(
#'   "F:/whaleshark_sdm/4D/TTD/Q1/Q1_min_depth.asc",
#'   "F:/whaleshark_sdm/4D/TTD/Q1/Q1_max_depth.asc",
#'   "F:/whaleshark_sdm/4D/TTD/Q1/Q1_depth_range.asc",
#'   "F:/whaleshark_sdm/4D/TTD/Q2/Q2_min_depth.asc",
#'   "F:/whaleshark_sdm/4D/TTD/Q2/Q2_max_depth.asc",
#'   "F:/whaleshark_sdm/4D/TTD/Q2/Q2_depth_range.asc",
#'   "F:/whaleshark_sdm/4D/TTD/Q3/Q3_min_depth.asc",
#'   "F:/whaleshark_sdm/4D/TTD/Q3/Q3_max_depth.asc",
#'   "F:/whaleshark_sdm/4D/TTD/Q3/Q3_depth_range.asc",
#'   "F:/whaleshark_sdm/4D/TTD/Q4/Q4_min_depth.asc",
#'   "F:/whaleshark_sdm/4D/TTD/Q4/Q4_max_depth.asc",
#'   "F:/whaleshark_sdm/4D/TTD/Q4/Q4_depth_range.asc"
#' )
#'
#' save_path_TPD <- "F:/whaleshark_sdm/4D/TPD"
#' save_path_TTD <- "F:/whaleshark_sdm/4D/TTD"
#'
#' # Run the function for preferred depth and tolerable depth
#' volum_TPD <- plot_stat_4D(TPD_quarter, save_path_TPD, time_unit = "quarter")
#' volum_TTD <- plot_stat_4D(TTD_quarter, save_path_TTD, time_unit = "quarter")
#' }
#' \donttest{
#' # Assuming raster files exist in a directory
#' file_list <- list.files("path/to/rasters", full.names = TRUE)
#' plot_stat_4D(file_list, save_path = "./output", time_unit = "quarter")
#' }
#' @export
plot_stat_4D <- function(file_list, save_path, time_unit = "quarter") {

  # --- Validate time_unit and set number of time points and labels ---
  if (time_unit == "quarter") {
    n_time <- 4
    time_labels <- c("Q1", "Q2", "Q3", "Q4")
  } else if (time_unit == "month") {
    n_time <- 12
    time_labels <- month.abb # Jan, Feb, ..., Dec
  } else if (time_unit == "season") {
    n_time <- 4
    time_labels <- c("Spring", "Summer", "Autumn", "Winter")
  } else {
    stop("time_unit must be one of 'quarter', 'month', or 'season'.")
  }

  n_files <- n_time * 3 # each time point has min, max, range
  if (length(file_list) != n_files) {
    stop(sprintf("For time_unit='%s', file_list must have %d files (got %d).",
                 time_unit, n_files, length(file_list)))
  }

  # Create save directory if it does not exist
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

  # Read rasters (no package loading)
  cat(sprintf("Reading %d raster files...\n", n_files))
  data_list <- lapply(file_list, function(f) {
    cat(" Reading", basename(f), "\n")
    raster::raster(f)
  })
  cat("All rasters loaded.\n")

  stats_labels <- c("Minimum", "Maximum", "Range")

  # ---- Volume calculation for each time point (using Range rasters) ----
  cat("\nComputing habitat volume for each", time_unit, "...\n")
  volume_results <- numeric(n_time)
  range_indices <- seq(3, n_files, by = 3) # positions of Range files
  for (i in 1:n_time) {
    idx <- range_indices[i]
    r <- data_list[[idx]]
    cat(" Processing", time_labels[i], "depth_range raster...\n")
    proj <- raster::projection(r)
    is_longlat <- grepl("longlat", proj, ignore.case = TRUE) || proj == ""
    if (is_longlat) {
      newproj <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
      r_proj <- tryCatch({
        raster::projectRaster(r, crs = newproj)
      }, error = function(e) {
        warning("Projection failed for ", time_labels[i], ". Using original raster with area approximation.")
        return(r)
      })
    } else {
      r_proj <- r
    }
    res_m <- raster::res(r_proj)
    cell_area <- res_m[1] * res_m[2]
    total_depth <- raster::cellStats(r_proj, 'sum', na.rm = TRUE)
    volume <- total_depth * cell_area
    volume_results[i] <- volume
    cat(" Volume:", format(volume, scientific = FALSE), "m³\n")
  }

  volum_data <- data.frame("Temporal_bin" = time_labels, "Volume" = volume_results)
  cat("\n=== Habitat Volume by", time_unit, "===\n")
  for (i in 1:n_time) cat(time_labels[i], ":", format(volume_results[i], scientific = FALSE), "cubic meters\n")
  cat("================================\n\n")

  # ---- Land background preparation (only for geographic projection) ----
  proj <- raster::projection(data_list[[1]])
  is_longlat <- grepl("longlat", proj, ignore.case = TRUE) || proj == ""
  land_union <- NULL
  if (is_longlat) {
    cat("Geographic projection detected. Preparing land background...\n")
    tryCatch({
      land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
      land_union <- sf::st_union(land)
      cat(" Land data prepared.\n")
    }, error = function(e) {
      warning("Could not download land data. Proceeding without land background.")
      cat(" Land data not available. Continuing without land background.\n")
    })
  } else {
    warning("Raster projection is not geographic (longlat). Land background disabled.")
  }

  # ---- Build data frames for boxplots and maps ----
  cat("Extracting raster values...\n")
  time_rep <- rep(time_labels, each = 3)
  stat_rep <- rep(stats_labels, times = n_time)

  values_list <- lapply(data_list, function(r) {
    vals <- raster::values(r)
    vals[!is.na(vals)]
  })
  df_box <- data.frame(
    value = unlist(values_list),
    time = factor(rep(time_rep, times = sapply(values_list, length)), levels = time_labels),
    statistic = factor(rep(stat_rep, times = sapply(values_list, length)), levels = stats_labels)
  )

  raster_dfs <- list()
  for (i in seq_along(data_list)) {
    r <- data_list[[i]]
    df <- raster::as.data.frame(r, xy = TRUE)
    colnames(df) <- c("x", "y", "value")
    df$time <- time_rep[i]
    df$statistic <- stat_rep[i]
    raster_dfs[[i]] <- df
  }
  df_map <- do.call(rbind, raster_dfs)

  # ---- Color limits for each statistic ----
  col_low <- "lightblue"; col_high <- "darkblue"
  limits_list <- list()
  for (st in stats_labels) {
    vals <- df_map[df_map$statistic == st, "value"]
    limits_list[[st]] <- range(vals, na.rm = TRUE)
  }

  # ---- Map plotting function (unchanged) ----
  plot_time_map <- function(stat_name, limits, n_time, time_labels) {
    cat(" Creating map panels for", stat_name, "...\n")
    df_sub <- df_map[df_map$statistic == stat_name, ]
    ncol <- ceiling(sqrt(n_time))
    nrow <- ceiling(n_time / ncol)
    while (ncol * nrow < n_time) nrow <- nrow + 1
    plots <- list()
    for (i in 1:n_time) {
      t <- time_labels[i]
      df_t <- df_sub[df_sub$time == t, ]
      raster_idx <- (i - 1) * 3 + which(stats_labels == stat_name)
      r_t <- data_list[[raster_idx]]
      xlim <- c(raster::xmin(r_t), raster::xmax(r_t))
      ylim <- c(raster::ymin(r_t), raster::ymax(r_t))
      p <- ggplot2::ggplot() +
        { if (!is.null(land_union)) ggplot2::geom_sf(data = land_union, fill = "gray80", color = NA, size = 0) } +
        ggplot2::geom_raster(data = df_t, ggplot2::aes(x = x, y = y, fill = value)) +
        ggplot2::scale_fill_gradient(low = col_low, high = col_high, limits = limits,
                                     name = stat_name, na.value = "transparent") +
        ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
        ggplot2::labs(title = t, x = NULL, y = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none",
                       panel.grid = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5, size = 10))
      plots[[i]] <- p
    }
    combined <- cowplot::plot_grid(plotlist = plots, ncol = ncol, nrow = nrow,
                                   labels = paste0("(", letters[1:n_time], ")"),
                                   label_size = 12, label_x = 0, label_y = 1)
    legend_p <- ggplot2::ggplot(df_sub, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_gradient(low = col_low, high = col_high, limits = limits,
                                   name = stat_name, na.value = "transparent") +
      ggplot2::theme(legend.position = "right")
    legend <- ggpubr::get_legend(legend_p)
    combined_with_legend <- cowplot::plot_grid(combined, legend, ncol = 2, rel_widths = c(1, 0.2))
    return(combined_with_legend)
  }

  p_min <- plot_time_map("Minimum", limits_list[["Minimum"]], n_time, time_labels)
  p_max <- plot_time_map("Maximum", limits_list[["Maximum"]], n_time, time_labels)
  p_range <- plot_time_map("Range", limits_list[["Range"]], n_time, time_labels)
  print(p_min); print(p_max); print(p_range)
  suffix <- paste0("_", time_unit, ".tif")
  ggplot2::ggsave(file.path(save_path, paste0("Minimum_map", suffix)), plot = p_min,
                  width = 8, height = 6.8, dpi = 300, device = "tiff", compression = "lzw")
  ggplot2::ggsave(file.path(save_path, paste0("Maximum_map", suffix)), plot = p_max,
                  width = 8, height = 6.8, dpi = 300, device = "tiff", compression = "lzw")
  ggplot2::ggsave(file.path(save_path, paste0("Range_map", suffix)), plot = p_range,
                  width = 8, height = 6.8, dpi = 300, device = "tiff", compression = "lzw")

  # ---- Boxplots and statistical tests (MODIFIED: font size = 8 for all text) ----
  cat("Generating boxplots with significance tests...\n")
  boxplots <- list()
  test_results <- list()

  for (st in stats_labels) {
    cat(" Processing", st, "...\n")
    df_stat <- df_box[df_box$statistic == st, ]
    kw <- stats::kruskal.test(value ~ time, data = df_stat)
    p_val <- kw$p.value

    if (is.na(p_val) || !is.finite(p_val)) {
      cat(" Kruskal-Wallis p-value is NaN/NA. Skipping post-hoc tests.\n")
      p_val <- 1
      p <- ggpubr::ggboxplot(df_stat, x = "time", y = "value",
                             fill = "time", palette = "Blues",
                             xlab = time_unit, ylab = st) +
        ggplot2::theme_minimal() +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "none",
                       text = ggplot2::element_text(size = 8))  # font size 8
      if (n_time > 6) {
        p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))
      }
      boxplots[[st]] <- p
      test_results[[st]] <- list(kruskal_test = kw, dunn_test = NULL)
      next
    }

    cat(" Kruskal-Wallis p-value:", format(p_val, digits = 4), "\n")
    p <- ggpubr::ggboxplot(df_stat, x = "time", y = "value",
                           fill = "time", palette = "Blues",
                           xlab = time_unit, ylab = st) +
      ggplot2::theme_minimal() +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 8))  # font size 8
    if (n_time > 6) {
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))
    }

    if (p_val < 0.05) {
      dunn <- df_stat |>
        rstatix::dunn_test(value ~ time, p.adjust.method = "bonferroni")
      cat(" Significant differences detected, running Dunn's post-hoc test.\n")
      groups <- levels(df_stat$time)
      p_mat <- matrix(1, nrow = length(groups), ncol = length(groups),
                      dimnames = list(groups, groups))
      for (j in 1:nrow(dunn)) {
        g1 <- dunn$group1[j]; g2 <- dunn$group2[j]
        p_mat[g1, g2] <- dunn$p.adj[j]
        p_mat[g2, g1] <- dunn$p.adj[j]
      }
      letters <- multcompView::multcompLetters(p_mat)$Letters
      signif_labels <- data.frame(time = names(letters), label = letters)
      y_max <- max(df_stat$value, na.rm = TRUE)
      y_offset <- y_max * 0.05
      p <- p + ggplot2::geom_text(data = signif_labels,
                                  ggplot2::aes(x = time, y = y_max + y_offset, label = label),
                                  inherit.aes = FALSE, size = 8)  # size 8 to match axis text
    } else {
      cat(" No significant differences (p >= 0.05).\n")
    }
    boxplots[[st]] <- p
    test_results[[st]] <- list(kruskal_test = kw,
                               dunn_test = if (p_val < 0.05) dunn else NULL)
  }

  combined_boxplot <- cowplot::plot_grid(plotlist = boxplots, ncol = 1, nrow = 3, labels = NULL)
  ggplot2::ggsave(file.path(save_path, paste0("boxplots_", time_unit, ".tif")), plot = combined_boxplot,
                  width = 8, height = 12, dpi = 300, device = "tiff", compression = "lzw")

  # ---- Save report (unchanged) ----
  report_file <- file.path(save_path, paste0("statistical_tests_", time_unit, ".txt"))
  sink(report_file)
  cat(sprintf("Statistical Test Report for %s Depth Statistics\n", toupper(time_unit)))
  cat("======================================================\n\n")
  for (st in stats_labels) {
    cat("###", toupper(st), "###\n")
    res <- test_results[[st]]
    cat("Kruskal-Wallis test:\n"); print(res$kruskal_test); cat("\n")
    if (!is.null(res$dunn_test)) {
      cat("Dunn's post-hoc test (Bonferroni adjusted):\n"); print(res$dunn_test); cat("\n")
    }
    cat("------------------------------------------------------------\n\n")
  }
  cat("\n### HABITAT VOLUME (", toupper(time_unit), ") ###\n")
  cat("Volume computed from depth_range rasters (projected to Mollweide, m³):\n")
  for (i in 1:n_time) cat(time_labels[i], ": ", format(volume_results[i], scientific = FALSE), " m³\n")
  cat("======================================================\n")
  sink()

  cat("All plots and report saved to:", save_path, "\n")
  cat("Done.\n")

  return(
    list(
      volum = volum_data,
      depth_range = p_range,
      depth_min = p_min,
      depth_max = p_max
    )
  )
}
