#' Plot depth statistics and compute habitat volume (quarterly, monthly, or seasonal)
#'
#' @param file_list Character vector. File paths to raster files.
#' @param save_path Character. Directory path to save plots and report.
#' @param time_type Character. One of \code{"quarter"}, \code{"month"}, or \code{"season"}.
#' @param presence_files Character vector (optional). Paths to binary presence rasters.
#' @param plot_overlay Logical. Whether to create volume-area overlay plot.
#' @param sample_fraction Numeric. Fraction of complete cells to sample for statistical tests. Default: 0.2 (20\%).
#'
#' @return A list with the following components:
#' \item{volum}{Data frame containing temporal bins and computed volume (km³).}
#' \item{area}{Data frame containing temporal bins and computed area (km²).}
#' \item{depth_range}{Combined plot of range maps.}
#' \item{depth_min}{Combined plot of minimum depth maps.}
#' \item{depth_max}{Combined plot of maximum depth maps.}
#' \item{overlay_plot}{A \code{ggplot} object of the volume-area overlay plot.}
#'
#' @examples
#' \donttest{
#' TPD_quarter <- list(...)
#' save_path_TPD <- "F:/whaleshark_sdm/4D/TPD"
#' quarter_files <- c(...)
#' result <- plot_stat_4D(TPD_quarter, save_path_TPD, time_type = "quarter",
#'                        presence_files = quarter_files, sample_fraction = 0.2)
#' }
#' @export
plot_stat_4D <- function(file_list, save_path, time_type = "quarter",
                         presence_files = NULL, plot_overlay = TRUE,
                         sample_fraction = 0.2) {

  # ---- Validate time_type and set time labels ----
  if (time_type == "quarter") {
    n_time <- 4
    time_labels <- c("Q1", "Q2", "Q3", "Q4")
  } else if (time_type == "month") {
    n_time <- 12
    time_labels <- month.abb
  } else if (time_type == "season") {
    n_time <- 4
    time_labels <- c("Spring", "Summer", "Autumn", "Winter")
  } else {
    stop("time_type must be one of 'quarter', 'month', or 'season'.")
  }

  n_files <- n_time * 3
  if (length(file_list) != n_files) {
    stop(sprintf("For time_type='%s', file_list must have %d files (got %d).",
                 time_type, n_files, length(file_list)))
  }

  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

  # ---- Validate sample_fraction ----
  if (sample_fraction <= 0 || sample_fraction > 1) {
    stop("sample_fraction must be between 0 and 1.")
  }

  # ---- Read rasters using terra ----
  cat(sprintf("Reading %d raster files...\n", n_files))
  data_list <- lapply(file_list, function(f) {
    cat("  Reading", basename(f), "\n")
    terra::rast(f)
  })
  cat("All rasters loaded.\n")

  stats_labels <- c("Minimum", "Maximum", "Range")

  # ---- Volume calculation ----
  cat("\nComputing habitat volume for each", time_type, "...\n")
  volume_results <- numeric(n_time)
  range_indices <- seq(3, n_files, by = 3)

  for (i in 1:n_time) {
    idx <- range_indices[i]
    r <- data_list[[idx]]
    standard_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    terra::crs(r) <- standard_crs
    cat("  Processing", time_labels[i], "depth_range raster...\n")

    is_longlat <- terra::is.lonlat(r)
    if (is_longlat) {
      newproj <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
      r_proj <- tryCatch({
        terra::project(r, newproj, method = "bilinear")
      }, error = function(e) {
        warning("Projection failed for ", time_labels[i], ". Using original raster with area approximation.")
        return(r)
      })
    } else {
      r_proj <- r
    }

    res_m <- terra::res(r_proj)
    cell_area <- res_m[1] * res_m[2]
    total_depth <- sum(terra::values(r_proj), na.rm = TRUE)
    volume <- total_depth * cell_area
    volume_results[i] <- volume / 1e9
    cat("    Volume:", format(volume, scientific = FALSE), "km³\n")
  }

  volum_data <- data.frame(Temporal_bin = time_labels, Volume_km3 = volume_results)

  cat("\n=== Habitat Volume by", time_type, "===\n")
  for (i in 1:n_time) {
    cat(time_labels[i], ":", format(volume_results[i], scientific = FALSE), "cubic kilometers\n")
  }
  cat("================================\n\n")

  # ---- Land background ----
  r_first <- data_list[[1]]
  is_longlat <- terra::is.lonlat(r_first)
  land_union <- NULL

  if (is_longlat) {
    cat("Geographic projection detected. Preparing land background...\n")
    tryCatch({
      land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
      land_union <- sf::st_union(land)
      cat("  Land data prepared.\n")
    }, error = function(e) {
      warning("Could not download land data. Proceeding without land background.")
    })
  } else {
    warning("Raster projection is not geographic (longlat). Land background disabled.")
  }

  # ---- Build data frames for boxplots and maps ----
  cat("Extracting raster values...\n")
  time_rep <- rep(time_labels, each = 3)
  stat_rep <- rep(stats_labels, times = n_time)
  values_list <- lapply(data_list, function(r) {
    vals <- terra::values(r)
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
    df <- as.data.frame(r, xy = TRUE)
    colnames(df) <- c("x", "y", "value")
    df$time <- time_rep[i]
    df$statistic <- stat_rep[i]
    raster_dfs[[i]] <- df
  }
  df_map <- do.call(rbind, raster_dfs)

  # ---- Color limits for maps ----
  col_low <- "lightblue"
  col_high <- "darkblue"
  limits_list <- list()
  for (st in stats_labels) {
    vals <- df_map[df_map$statistic == st, "value"]
    limits_list[[st]] <- range(vals, na.rm = TRUE)
  }

  # ---- Map plotting function ----
  plot_time_map <- function(stat_name, limits, n_time, time_labels) {
    cat("  Creating map panels for", stat_name, "...\n")
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
      ext <- terra::ext(r_t)
      xlim <- c(ext[1], ext[2])
      ylim <- c(ext[3], ext[4])

      p <- ggplot2::ggplot() +
        { if (!is.null(land_union)) ggplot2::geom_sf(data = land_union, fill = "gray80", color = NA, linewidth = 0) } +
        suppressWarnings(ggplot2::geom_raster(data = df_t, ggplot2::aes(x = x, y = y, fill = value))) +
        ggplot2::scale_fill_gradient(low = col_low, high = col_high, limits = limits,
                                     name = stat_name, na.value = "transparent") +
        ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
        ggplot2::labs(title = t, x = NULL, y = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none",
                       panel.grid = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5, size = 8))
      plots[[i]] <- p
    }

    combined <- suppressWarnings(cowplot::plot_grid(plotlist = plots, ncol = ncol, nrow = nrow,
                                                    labels = paste0("(", letters[1:n_time], ")"),
                                                    label_size = 12, label_x = 0, label_y = 0.98))
    legend_p <- ggplot2::ggplot(df_sub, ggplot2::aes(x = x, y = y, fill = value)) +
      suppressWarnings(ggplot2::geom_raster()) +
      ggplot2::scale_fill_gradient(low = col_low, high = col_high, limits = limits,
                                   name = stat_name, na.value = "transparent") +
      ggplot2::theme(legend.position = "right")
    legend <- suppressWarnings(cowplot::get_legend(legend_p))
    combined_with_legend <- suppressWarnings(cowplot::plot_grid(combined, legend, ncol = 2, rel_widths = c(1, 0.2)))
    return(combined_with_legend)
  }

  p_min <- plot_time_map("Minimum", limits_list[["Minimum"]], n_time, time_labels)
  p_max <- plot_time_map("Maximum", limits_list[["Maximum"]], n_time, time_labels)
  p_range <- plot_time_map("Range", limits_list[["Range"]], n_time, time_labels)

  suffix <- paste0("_", time_type, ".tif")
  suppressWarnings(ggplot2::ggsave(file.path(save_path, paste0("Minimum_map", suffix)),
                                   plot = p_min, width = 8, height = 6.8, dpi = 300, device = "tiff", compression = "lzw"))
  suppressWarnings(ggplot2::ggsave(file.path(save_path, paste0("Maximum_map", suffix)),
                                   plot = p_max, width = 8, height = 6.8, dpi = 300, device = "tiff", compression = "lzw"))
  suppressWarnings(ggplot2::ggsave(file.path(save_path, paste0("Range_map", suffix)),
                                   plot = p_range, width = 8, height = 6.8, dpi = 300, device = "tiff", compression = "lzw"))

  # ---- Boxplots and statistical tests (using the SAME random cells across time) ----
  # Step 1: Find cells with complete data across ALL time periods
  # Step 2: Randomly sample a fraction of these complete cells (default: 20%)
  # Step 3: Use these same-position cells for all statistical tests
  cat("Generating boxplots with significance tests using", sample_fraction * 100, "% of complete cells...\n")

  boxplots <- list()
  test_results <- list()

  for (st in stats_labels) {
    cat("  Processing", st, "...\n")

    # Extract data for the current statistic (with coordinates)
    df_map_st <- df_map[df_map$statistic == st, ]

    if (nrow(df_map_st) == 0) {
      cat("    No data for statistic", st, "\n")
      next
    }

    # ---- Step 1: Find cells that have data in ALL time periods ----
    # Get unique cell coordinates
    all_cells <- unique(df_map_st[, c("x", "y")])

    # For each cell, check if it exists in all time periods
    complete_cells <- all_cells
    for (t in time_labels) {
      df_t <- df_map_st[df_map_st$time == t, c("x", "y")]
      complete_cells <- merge(complete_cells, df_t, by = c("x", "y"))
    }

    if (nrow(complete_cells) == 0) {
      cat("    No cells with complete data across all time periods.\n")
      next
    }

    cat("    Found", nrow(complete_cells), "cells with complete data across all", length(time_labels), "time periods\n")

    # ---- Step 2: Randomly sample a fraction of these complete cells ----
    n_available <- nrow(complete_cells)
    n_draw <- max(1, round(n_available * sample_fraction))
    n_draw <- min(n_draw, n_available)

    set.seed(123)
    sampled_indices <- sample(1:n_available, size = n_draw, replace = FALSE)
    sampled_positions <- complete_cells[sampled_indices, c("x", "y")]

    cat("    Sampled", n_draw, "cell positions (", round(n_draw/n_available*100, 1), "% of complete cells)\n")

    # ---- Step 3: Extract values at the SAME positions across all time periods ----
    sampled_list <- list()
    for (t in time_labels) {
      df_t <- df_map_st[df_map_st$time == t, c("x", "y", "value")]

      # Merge to get values at the sampled positions
      df_merged <- merge(sampled_positions, df_t, by = c("x", "y"), all.x = TRUE)

      # Since we only use complete cells, there should be no NAs
      na_count <- sum(is.na(df_merged$value))
      if (na_count > 0) {
        cat("    Warning:", na_count, "cells have NA in", t, "(should be 0)\n")
      }

      df_merged$time <- t
      sampled_list[[t]] <- df_merged[, c("x", "y", "value", "time")]
    }

    # Combine all time periods
    df_stat <- do.call(rbind, sampled_list)

    # Remove any rows with NA (just in case)
    df_stat <- df_stat[complete.cases(df_stat), ]

    if (nrow(df_stat) == 0) {
      cat("    No complete data across all time periods after filtering.\n")
      next
    }

    # Ensure time is a factor with proper levels
    df_stat$time <- factor(df_stat$time, levels = time_labels)

    cat("    Final sample size:", nrow(df_stat), "cells across", length(time_labels), "groups\n")

    # ---- Kruskal-Wallis test ----
    kw <- stats::kruskal.test(value ~ time, data = df_stat)
    p_val <- kw$p.value

    if (is.na(p_val) || !is.finite(p_val)) {
      cat("    Kruskal-Wallis p-value is NaN/NA. Skipping post-hoc tests.\n")
      p_val <- 1

      p <- ggplot2::ggplot(df_stat, ggplot2::aes(x = time, y = value, fill = time)) +
        ggplot2::geom_boxplot() +
        ggplot2::scale_fill_brewer(palette = "Blues") +
        ggplot2::labs(x = NULL, y = st) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          legend.position = "none",
          text = ggplot2::element_text(family = "sans", size = 18),
          axis.title = ggplot2::element_text(family = "sans", size = 18),
          axis.text = ggplot2::element_text(family = "sans", size = 18),
          axis.text.x = if (n_time > 6) ggplot2::element_text(angle = 45, hjust = 1, family = "sans", size = 18)
          else ggplot2::element_text(family = "sans", size = 18),
          plot.title = ggplot2::element_text(family = "sans", size = 18),
          strip.text = ggplot2::element_text(family = "sans", size = 18)
        )

      if (n_time > 6) {
        p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 18))
      }

      boxplots[[st]] <- p
      test_results[[st]] <- list(kruskal_test = kw, dunn_test = NULL,
                                 complete_cells_found = nrow(complete_cells),
                                 sample_size = nrow(df_stat))
      next
    }

    cat("    Kruskal-Wallis p-value:", format(p_val, digits = 4), "\n")

    # ---- Base boxplot ----
    p <- ggplot2::ggplot(df_stat, ggplot2::aes(x = time, y = value, fill = time)) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_fill_brewer(palette = "Blues") +
      ggplot2::labs(x = NULL, y = st) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(family = "sans", size = 18),
        axis.title = ggplot2::element_text(family = "sans", size = 18),
        axis.text = ggplot2::element_text(family = "sans", size = 18),
        axis.text.x = if (n_time > 6) ggplot2::element_text(angle = 45, hjust = 1, family = "sans", size = 18)
        else ggplot2::element_text(family = "sans", size = 18),
        plot.title = ggplot2::element_text(family = "sans", size = 18),
        strip.text = ggplot2::element_text(family = "sans", size = 18)
      )

    if (n_time > 6) {
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 18))
    }

    # ---- Post-hoc tests if KW is significant ----
    if (p_val < 0.05) {
      dunn <- rstatix::dunn_test(df_stat, value ~ time, p.adjust.method = "bonferroni")
      cat("    Significant differences detected, running Dunn's post-hoc test.\n")

      # Get the groups that actually appear in the data
      groups_present <- levels(df_stat$time)

      # Create p_mat only for groups that exist
      n_groups <- length(groups_present)
      p_mat <- matrix(1, nrow = n_groups, ncol = n_groups,
                      dimnames = list(groups_present, groups_present))

      for (j in 1:nrow(dunn)) {
        g1 <- dunn$group1[j]
        g2 <- dunn$group2[j]
        # Only assign if both groups exist in our matrix
        if (g1 %in% groups_present && g2 %in% groups_present) {
          p_mat[g1, g2] <- dunn$p.adj[j]
          p_mat[g2, g1] <- dunn$p.adj[j]
        }
      }

      # Get compact letter display
      letters <- multcompView::multcompLetters(p_mat)$Letters

      # Create signif_labels with time as factor matching df_stat
      signif_labels <- data.frame(
        time = factor(names(letters), levels = levels(df_stat$time)),
        label = letters
      )

      # Ensure we only keep labels for groups that exist
      signif_labels <- signif_labels[signif_labels$time %in% levels(df_stat$time), ]

      y_max <- max(df_stat$value, na.rm = TRUE)
      y_offset <- y_max * 0.05

      p <- p + ggplot2::geom_text(
        data = signif_labels,
        ggplot2::aes(x = time, y = y_max + y_offset, label = label),
        inherit.aes = FALSE,
        size = 8
      )
    } else {
      cat("    No significant differences (p >= 0.05).\n")
    }

    boxplots[[st]] <- p
    test_results[[st]] <- list(
      kruskal_test = kw,
      dunn_test = if (p_val < 0.05) dunn else NULL,
      complete_cells_found = nrow(complete_cells),
      sample_size = nrow(df_stat)
    )
  }

  combined_boxplot <- suppressWarnings(cowplot::plot_grid(plotlist = boxplots, ncol = 1, nrow = 3, labels = NULL))
  suppressWarnings(ggplot2::ggsave(file.path(save_path, paste0("boxplots_", time_type, ".tif")),
                                   plot = combined_boxplot, width = 8, height = 12, dpi = 300,
                                   device = "tiff", compression = "lzw"))

  # ---- Save statistical report ----
  report_file <- file.path(save_path, paste0("statistical_tests_", time_type, ".txt"))
  sink(report_file)
  cat(sprintf("Statistical Test Report for %s Depth Statistics\n", toupper(time_type)))
  cat("======================================================\n\n")
  cat("Note: All statistical tests were performed on the SAME randomly sampled cell positions across all time periods.\n")
  cat("  Sampling method: 1) Identified cells with complete data in ALL time periods\n")
  cat("                   2) Randomly sampled ", sample_fraction * 100, "% of these complete cells\n")
  cat("                   3) Extracted values from the SAME positions across all time periods\n\n")

  for (st in stats_labels) {
    cat("###", toupper(st), "###\n")
    res <- test_results[[st]]
    if (!is.null(res$complete_cells_found)) {
      cat("Complete cells (data in all periods):", res$complete_cells_found, "\n")
    }
    if (!is.null(res$sample_size)) {
      cat("Final sample size used:", res$sample_size, "\n\n")
    }
    cat("Kruskal-Wallis test:\n")
    print(res$kruskal_test)
    cat("\n")
    if (!is.null(res$dunn_test)) {
      cat("Dunn's post-hoc test (Bonferroni adjusted):\n")
      print(res$dunn_test)
      cat("\n")
    }
    cat("---------------------------------------------------------------------\n\n")
  }

  cat("\n### HABITAT VOLUME (", toupper(time_type), ") ###\n")
  cat("Volume computed from depth_range rasters (projected to Mollweide, km³):\n")
  for (i in 1:n_time) {
    cat(time_labels[i], ":", format(volume_results[i], scientific = FALSE), " km³\n")
  }
  cat("======================================================\n")
  sink()

  # ---- Volume-Area overlay plot ----
  overlay_plot <- NULL

  if (!is.null(presence_files) && plot_overlay) {
    if (length(presence_files) != n_time) {
      warning("length(presence_files) != n_time. Skipping overlay plot.")
    } else {
      cat("\nCalculating area from presence rasters...\n")

      ref_raster <- data_list[[1]]
      ref_crs <- terra::crs(ref_raster)

      presence_rast <- lapply(presence_files, terra::rast)
      moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
      area_km2 <- numeric(n_time)

      for (i in 1:n_time) {
        r <- presence_rast[[i]]

        standard_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        terra::crs(r) <- standard_crs

        r_proj <- terra::project(r, moll_crs, method = "near")
        cells_1 <- which(terra::values(r_proj) == 1)
        res_m <- terra::res(r_proj)
        cell_area_m2 <- res_m[1] * res_m[2]
        total_area_km2 <- length(cells_1) * cell_area_m2 / 1e6

        area_km2[i] <- total_area_km2
        cat("  ", time_labels[i], "area:", format(area_km2[i], scientific = FALSE), "km²\n")
      }

      area_data <- data.frame(Temporal_bin = time_labels, Area_km2 = area_km2)

      choose_factor <- function(x) {
        if (all(x == 0)) return(1)
        max_val <- max(x, na.rm = TRUE)
        e <- floor(log10(max_val))
        factor <- 10^e
        while (max_val / factor > 10) factor <- factor * 10
        while (max_val / factor < 1 && factor > 1) factor <- factor / 10
        return(factor)
      }

      vol_factor <- choose_factor(volum_data$Volume_km3)
      area_factor <- choose_factor(area_data$Area_km2)

      vol_display <- volum_data$Volume_km3 / vol_factor
      area_display <- area_data$Area_km2 / area_factor
      scale_ratio <- max(vol_display) / max(area_display)
      area_transformed <- area_display * scale_ratio

      plot_data <- data.frame(
        Temporal_bin = factor(time_labels, levels = time_labels),
        Volume = vol_display,
        Area = area_transformed
      )

      vol_exp <- log10(vol_factor)
      area_exp <- log10(area_factor)

      vol_label <- bquote("Habitat volume (×" ~ 10^.(vol_exp) ~ " km"^3 ~ ")")
      area_label <- bquote("Area (×" ~ 10^.(area_exp) ~ " km"^2 ~ ")")

      p_overlay <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Temporal_bin)) +
        ggplot2::geom_col(ggplot2::aes(y = Volume), fill = "black",
                          position = ggplot2::position_dodge(0.9), width = 0.7) +
        ggplot2::geom_line(ggplot2::aes(y = Area, group = 1), color = "red", linewidth = 1) +
        ggplot2::geom_point(ggplot2::aes(y = Area), color = "red", size = 3) +
        ggplot2::scale_y_continuous(
          name = vol_label,
          sec.axis = ggplot2::sec_axis(~ . / scale_ratio, name = area_label,
                                       labels = function(x) format(x, digits = 2))
        ) +
        ggplot2::labs(x = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.line = ggplot2::element_line(),
          axis.title.y.right = ggplot2::element_text(color = "red"),
          axis.text.y.right = ggplot2::element_text(color = "red"),
          legend.position = "none"
        ) +
        ggplot2::scale_fill_manual(values = c("Volume" = "gray70"))

      suppressWarnings(ggplot2::ggsave(file.path(save_path, paste0("volume_area_overlay_", time_type, ".tif")),
                                       plot = p_overlay, width = 6, height = 5, dpi = 300,
                                       device = "tiff", compression = "lzw"))
      overlay_plot <- p_overlay
      cat("Overlay plot saved.\n")
    }
  }

  cat("All plots and report saved to:", save_path, "\n")
  cat("Done.\n")

  return(list(
    volum = volum_data,
    area = if (exists("area_data")) area_data else NULL,
    depth_range = p_range,
    depth_min = p_min,
    depth_max = p_max,
    overlay_plot = overlay_plot
  ))
}
