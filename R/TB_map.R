#' Plot SDM predictions (probability and binary maps) for multiple time periods
#'
#' This function reads probability and binary raster files for a set of time
#' periods (e.g., quarters, months, seasons) from a predefined folder structure
#' and creates a composite figure. Each time period occupies one row: left column
#' shows the probability map, right column shows the binary map. Subplots are
#' labeled a) through h) (or up to 24) outside the top-left corner of each plot.
#' Legends are combined horizontally below the main grid. For monthly data (12
#' months), the function automatically splits into three separate figures, each
#' covering four months.
#'
#' Expected folder structure (for a given time period folder, e.g., "Q1"):
#'   - {root_path}/{period_folder}/predicted_probability.asc   (probability map)
#'   - {root_path}/{period_folder}/predicted_presence.asc      (binary map)
#'
#' @param root_path Character. Path to the root directory containing time period subfolders.
#' @param time_type Character. Type of time periods: "quarter", "month", "season".
#'   Default "quarter". Ignored if `time_folders` is provided.
#' @param time_folders Optional character vector of folder names (e.g., c("Q1","Q2","Q3","Q4")).
#'   If provided, overrides `time_type`.
#' @param time_labels Optional character vector of labels for the time periods
#'   (used in plot row labels). Default uses folder names.
#' @param output_width Numeric. Width of the output figure in inches. Default: 12.
#' @param output_height Numeric. Height of the output figure (per figure) in inches.
#'   For monthly splits, each figure will have this height. Default: 16.
#' @param output_filename_base Character. Base name for output TIFF files.
#'   For monthly data, additional suffixes are added. Default: "sdm_plot".
#' @param split_monthly Logical. For monthly data, if TRUE (default), split into
#'   three figures (Jan-Apr, May-Aug, Sep-Dec). If FALSE, produce a single figure
#'   with 12 rows (very tall).
#'
#' @return No return value. The plot(s) are displayed and saved as TIFF files in `root_path`.
#'
#' @importFrom terra rast is.lonlat xmin xmax ymin ymax as.data.frame
#' @importFrom ggplot2 ggplot geom_sf geom_raster aes scale_fill_gradient scale_fill_manual
#'   coord_sf theme_bw theme element_blank labs annotate guide_colorbar guide_legend
#'   unit margin
#' @importFrom rnaturalearth ne_countries
#' @importFrom sf st_union
#' @importFrom cowplot plot_grid ggdraw draw_plot
#' @importFrom grDevices tiff dev.off
#' @export
TB_map <- function(root_path = "F:/whaleshark_sdm/result",
                   time_type = "quarter",
                   time_folders = NULL,
                   time_labels = NULL,
                   output_width = 12,
                   output_height = 16,
                   output_filename_base = "sdm_plot",
                   split_monthly = TRUE) {

  # ----- 1. Determine time periods -----
  if (!is.null(time_folders)) {
    periods <- time_folders
  } else {
    if (time_type == "quarter") {
      periods <- c("Q1", "Q2", "Q3", "Q4")
    } else if (time_type == "month") {
      periods <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    } else if (time_type == "season") {
      periods <- c("Sp", "Sm", "Am", "Wt")   # Spring, Summer, Autumn, Winter
    } else {
      stop("time_type must be 'quarter', 'month', or 'season'")
    }
  }

  if (is.null(time_labels)) {
    time_labels <- periods
  } else {
    if (length(time_labels) != length(periods)) {
      stop("length of time_labels must equal number of periods")
    }
  }

  # ----- 2. Prepare land background (using terra to check projection) -----
  land_union <- NULL
  first_period <- periods[1]
  first_prob_path <- file.path(root_path, first_period, "predicted_probability.asc")
  if (file.exists(first_prob_path)) {
    r_test <- terra::rast(first_prob_path)
    is_longlat <- terra::is.lonlat(r_test)
    if (is_longlat) {
      tryCatch({
        land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
        land_union <- sf::st_union(land)
        message("Land background prepared.")
      }, error = function(e) {
        warning("Could not download land data. Proceeding without land background.")
        land_union <- NULL
      })
    } else {
      warning("Raster projection is not geographic. Land background disabled.")
    }
  } else {
    warning("First probability raster not found; land background disabled.")
  }

  # ----- 3. Helper: create a single map (without legend) using terra -----
  create_map <- function(raster_path, period_label, map_type, land_union) {
    if (is.null(raster_path) || !file.exists(raster_path)) {
      warning("Raster file missing for ", period_label, " ", map_type)
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Data missing", size = 5) +
               ggplot2::theme_void())
    }

    r <- terra::rast(raster_path)
    # Convert to data frame with coordinates
    df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
    colnames(df) <- c("x", "y", "value")
    df <- df[!is.na(df$value), ]

    xmin <- terra::xmin(r)
    xmax <- terra::xmax(r)
    ymin <- terra::ymin(r)
    ymax <- terra::ymax(r)

    p <- ggplot2::ggplot()
    if (!is.null(land_union)) {
      p <- p + ggplot2::geom_sf(data = land_union, fill = "gray80", color = NA)
    }

    if (map_type == "Probability") {
      p <- p +
        ggplot2::geom_raster(data = df, ggplot2::aes(x = x, y = y, fill = value)) +
        ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue",
                                     name = "Probability", na.value = "transparent") +
        ggplot2::coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none",
                       panel.grid = ggplot2::element_blank(),
                       plot.title = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(5, 0, 0, 0, "pt")) +
        ggplot2::labs(x = NULL, y = NULL)
    } else { # Binary
      df$value <- factor(df$value, levels = c(0, 1))
      p <- p +
        ggplot2::geom_raster(data = df, ggplot2::aes(x = x, y = y, fill = value)) +
        ggplot2::scale_fill_manual(values = c("0" = "lightblue", "1" = "darkblue"),
                                   name = "Binary",
                                   labels = c("Absence", "Presence"),
                                   na.value = "transparent",
                                   drop = FALSE) +
        ggplot2::coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none",
                       panel.grid = ggplot2::element_blank(),
                       plot.title = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(5, 0, 0, 0, "pt")) +
        ggplot2::labs(x = NULL, y = NULL)
    }
    return(p)
  }

  # ----- 4. Helper: generate legends (unchanged, uses ggplot2) -----
  prob_data <- data.frame(x = 1:100, y = 1, value = seq(0, 1, length.out = 100))
  prob_legend_plot <- ggplot2::ggplot(prob_data, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Probability") +
    ggplot2::theme(legend.position = "bottom",
                   legend.direction = "horizontal",
                   legend.title.position = "bottom",
                   legend.title = ggplot2::element_text(hjust = 0.5, size = 10),
                   legend.text = ggplot2::element_text(hjust = 0.5, size = 9),
                   legend.key.width = ggplot2::unit(1.5, "cm"),
                   legend.key.height = ggplot2::unit(0.4, "cm"),
                   plot.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "bottom", title.hjust = 0.5,
                                                   barwidth = grid::unit(3, "cm"),
                                                   barheight = grid::unit(0.4, "cm")))
  prob_grob <- ggplot2::ggplotGrob(prob_legend_plot)
  prob_legend <- prob_grob$grobs[[which(sapply(prob_grob$grobs, function(x) x$name) == "guide-box")[1]]]

  binary_data <- data.frame(x = c(0, 0), y = c(0, 0), value = factor(c(0, 1), levels = c(0, 1)))
  binary_legend_plot <- ggplot2::ggplot(binary_data, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile(alpha = 1, width = 0, height = 0) +
    ggplot2::scale_fill_manual(values = c("0" = "lightblue", "1" = "darkblue"),
                               name = "Binary", labels = c("Absence", "Presence")) +
    ggplot2::theme(legend.position = "bottom",
                   legend.direction = "horizontal",
                   legend.title.position = "bottom",
                   legend.title = ggplot2::element_text(hjust = 0.5, size = 10),
                   legend.text = ggplot2::element_text(hjust = 0.5, size = 9),
                   legend.spacing.x = ggplot2::unit(0.3, "cm"),
                   legend.key.width = ggplot2::unit(0.8, "cm"),
                   plot.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(title.position = "bottom", title.hjust = 0.5,
                                                 nrow = 1, byrow = TRUE,
                                                 override.aes = list(fill = c("lightblue", "darkblue"))))
  binary_grob <- ggplot2::ggplotGrob(binary_legend_plot)
  binary_legend <- binary_grob$grobs[[which(sapply(binary_grob$grobs, function(x) x$name) == "guide-box")[1]]]

  #legend_panel <- gridExtra::arrangeGrob(prob_legend, binary_legend, ncol = 2, widths = grid::unit(rep(1, 2), "null"))
  legend_panel <- cowplot::plot_grid(prob_legend, binary_legend, ncol = 2, rel_widths = c(1, 1))
  # ----- 5. Define function to create a figure for a subset of periods -----
  create_figure <- function(period_subset, label_subset, fig_suffix = "") {
    n_periods <- length(period_subset)
    plots <- list()
    for (i in 1:n_periods) {
      period <- period_subset[i]
      label <- label_subset[i]

      prob_path <- file.path(root_path, period, "predicted_probability.asc")
      binary_path <- file.path(root_path, period, "predicted_presence.asc")

      prob_plot <- create_map(prob_path, label, "Probability", land_union)
      binary_plot <- create_map(binary_path, label, "Binary", land_union)

      plots[[2*i - 1]] <- prob_plot
      plots[[2*i]] <- binary_plot
    }

    main_grid <- cowplot::plot_grid(plotlist = plots, ncol = 2, nrow = n_periods,
                                    labels = paste0(letters[1:(2*n_periods)], ")"),
                                    label_size = 12, label_x = 0, label_y = 1)

    combined_plot <- cowplot::ggdraw() +
      cowplot::draw_plot(main_grid, x = 0, y = 0.12, width = 1, height = 0.88) +
      cowplot::draw_plot(legend_panel, x = 0, y = 0.02, width = 1, height = 0.1)

    output_file <- file.path(root_path, paste0(output_filename_base, fig_suffix, ".tif"))
    ggplot2::ggsave(filename = output_file, plot = combined_plot,
                    width = output_width, height = output_height, dpi = 300, device = "tiff",
                    compression = "lzw")
    return(combined_plot)
  }

  # ----- 6. Generate plots according to time_type -----
  if (time_type == "month" && split_monthly && length(periods) == 12) {
    group1_idx <- 1:4
    group2_idx <- 5:8
    group3_idx <- 9:12
    create_figure(periods[group1_idx], time_labels[group1_idx], "_JanApr")
    message("Plot saved to: ", output_filename_base,"_JanApr.tif")
    create_figure(periods[group2_idx], time_labels[group2_idx], "_MayAug")
    message("Plot saved to: ", output_filename_base,"_MayAug.tif")
    create_figure(periods[group3_idx], time_labels[group3_idx], "_SepDec")
    message("Plot saved to: ", output_filename_base,"_SepDec.tif")
  } else {
    create_figure(periods, time_labels, "")
    message("Plot saved to: ", output_filename_base,".tif")
  }

  # return figures
  if (time_type == "month" && split_monthly && length(periods) == 12) {
    fig1 <- create_figure(periods[group1_idx], time_labels[group1_idx], "_JanApr")
    fig2 <- create_figure(periods[group2_idx], time_labels[group2_idx], "_MayAug")
    fig3 <- create_figure(periods[group3_idx], time_labels[group3_idx], "_SepDec")
    return(list(JanApr = fig1, MayAug = fig2, SepDec = fig3))
  } else {
    fig <- create_figure(periods, time_labels, "")
    return(fig)
  }
  message("All plots completed.")
}
