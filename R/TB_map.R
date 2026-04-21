#' Plot quarterly SDM predictions (probability and binary maps) in a 4x2 layout
#'
#' This function reads probability and binary raster files for four quarters (Q1-Q4)
#' from a predefined folder structure and creates a 4x2 composite figure. Each quarter
#' occupies one row: left column shows the probability map, right column shows the
#' binary map. Subplots are labeled a) through h) outside the top-left corner of each
#' plot. Legends are combined horizontally below the main grid. The figure is displayed
#' in R and saved as a TIFF file with customizable dimensions.
#'
#' The expected folder structure is:
#' \itemize{
#'   \item \code{root_path/Q1/Ensemble/SUP/*.tif} (probability map)
#'   \item \code{root_path/Q1/Ensemble/SUP/MAX_TSS/*.tif} (binary map)
#'   \item (similar for Q2, Q3, Q4)
#' }
#' If a directory is missing or contains no .tif file, a warning is issued and a
#' placeholder plot is used.
#'
#' @param root_path Character. Path to the root directory containing Q1, Q2, Q3, Q4 subfolders.
#'                  Default: "F:/whaleshark_sdm/result"
#' @param output_width Numeric. Width of the output figure in inches. Default: 12
#' @param output_height Numeric. Height of the output figure in inches. Default: 16
#' @param output_filename Character. Name of the output TIFF file (without extension).
#'                        Default: "sdm_quarterly"
#'
#' @return No return value. The plot is printed in the current graphics device and
#'         saved as a TIFF file in \code{root_path}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' TB_map(
#'   root_path = "F:/whaleshark_sdm/result",
#'   output_width = 12,
#'   output_height = 16，
#'   output_filename = "sdm_quarterly"
#' )
#' }
TB_map <- function(root_path = "F:/whaleshark_sdm/result",
                   output_width = 12, output_height = 16,
                   output_filename = "sdm_quarterly") {

  # Check root directory
  if (!dir.exists(root_path)) {
    stop("Root directory does not exist: ", root_path)
  }

  # Prepare land background (medium resolution, merged into single polygon)
  land_union <- NULL
  # Check projection of the first raster to see if it's geographic
  first_raster_path <- list.files(file.path(root_path, "Q1", "Ensemble", "SUP"),
                                  pattern = "\\.tif$", full.names = TRUE)[1]
  if (!is.na(first_raster_path) && file.exists(first_raster_path)) {
    r_test <- raster::raster(first_raster_path)
    proj <- raster::projection(r_test)
    is_longlat <- grepl("longlat", proj, ignore.case = TRUE) || proj == ""
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
  }

  # Quarters
  quarters <- c("Q1", "Q2", "Q3", "Q4")

  # Helper function to get the first .tif file in a directory
  get_tif_file <- function(dir_path) {
    if (!dir.exists(dir_path)) {
      warning("Directory does not exist: ", dir_path)
      return(NULL)
    }
    tif_files <- list.files(dir_path, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
    if (length(tif_files) == 0) {
      warning("No .tif file found in: ", dir_path)
      return(NULL)
    }
    if (length(tif_files) > 1) {
      warning("Multiple .tif files found in ", dir_path, ". Using first: ", basename(tif_files[1]))
    }
    return(tif_files[1])
  }

  # Function to read raster and return a ggplot object WITHOUT legend
  create_map <- function(raster_path, quarter, type, land_union) {
    if (is.null(raster_path) || !file.exists(raster_path)) {
      warning("Raster file missing for ", quarter, " ", type)
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Data missing", size = 5) +
               ggplot2::theme_void())
    }
    r <- raster::raster(raster_path)
    df <- raster::as.data.frame(r, xy = TRUE)
    colnames(df) <- c("x", "y", "value")

    # Remove NA values to avoid NA in legends
    df <- df[!is.na(df$value), ]

    # Get extent for cropping
    xmin <- raster::xmin(r)
    xmax <- raster::xmax(r)
    ymin <- raster::ymin(r)
    ymax <- raster::ymax(r)

    # Base plot with land background
    p <- ggplot2::ggplot()
    if (!is.null(land_union)) {
      p <- p + ggplot2::geom_sf(data = land_union, fill = "gray80", color = NA)
    }

    if (type == "Probability") {
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
      # Convert to factor, ensuring only 0 and 1
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

  # Build plot list (all without legends)
  plots <- list()
  idx <- 1
  for (q in quarters) {
    # Probability map
    prob_dir <- file.path(root_path, q, "Ensemble", "SUP")
    prob_file <- get_tif_file(prob_dir)
    prob_plot <- create_map(prob_file, q, "Probability", land_union)
    plots[[idx]] <- prob_plot
    idx <- idx + 1

    # Binary map
    binary_dir <- file.path(root_path, q, "Ensemble", "SUP", "MAX_TSS")
    binary_file <- get_tif_file(binary_dir)
    binary_plot <- create_map(binary_file, q, "Binary", land_union)
    plots[[idx]] <- binary_plot
    idx <- idx + 1
  }

  # Combine the 8 plots in a 4x2 grid (without legends)
  main_grid <- cowplot::plot_grid(plotlist = plots, ncol = 2, nrow = 4,
                                  labels = paste0(letters[1:8], ")"),
                                  label_size = 12, label_x = 0, label_y = 1)

  # --- Create legend figures using visible (but tiny) elements to ensure proper colors ---
  # Probability legend (horizontal color bar)
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
                                                   barwidth = unit(3, "cm"), barheight = unit(0.4, "cm")))
  prob_grob <- ggplot2::ggplotGrob(prob_legend_plot)
  prob_legend <- prob_grob$grobs[[which(sapply(prob_grob$grobs, function(x) x$name) == "guide-box")[1]]]

  # Binary legend: use geom_tile with alpha = 0 so rectangles are invisible but legend still gets fill colors
  binary_data <- data.frame(x = c(0, 0), y = c(0, 0), value = factor(c(0, 1), levels = c(0, 1)))
  binary_legend_plot <- ggplot2::ggplot(binary_data, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile(alpha = 1, width = 0, height = 0) +  # invisible tiles
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

  # Combine the two legend grobs horizontally
  legend_panel <- gridExtra::arrangeGrob(prob_legend, binary_legend, ncol = 2, widths = grid::unit(rep(1, 2), "null"))

  # Use ggdraw to place main grid above and legend panel below
  combined_plot <- cowplot::ggdraw() +
    cowplot::draw_plot(main_grid, x = 0, y = 0.12, width = 1, height = 0.88) +
    cowplot::draw_grob(legend_panel, x = 0, y = 0.02, width = 1, height = 0.1)

  # Display in R
  print(combined_plot)

  # Save as TIFF
  output_file <- file.path(root_path, paste0(output_filename, ".tif"))
  ggplot2::ggsave(filename = output_file, plot = combined_plot,
                  width = output_width, height = output_height, dpi = 300, device = "tiff",
                  compression = "lzw")

  message("Plot saved to: ", output_file)
}
