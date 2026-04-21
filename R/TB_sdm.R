#' Building SDM for each temporal bin and create related figures
#'
#' @param env_path character. Directory path to the environmental predictors used to calibrate the model.
#' @param prj_path character. Directory path to the environmental predictors used to project the distribution
#'   in another temporal or spatial extents. At least two separated folders for different scenarios.
#' @param occ_path character. Directory path to occurrence table (tab-delimited text).
#' @param eval_occ character. Directory path to independent occurrence data for testing the model.
#' @param res_path character. Directory path to store the model results.
#' @param species character. Name of the column of species name in the occurrence table.
#' @param x character. Name of the column of longitude in the occurrence table.
#' @param y character. Name of the column of latitude in the occurrence table.
#' @param colin_var character. Method to reduce variable collinearity (default NULL). See ENMTML package.
#' @param sp_accessible_area character. Restrict for each species the accessible area (default NULL,
#'   all spatial extent defined by the environmental layer will be used). See ENMTML.
#' @param pseudoabs_method character. Pseudo-absence allocation method. See ENMTML documentation for details.
#' @param pres_abs_ratio numeric. Presence-Absence ratio (values between 0 and 1).
#' @param part character. Partition method for model's validation. Only one method can be chosen.
#' @param algorithm character. Algorithm(s) to construct ecological niche models. See ENMTML.
#' @param thr character. Threshold(s) used for presence-absence predictions. See ENMTML.
#' @param ensemble character. Method(s) used to ensemble different algorithms. See ENMTML.
#'
#' @return A list containing:
#'   \item{var_imp}{Data frame of variable importance across algorithms.}
#'   \item{var_plot}{ggplot2 object of variable importance boxplot.}
#'   \item{res_cur}{Data frame of response curve data (from sdm::rcurve).}
#'   \item{rcur_plot}{List of response curve plots (from sdm::rcurve).}
#' @import ENMTML raster sdm ggplot2
#' @export
#'
#' @examples
#' \donttest{
#' # Define paths for environmental predictors, projection, occurrences, and results
#' d_env <- "F:/whaleshark_sdm/bin/Q1"
#' d_prj <- "F:/whaleshark_sdm/bin/proj"
#' d_occ <- "F:/whaleshark_sdm/occ/R_typus_q1.txt"
#' d_res <- "F:/whaleshark_sdm/result_sdm/Q1"
#'
#' # Run SDM for the first quarter
#' sdm_q1 <- TB_sdm(
#'   env_path = d_env,
#'   prj_path = d_prj,
#'   occ_path = d_occ,
#'   eval_occ = NULL,
#'   res_path = d_res,
#'   species = "Species",
#'   x = "x",
#'   y = "y",
#'   colin_var = NULL,
#'   sp_accessible_area = NULL,
#'   pseudoabs_method = c(method = "ENV_CONST"),
#'   pres_abs_ratio = 0.3,
#'   part = c(method = "BLOCK"),
#'   algorithm = c("BRT", "SVM", "RDF", "MXD"),
#'   thr = c(type = "MAX_TSS"),
#'   ensemble = c(method = "SUP", metric = "TSS")
#' )
#' }
TB_sdm <- function(env_path, prj_path = NULL, occ_path,
                   eval_occ = NULL, res_path, species, x, y,
                   colin_var, sp_accessible_area, pseudoabs_method,
                   pres_abs_ratio, part, algorithm, thr, ensemble) {

  # ---- Run ENMTML model ----
  ENMTML::ENMTML(
    pred_dir = env_path,
    proj_dir = prj_path,
    occ_file = occ_path,
    eval_occ = eval_occ,
    result_dir = res_path,
    sp = species,
    x = x,
    y = y,
    colin_var = colin_var,
    sp_accessible_area = sp_accessible_area,
    pseudoabs_method = pseudoabs_method,
    pres_abs_ratio = pres_abs_ratio,
    part = part,
    imp_var = TRUE,
    algorithm = algorithm,
    thr = thr,
    msdm = NULL,
    ensemble = ensemble
  )

  # ---- Read evaluation table ----
  eval_file <- file.path(res_path, "Evaluation_Table.txt")
  sh <- utils::read.table(eval_file, header = TRUE)
  cat("Model evaluation table:\n")
  print(sh[c(2, 5, 7, 11, 12, 14, 18)])  # AUC, TSS, OR, etc.

  # ---- Load ensemble rasters ----
  PB_sup <- raster::raster(file.path(res_path, "Ensemble", "SUP", "Rhincodon_typus.tif"))
  PB_ab <- raster::raster(file.path(res_path, "Ensemble", "SUP", "MAX_TSS", "Rhincodon_typus.tif"))

  # ---- Plot predicted maps ----
  ncols <- raster::ncol(PB_sup)
  nrows <- raster::nrow(PB_sup)
  width_inches <- 10
  single_height <- width_inches / (ncols / nrows)
  total_height <- single_height * 2
  dpi <- 300

  grDevices::tiff(filename = file.path(res_path, "predicted_maps.tif"),
                  width = width_inches * dpi,
                  height = total_height * dpi,
                  res = dpi,
                  compression = "lzw",
                  pointsize = 12)

  graphics::par(mfrow = c(2, 1))
  graphics::par(mar = c(4, 4, 4, 1))
  raster::plot(PB_sup)
  graphics::mtext("a) Habitat suitability", side = 3, line = 0, adj = 0, cex = 1, font = 2)

  graphics::par(mar = c(4, 4, 1, 1))
  raster::plot(PB_ab, legend = FALSE)
  graphics::mtext("b) Presence map", side = 3, line = 0, adj = 0, cex = 1, font = 2)

  grDevices::dev.off()

  # ---- Plot presence/absence points on predicted map ----
  graphics::par(mfrow = c(1, 1), mar = c(3, 3, 1, 1), las = 1)
  raster::plot(PB_ab,
               main = "Presences (red) vs. pseudo-absences (black) in predicted presence map",
               cex.main = 0.5)
  occ_blocks <- utils::read.table(file.path(res_path, "BLOCK", "OccBlocks.txt"), header = TRUE)
  pres <- subset(occ_blocks[2:3], occ_blocks$PresAbse == 1)
  graphics::points(pres$x, pres$y, col = "red", pch = 3, cex = 0.5)
  abse <- subset(occ_blocks[2:3], occ_blocks$PresAbse == 0)
  graphics::points(abse$x, abse$y, col = "black", pch = 1, cex = 0.5)

  # ---- Read variable importance scores ----
  cat("Start reading variable importance data......\n")
  result_imp <- list()
  for (algo in algorithm) {
    imp_file <- file.path(res_path, "Algorithm", algo,
                          "Response Curves & Variable Importance",
                          "VariableImportance.txt")
    if (file.exists(imp_file)) {
      df <- utils::read.table(imp_file, header = TRUE, row.names = NULL)
      df <- df[, 2:5]
      colnames(df) <- c("Sp", "Al", "Var", "Im")
      result_imp[[algo]] <- df
    }
  }
  var_data <- do.call(rbind, result_imp)
  cat("Variable importance table:\n")
  str(var_data)

  # ---- Variable importance boxplot ----
  var_order <- stats::aggregate(Im ~ Var, data = var_data, FUN = mean)
  var_order <- var_order[order(var_order$Im, decreasing = FALSE), ]
  var_data$Var <- factor(var_data$Var, levels = var_order$Var)

  p1 <- ggplot2::ggplot(var_data, ggplot2::aes(x = Im, y = Var)) +
    ggplot2::geom_boxplot(fill = "lightblue", color = "darkblue", alpha = 0.7) +
    ggplot2::labs(x = "Importance", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 9)) +
    ggplot2::theme_classic() +
    ggplot2::scale_fill_viridis_d()

  n_vars <- length(unique(var_data$Var))
  n_algorithms <- length(unique(var_data$Al))
  base_width <- 4
  base_height <- 4
  if (n_vars > 8) base_height <- base_height + (n_vars - 8) * 1

  output_file <- file.path(res_path, "Variable_Importance.tif")
  cat("Done......variable importance plot was successfully created!\n")
  ggplot2::ggsave(filename = output_file, plot = p1,
                  width = base_width, height = base_height,
                  dpi = 300, units = "in", device = "tiff")

  # ---- Create response curves using sdm package ----
  cat("Create response curves based on package 'sdm'.......\n")
  files <- list.files(env_path, pattern = "\\.asc$", full.names = TRUE)
  env_data <- raster::stack(files)
  raster::projection(env_data) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

  occ_blocks <- utils::read.table(file.path(res_path, "BLOCK", "OccBlocks.txt"), header = TRUE)
  env <- raster::extract(env_data, occ_blocks[, c("x", "y")])
  env <- as.data.frame(env)
  env$pa <- occ_blocks$PresAbse

  d1 <- sdm::sdmData(pa ~ ., train = env)
  sdm_model <- sdm::sdm(pa ~ ., data = d1,
                        methods = c('brt', 'maxent', 'svm', 'rf'))

  # Extract response curves for all variables
  # Note: sdm::rcurve with gg=TRUE returns a list; we take the data component
  rcur_list <- sdm::rcurve(sdm_model, names(env_data), mean = TRUE, confidence = TRUE, gg = T)
  # rcur_list is a list of lists, one per variable; each has a $data element
  rcur_data <- rcur_list$data
  # The plot object (list of ggplot objects) is stored in rcur_list as well, but we keep the first element as representative?
  # Original code returned p2 (the whole list) - we return the list.
  p2 <- rcur_list  # keep the full list of response curve plots

  cat("Done......response curves were successfully created!\n")

  return(list(
    var_imp = var_data,
    var_plot = p1,
    res_cur = rcur_data,
    rcur_plot = p2
  ))
}
