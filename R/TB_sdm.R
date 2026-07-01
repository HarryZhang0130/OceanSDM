#' Build Species Distribution Model using `sdm` and `blockCV` for spatial cross-validation
#'
#' @param env_path character. Directory containing environmental rasters (.asc or .tif).
#' @param occ_path character. Path to occurrence file (tab‑delimited).
#' @param res_path character. Output directory.
#' @param species character. Column name for species (not used, kept for compatibility).
#' @param x character. Column name for longitude.
#' @param y character. Column name for latitude.
#' @param pre_abs_ratio numeric. Ratio of presences to pseudo-absences points.
#' @param collinearity_thresh numeric. Correlation threshold.
#' @param remove_collinear logical. Remove correlated variables.
#' @param max_dist numeric. Maximum distance (unit: meters) to be used in calculating Moran's I. Default is NULL
#' @param block_size numeric. Block size (unit: meters) can be user defined if you did not want to find it through the time-consuming auto-process based on Moran's I curve
#' @param algorithms character vector. sdm methods.
#' @param cv_method character. type of assignment of spatial blocks into folds. Can be random (default), systematic, or checkerboard.
#' @param n_folds integer. Number of folds for spatial cross-validation.
#' @param var.importance logical.
#' @param response.curve logical.
#' @param ensemble.method character. "weighted" or "unweighted".
#' @param ensemble.stat character. "AUC" or "TSS".
#' @param expr character. Select which models will be used for ensemble based on AUC and/or TSS, e.g. expr = 'auc > 0.8 & tss > 0.6'.
#' @param ncore integer. Number of cores.
#' @param seed integer. Random seed.
#' @param proj_path character. Optional projection rasters directory.
#'
#' @return A list containing:
#'   \item{final_model}{A final sdm model object}
#'   \item{ensemble}{Raster object of the predicted habitat suitability by the ensemble model}
#'   \item{evaluation}{Model evaluation table}
#'   \item{var_imp}{Variable importance table}
#'   \item{response_curves}{Data frame of response curve data (from sdm::rcurve)}
#'   \item{optimal_block}{Identified optimal block size}
#'   \item{collinearity_report}{Collinearity analysis report}
#' @import Metrics sf sdm blockCV ggplot2
#' @export
#'
#' @examples
#' \donttest{
#' d_env <- "F:/whaleshark_sdm/bin/Q1"
#' d_prj <- "F:/whaleshark_sdm/bin/proj"
#' d_occ <- "F:/whaleshark_sdm/occ/R_typus_q1.txt"
#' d_res <- "F:/whaleshark_sdm/result_sdm/Q1"
#' sdm_q1 <- TB_sdm(
#'   env_path = d_env,
#'   occ_path = d_occ,
#'   res_path = d_res,
#'   species = "Species",
#'   x = "x",
#'   y = "y",
#'   pre_abs_ratio = 1,
#'   collinearity_thresh = 0.7,
#'   remove_collinear = TRUE,
#'   algorithms = c("brt", "rf", "maxent", "svm"),
#'   cv_method = "random",
#'   n_folds = 5,
#'   var.importance = TRUE,
#'   response.curve = TRUE,
#'   ensemble.method = "weighted",
#'   ensemble.stat = "TSS",
#'   expr = 'auc > 0.8 & tss > 0.6',
#'   ncore = 1,
#'   seed = 123,
#'   proj_path = d_prj)
#' }
TB_sdm <- function(env_path,
                   occ_path,
                   res_path,
                   species,
                   x,
                   y,
                   pre_abs_ratio = 1,
                   collinearity_thresh = 0.7,
                   remove_collinear = FALSE,
                   block_size = NULL,
                   algorithms = c("glm", "gam", "brt", "rf", "svm"),
                   cv_method = "random",
                   n_folds = 5,
                   var.importance = TRUE,
                   response.curve = TRUE,
                   ensemble.method = "weighted",
                   ensemble.stat = "AUC",
                   expr = 'auc > 0.8 & tss > 0.6',
                   ncore = 1,
                   seed = 123,
                   proj_path = NULL) {

  ############### Helper: custom cv_spatial_autocor() ############
  cv_spatial_autocor <- function(x, column, r = NULL, max_dist = NULL, n_dist = 50000, plot = TRUE) {
    # Extract coordinates and values
    coords <- sf::st_coordinates(x)
    values <- x[[column]]

    # Compute pairwise distances
    dist_mat <- as.matrix(stats::dist(coords))

    # Determine maximum distance if not provided
    if (is.null(max_dist)) {
      bbox <- sf::st_bbox(x)
      max_dist <- sqrt((bbox$xmax - bbox$xmin)^2 + (bbox$ymax - bbox$ymin)^2) / 2
    }

    # Create distance bins
    breaks <- seq(0, max_dist, length.out = n_dist + 1)
    distances <- (breaks[-1] + breaks[-(n_dist+1)]) / 2

    # Function to calculate Moran's I for points within a distance threshold
    moran_at_distance <- function(d) {
      # Create binary spatial weights matrix: 1 if distance <= d, 0 otherwise (excluding self)
      w <- ifelse(dist_mat <= d & dist_mat > 0, 1, 0)
      # Row‑standardise weights
      w <- w / rowSums(w)
      w[is.na(w)] <- 0
      # Manual calculation: I = (n / S0) * (z' W z) / (z' z)
      n <- length(values)
      z <- values - mean(values, na.rm = TRUE)
      S0 <- sum(w)
      if (S0 == 0) return(NA)
      numerator <- sum(z * (w %*% z))
      denominator <- sum(z^2)
      I <- (n / S0) * numerator / denominator
      return(I)
    }

    # Compute Moran's I for each distance bin
    autocorr <- sapply(distances, moran_at_distance)

    # Find the distance where autocorrelation drops below threshold (e.g., 0.1)
    threshold <- 0.1
    idx <- which(autocorr < threshold)
    if (length(idx) > 0) {
      recommended <- distances[min(idx)]
    } else {
      recommended <- max_dist * 0.5 # fallback
    }

    # Plot if requested
    if (plot) {
      require(ggplot2)
      p <- ggplot2::ggplot(data.frame(d = distances, I = autocorr), ggplot2::aes(x = d, y = I)) +
        ggplot2::geom_line(color = "steelblue", linewidth = 1) +
        ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
        ggplot2::geom_vline(xintercept = recommended, linetype = "dotted", color = "darkgreen") +
        ggplot2::labs(x = "Distance (map units)", y = "Moran's I",
                      title = "Spatial autocorrelation decay") +
        ggplot2::theme_minimal()
      print(p)
    }

    return(list(
      recommended_size = recommended,
      distances = distances,
      autocorr = autocorr,
      plot = if (plot) p else NULL
    ))
  }

  ########### Helper: custom Boyce index ############
  calc_boyce <- function(pres_prob, bg_prob, n_bins = 10) {
    if (length(pres_prob) == 0 || length(bg_prob) == 0) return(NA_real_)
    pres_prob <- pres_prob[!is.na(pres_prob)]
    bg_prob <- bg_prob[!is.na(bg_prob)]
    if (length(pres_prob) < 2 || length(bg_prob) < 2) return(NA_real_)
    bin_breaks <- seq(min(bg_prob), max(bg_prob), length.out = n_bins + 1)
    bin_centers <- (bin_breaks[-1] + bin_breaks[-(n_bins+1)]) / 2
    F_ratio <- sapply(1:n_bins, function(b) {
      bg_in <- sum(bg_prob >= bin_breaks[b] & bg_prob < bin_breaks[b+1])
      pres_in <- sum(pres_prob >= bin_breaks[b] & pres_prob < bin_breaks[b+1])
      if (bg_in == 0) return(NA)
      (pres_in / length(pres_prob)) / (bg_in / length(bg_prob))
    })
    valid <- !is.na(F_ratio)
    if (sum(valid) < 2) return(NA_real_)
    cor(bin_centers[valid], F_ratio[valid], method = "spearman", use = "pairwise.complete.obs")
  }

  ######### 0. Setup ######################
  if (!dir.exists(res_path)) dir.create(res_path, recursive = TRUE)
  set.seed(seed)
  if (ncore > 1) {
    doParallel::registerDoParallel(cores = ncore)
    on.exit(doParallel::stopImplicitCluster())
  }

  ############# 1. Read environmental rasters and sanitise names ############
  message("Reading environmental rasters...")
  env_files <- list.files(env_path, pattern = "\\.(asc|tif)$", full.names = TRUE)
  if (length(env_files) == 0) stop("No .asc or .tif files in env_path")
  env_stack <- terra::rast(env_files)
  env_names <- names(env_stack)

  ############# 2. Read occurrence data (only presence) #############
  message("Reading occurrence data...")
  occ_raw <- utils::read.table(occ_path, header = TRUE, sep = "\t")
  occ_raw <- occ_raw[, c("x", "y")]
  occ <- unique(occ_raw)
  message("Number of unique presence records: ", nrow(occ))

  ############# 3. Generate environmentally weighted pseudo-absences (eDist) and merge #############
  message("Generating pseudo-absences using sdm::background with method = 'eDist'...")
  bg <- sdm::background(
    x = env_stack,
    n = nrow(occ)/pre_abs_ratio,
    method = 'eDist',
    sp = occ
  )
  bg_pts <- bg[, c("x", "y")]
  bg_pts$pr_ab <- 0
  occ_df <- occ[, c("x", "y")]
  occ_df$pr_ab <- 1
  all_data <- rbind(occ_df, bg_pts)
  all_data <- all_data[complete.cases(all_data), ]
  message(sprintf("Total points (presence + pseudo-absence): %d (presence: %d, pseudo-absence: %d)",
                  nrow(all_data), nrow(occ_df), nrow(bg_pts)))

  ############# 4. Collinearity diagnostics #############
  message("Performing collinearity diagnostics...")
  env_vals <- terra::extract(env_stack, occ[, c("x", "y")], ID = FALSE)
  env_vals <- env_vals[complete.cases(env_vals), ]
  if (nrow(env_vals) < 2) {
    collinearity_report <- data.frame(Var1 = character(), Var2 = character(), Correlation = numeric())
  } else {
    cor_mat <- cor(env_vals, use = "pairwise.complete.obs")
    high_cor <- which(abs(cor_mat) > collinearity_thresh & abs(cor_mat) < 1, arr.ind = TRUE)
    if (nrow(high_cor) > 0) {
      high_cor <- high_cor[high_cor[,1] < high_cor[,2], , drop = FALSE]
      collinearity_report <- data.frame(
        Var1 = rownames(cor_mat)[high_cor[,1]],
        Var2 = colnames(cor_mat)[high_cor[,2]],
        Correlation = cor_mat[high_cor]
      )
      collinearity_report <- collinearity_report[order(-abs(collinearity_report$Correlation)), ]
      message(sprintf("Found %d variable pairs with |cor| > %.2f.", nrow(collinearity_report), collinearity_thresh))
    } else {
      collinearity_report <- data.frame(Var1 = character(), Var2 = character(), Correlation = numeric())
      message("No variable pairs with |cor| > ", collinearity_thresh, " found.")
    }
  }
  utils::write.csv(collinearity_report, file.path(res_path, "collinearity_report.csv"), row.names = FALSE)

  if (remove_collinear && nrow(collinearity_report) > 0) {
    cat("\n=== Collinearity Report ===\n")
    print(collinearity_report)
    cat("===========================\n\n")
    vars_to_remove <- c()
    for (i in 1:nrow(collinearity_report)) {
      v1 <- collinearity_report$Var1[i]
      v2 <- collinearity_report$Var2[i]
      if (v1 %in% vars_to_remove || v2 %in% vars_to_remove) next
      cat(sprintf("Pair %d: %s (cor = %.3f) vs %s\n", i, v1, collinearity_report$Correlation[i], v2))
      cat(" Choose which variable to keep:\n")
      cat(sprintf(" 1: %s\n", v1))
      cat(sprintf(" 2: %s\n", v2))
      choice <- readline(prompt = " Enter number (1 or 2): ")
      while (!(choice %in% c("1", "2"))) {
        choice <- readline(prompt = "Invalid input. Please enter 1 or 2: ")
      }
      if (choice == "1") {
        vars_to_remove <- c(vars_to_remove, v2)
      } else {
        vars_to_remove <- c(vars_to_remove, v1)
      }
    }
    vars_to_remove <- unique(vars_to_remove)
    message("Removing collinear variables: ", paste(vars_to_remove, collapse = ", "))
    keep_vars <- setdiff(env_names, vars_to_remove)
    if (length(keep_vars) == 0) stop("All variables would be removed; adjust threshold or set remove_collinear = FALSE.")
    env_stack <- env_stack[[keep_vars]]
    env_names <- keep_vars
    message("Remaining variables: ", paste(env_names, collapse = ", "))
  } else if (remove_collinear && nrow(collinearity_report) == 0) {
    message("No high correlations detected -- no variables removed.")
  } else {
    message("Variables are kept unchanged for modeling (remove_collinear = FALSE).")
  }

  ############# 5. Create sdmData object with correct formula #############
  message("Creating full sdmData object...")
  env_full <- terra::extract(env_stack, all_data[, c("x", "y")])
  env_full <- as.data.frame(env_full)
  env_full$pr_ab <- all_data$pr_ab
  env_full <- env_full[, -1] # remove ID column
  final_data <- sdm::sdmData(pr_ab ~ ., train = env_full)

  ############# 6. Spatial blocking using blockCV (if requested) #############
  ## covert all_data to a sf object
  all_sf <- sf::st_as_sf(all_data, coords = c("x", "y"), crs = 4326)
  ## projection
  center_lon <- mean(all_data$x, na.rm = TRUE)
  utm_zone <- floor((center_lon + 180) / 6) + 1
  utm_crs <- paste0("+proj=utm +zone=", utm_zone, " +datum=WGS84 +units=m")
  all_proj <- sf::st_transform(all_sf, crs = utm_crs)

  ## calculate block size
  if(is.null(block_size)){
    message("\n Start to find optimal block size based on Moran's I of species data ...")
    result <- cv_spatial_autocor(x = all_proj, column = "pr_ab", max_dist = max_dist, n_dist = 50000, plot = TRUE)
    block_size <- result$recommended_size
    message("Optimal block size = ", paste(block_size), " meters")
  }
  ## create spatial folds
  message("Creating spatial folds using blockCV...")
  if (!requireNamespace("blockCV", quietly = TRUE)) {
    stop("Package 'blockCV' is required for spatial blocking.")
  }
  pts_sf <- sf::st_as_sf(all_data[, c("x", "y")], coords = c("x", "y"), crs = terra::crs(env_stack))
  pts_sf$pr_ab <- all_data$pr_ab
  set.seed(seed)

  if(cv_method == "random"){
    blocks <- blockCV::cv_spatial(
      x = pts_sf,
      column = "pr_ab",
      r = env_stack[[1]],
      size = block_size,
      k = n_folds,
      selection = "random",
      iteration = 50,
      progress = FALSE,
      plot = FALSE
    )
    fold_ids <- blocks$folds_ids
    fold_list <- lapply(1:n_folds, function(f) which(fold_ids == f))
    message("Spatial folds created. Fold sizes: ", paste(sapply(fold_list, length), collapse = ", "))
  }

  if(cv_method == "systematic"){
    blocks <- blockCV::cv_spatial(
      x = pts_sf,
      column = "pr_ab",
      r = env_stack[[1]],
      size = block_size,
      k = n_folds,
      selection = "systematic",
      progress = FALSE,
      plot = FALSE
    )
    fold_ids <- blocks$folds_ids
    fold_list <- lapply(1:n_folds, function(f) which(fold_ids == f))
    message("Spatial folds created. Fold sizes: ", paste(sapply(fold_list, length), collapse = ", "))
  }

  if(cv_method == "checkerboard"){
    blocks <- blockCV::cv_spatial(
      x = pts_sf,
      column = "pr_ab",
      r = env_stack[[1]],
      hexagon = FALSE,
      size = block_size,
      k = n_folds,
      selection = "checkerboard",
      progress = FALSE,
      plot = FALSE
    )
    fold_ids <- blocks$folds_ids
    fold_list <- lapply(1:n_folds, function(f) which(fold_ids == f))
    message("Spatial folds created. Fold sizes: ", paste(sapply(fold_list, length), collapse = ", "))
  }

  ############## 7. Manual cross-validation loop #############
  message("Starting manual cross-validation...")
  algo_names <- algorithms
  perf_list <- list()
  for (algo in algo_names) {
    perf_list[[algo]] <- list(AUC = numeric(n_folds), TSS = numeric(n_folds), Boyce = numeric(n_folds))
  }

  for (fold in 1:n_folds) {
    message(" Processing fold ", fold, "/", n_folds)
    train_idx <- which(fold_ids != fold)
    test_idx <- which(fold_ids == fold)

    train_data <- all_data[train_idx, ]
    test_data <- all_data[test_idx, ]
    # ---- Prepare training data ----
    env_train <- terra::extract(env_stack, train_data[, c("x", "y")])
    env_train <- as.data.frame(env_train)
    env_train <- env_train[, -1, drop = FALSE]  # remove ID
    complete_train <- complete.cases(env_train)
    env_train <- env_train[complete_train, , drop = FALSE] # remove na rows
    train_data <- train_data[complete_train, , drop = FALSE]  # remove na rows
    if (nrow(train_data) == 0) next
    train_sdmdata <- sdm::sdmData(pr_ab ~ ., train = cbind(train_data[, "pr_ab", drop = FALSE], env_train))

    # ---- Prepare test data ----
    env_test <- terra::extract(env_stack, test_data[, c("x", "y")])
    env_test <- as.data.frame(env_test)
    env_test <- env_test[, -1, drop = FALSE]  # remove ID
    complete_test <- complete.cases(env_test)
    env_test <- env_test[complete_test, , drop = FALSE] # remove na rows
    test_data <- test_data[complete_test, , drop = FALSE]  # remove na rows
    if (nrow(test_data) == 0) next


    fold_model <- suppressWarnings(sdm::sdm(
      formula = pr_ab ~ .,
      data = train_sdmdata,
      methods = algorithms,
      ncore = ncore,
      seed = seed
    ))

    for (algo in algo_names) {
      mod_ids <- which(sdm::getModelInfo(fold_model)$method == algo)
      if (length(mod_ids) == 0) {
        warning("Algorithm ", algo, " not found in fold ", fold)
        next
      }

      pred_test <- tryCatch({
        predict(fold_model, newdata = env_test, id = mod_ids, type = "prob")
      }, error = function(e) NULL)
      if (is.null(pred_test)) {
        warning("Prediction failed for ", algo, " in fold ", fold)
        next
      }
      colnames(pred_test) <- "Prob"
      pred <- pred_test$Prob
      obs <- test_data$pr_ab
      valid <- !is.na(pred)
      pred <- pred[valid]
      obs <- obs[valid]
      if (length(unique(obs)) > 1 && !all(is.na(pred_test))) {
        auc_val <- Metrics::auc(obs, pred)
        thresholds <- seq(0, 1, by = 0.01)
        tss_vals <- sapply(thresholds, function(th) {
          pred_class <- ifelse(pred >= th, 1, 0)
          tp <- sum(pred_class == 1 & obs == 1)
          fp <- sum(pred_class == 1 & obs == 0)
          fn <- sum(pred_class == 0 & obs == 1)
          tn <- sum(pred_class == 0 & obs == 0)
          sens <- ifelse((tp+fn)==0, 0, tp/(tp+fn))
          spec <- ifelse((tn+fp)==0, 0, tn/(tn+fp))
          sens + spec - 1
        })
        tss_val <- max(tss_vals, na.rm = TRUE)
        boyce_val <- calc_boyce(pred[obs == 1], pred)
      } else {
        auc_val <- tss_val <- boyce_val <- NA
      }
      perf_list[[algo]]$AUC[fold] <- auc_val
      perf_list[[algo]]$TSS[fold] <- tss_val
      perf_list[[algo]]$Boyce[fold] <- boyce_val
    }
  }

  ############## 8. Summarise cross-validation results #############
  eval_summary <- data.frame()
  for (algo in algo_names) {
    df <- data.frame(
      methodID = algo,
      AUC_mean = mean(perf_list[[algo]]$AUC, na.rm = TRUE),
      AUC_sd = sd(perf_list[[algo]]$AUC, na.rm = TRUE),
      TSS_mean = mean(perf_list[[algo]]$TSS, na.rm = TRUE),
      TSS_sd = sd(perf_list[[algo]]$TSS, na.rm = TRUE),
      Boyce_mean = mean(perf_list[[algo]]$Boyce, na.rm = TRUE),
      Boyce_sd = sd(perf_list[[algo]]$Boyce, na.rm = TRUE)
    )
    eval_summary <- rbind(eval_summary, df)
  }

  ############## 9. Fit final model on all data #############
  message("Fitting final model on all data...")
  final_model <- suppressWarnings(sdm::sdm(
    formula = pr_ab ~ .,
    data = final_data,
    methods = algorithms,
    ncore = ncore,
    seed = seed
  ))

  ############## 10. Ensemble prediction #############
  message("Building ensemble prediction...")
  ens_setting <- list(
    method = ensemble.method,
    stat = ensemble.stat,
    wtest = "training",
    expr = expr
  )
  env_file<-file.path(res_path, "ensemble_suitability.asc")
  ensemble_raster <- suppressWarnings(sdm::ensemble(
    final_model,
    newdata = env_stack,
    setting = ens_setting
  ))
  names(ensemble_raster) <- "ensemble_suitability"
  terra::writeRaster(ensemble_raster,
                     filename = file.path(res_path, "predicted_probability.asc"),
                     overwrite = TRUE, NAflag = -9999, filetype = "AAIGrid")

  best_threshold <- suppressWarnings(sdm::threshold(final_model,
                                   id = 'ensemble',
                                   setting = ens_setting,
                                   opt = 2))
  cat("Best threshold based on MAX_TSS method =",best_threshold,"\n")
  binary_raster <- ensemble_raster >= best_threshold
  binary_raster <- terra::app(binary_raster, fun = function(x) as.integer(x))
  binary_raster <- terra::classify(binary_raster, cbind(NA, -9999))

  names(binary_raster) <- "Presence_SUP_MAXTSS"

    terra::writeRaster(binary_raster,
                     filename = file.path(res_path, "predicted_presence.asc"),
                     overwrite = TRUE, NAflag = -9999, filetype = "AAIGrid")

  ## Evaluate ensemble model using spatial CV folds
  message("Evaluate ensemble model ...")
  auc_vals <- numeric(n_folds)
  tss_vals <- numeric(n_folds)
  boyce_vals <- numeric(n_folds)
  for (fold in 1:n_folds) {
    cat("Processing fold ",fold,"/",n_folds,"\n")
    train_idx <- which(fold_ids != fold)
    test_idx <- which(fold_ids == fold)
    train_data <- all_data[train_idx, ]
    test_data <- all_data[test_idx, ]
    cat("Prepare training data...\n")
    env_train <- terra::extract(env_stack, train_data[, c("x", "y")])
    env_train <- as.data.frame(env_train)
    env_train <- env_train[, -1, drop = FALSE]
    complete_train <- complete.cases(env_train)
    env_train <- env_train[complete_train, , drop = FALSE]
    train_data <- train_data[complete_train, , drop = FALSE]
    if (nrow(train_data) == 0) next
    cat("Create sdmdata to build the model for fold = ",fold,"\n")
    train_sdmdata <- sdm::sdmData(pr_ab ~ ., train = cbind(train_data[, "pr_ab", drop = FALSE], env_train))

    fold_model <- suppressWarnings(sdm::sdm(
      formula = pr_ab ~ .,
      data = train_sdmdata,
      methods = algorithms,
      ncore = ncore,
      seed = seed
    ))
    cat("Prepare test data...\n")
    env_test <- terra::extract(env_stack, test_data[, c("x", "y")])
    env_test <- as.data.frame(env_test)
    env_test <- env_test[, -1, drop = FALSE]
    complete_test <- complete.cases(env_test)
    env_test <- env_test[complete_test, , drop = FALSE]
    test_data <- test_data[complete_test, , drop = FALSE]
    if (nrow(test_data) == 0) next

    ensemble_pred <- suppressWarnings(sdm::ensemble(
      fold_model,
      newdata = env_test,
      setting = ens_setting,
      overwrite = TRUE
    ))
    colnames(ensemble_pred) <- "Prob"
    pred <- ensemble_pred$Prob
    obs <- test_data$pr_ab
    valid <- !is.na(pred)
    pred <- pred[valid]
    obs <- obs[valid]

    #
    if (length(unique(obs)) > 1 && length(pred) > 1 && !all(is.na(pred))) {
      #
      auc_val <- tryCatch(Metrics::auc(obs, pred), error = function(e) NA)
      if (is.na(auc_val)) {
        auc_vals[fold] <- NA
        tss_vals[fold] <- NA
        boyce_vals[fold] <- NA
        next
      }
      auc_vals[fold] <- auc_val

      #
      thresholds <- seq(0, 1, by = 0.01)
      tss_val <- sapply(thresholds, function(th) {
        pred_class <- ifelse(pred >= th, 1, 0)
        tp <- sum(pred_class == 1 & obs == 1)
        fp <- sum(pred_class == 1 & obs == 0)
        fn <- sum(pred_class == 0 & obs == 1)
        tn <- sum(pred_class == 0 & obs == 0)
        sens <- ifelse((tp+fn)==0, 0, tp/(tp+fn))
        spec <- ifelse((tn+fp)==0, 0, tn/(tn+fp))
        sens + spec - 1
      })
      if (all(is.na(tss_val))) {
        tss_vals[fold] <- NA
      } else {
        tss_vals[fold] <- max(tss_val, na.rm = TRUE)
      }

      #
      boyce_vals[fold] <- tryCatch(
        calc_boyce(pred[obs == 1], pred),
        error = function(e) NA
      )
    } else {
      auc_vals[fold] <- NA
      tss_vals[fold] <- NA
      boyce_vals[fold] <- NA
    }
  }

  #
  auc_mean <- mean(auc_vals, na.rm = TRUE)
  auc_sd   <- sd(auc_vals, na.rm = TRUE)
  tss_mean <- mean(tss_vals, na.rm = TRUE)
  tss_sd   <- sd(tss_vals, na.rm = TRUE)
  boyce_mean <- mean(boyce_vals, na.rm = TRUE)
  boyce_sd   <- sd(boyce_vals, na.rm = TRUE)

  #
  ensemble_row <- data.frame(
    methodID = "Ensemble",
    AUC_mean = ifelse(is.nan(auc_mean), NA, auc_mean),
    AUC_sd   = ifelse(is.nan(auc_sd), NA, auc_sd),
    TSS_mean = ifelse(is.infinite(tss_mean) | is.nan(tss_mean), NA, tss_mean),
    TSS_sd   = ifelse(is.infinite(tss_sd) | is.nan(tss_sd), NA, tss_sd),
    Boyce_mean = ifelse(is.nan(boyce_mean), NA, boyce_mean),
    Boyce_sd   = ifelse(is.nan(boyce_sd), NA, boyce_sd)
  )
  eval_summary <- rbind(eval_summary, ensemble_row)
  message("Evaluation of the models ...")
  print(eval_summary)
  utils::write.csv(eval_summary, file.path(res_path, "Evaluation_Table.csv"), row.names = FALSE)

  ############# 11. Variable importance #############
  var_imp_raw <- suppressWarnings(sdm::getVarImp(final_model, id = 'ensemble', setting = ens_setting))
  var_data <- as.data.frame(var_imp_raw@varImportance)
  message("Variable importance:...")
  print(var_data)
  utils::write.csv(var_data, file.path(res_path, "VarImportance_Table.csv"), row.names = FALSE)

  ############## 12. Response curves #############
  if(response.curve){
    message("Generating response curves...")
    rcur_list <- suppressWarnings(sdm::rcurve(final_model, env_names, mean = TRUE, confidence = TRUE, gg = TRUE))
    rcur_data <- rcur_list$data
    p2 <- rcur_list
    utils::write.csv(rcur_data, file.path(res_path, "rcurve_Table.csv"), row.names = FALSE)
    message("Done...response curves were successfully created!\n")
  }

  ############## 13. Projection (if proj_path provided) #############
  proj_rasters <- list()
  if (!is.null(proj_path) && dir.exists(proj_path)) {
    message("Projecting to new environmental conditions...")
    proj_files <- list.files(proj_path, pattern = "\\.(asc|tif)$", full.names = TRUE)
    if (length(proj_files) > 0) {
      proj_stack <- terra::rast(proj_files)
      names(proj_stack) <- make.names(names(proj_stack))
      common_vars <- intersect(env_names, names(proj_stack))
      if (length(common_vars) != length(env_names)) {
        warning("Projection rasters have different variable names; using intersection.")
        proj_stack <- proj_stack[[common_vars]]
      }
      pred_proj <- sdm::ensemble(
        final_model,
        newdata = proj_stack,
        setting = ens_setting
      )
      names(pred_proj) <- "Habitat_suitability_SUP_proj"
      binary_proj <- pred_proj >= best_threshold
      binary_proj <- terra::app(binary_proj, fun = function(x) as.integer(x))
      binary_proj <- terra::classify(binary_proj, cbind(NA, -9999))

      names(binary_proj) <- "Presence_SUP_MAXTSS_proj"
      terra::writeRaster(pred_proj,
                         filename = file.path(res_path, "Proj_predicted_suitability.asc"),
                         overwrite = TRUE, NAflag = -9999, filetype = "AAIGrid")
      terra::writeRaster(binary_proj,
                         filename = file.path(res_path, "Proj_predicted_presence_proj.asc"),
                         overwrite = TRUE, NAflag = -9999, filetype = "AAIGrid")
      proj_rasters <- list(projection_suitability = pred_proj,
                           projection_presence = binary_proj)
    } else {
      warning("No .asc or .tif files found in proj_path")
    }
  }

  ########### return ############
  return(list(
    final_model = final_model,
    ensemble = ensemble_raster,
    evaluation = eval_summary,
    var_imp = var_data,
    res_cur = rcur_data,
    collinearity_report = collinearity_report
  ))
}
