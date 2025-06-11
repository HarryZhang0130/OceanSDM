#' Title This function was used to fit temperature-bin data from satellite
#'     studies on oceanic sharks
#'
#' @param ws the temperature-bin data with three columns: freq,lower,upper.
#'     freq, the frequency of occurrences; lower, the lower bound of the
#'     temperature bin; upper, the upper bound of the temperature bin.
#'     See an example data "ws.RData" in the folder "data".
#' @param factor the conversion factor that was used to transform original
#'     temperature to a value between 0 and 1.
#'
#' @return This function returns the model fit, AIC values, the estimated
#'     preferred and tolerable temperatures based on the best fitted model，
#'     and a ggplot of the model fits.
#' @export
#'
#' @examples
#' data(ws)
#' fit_temp(ws,30)
fit_temp <- function(ws,factor) {
  # 加载必要的包
  library(fitdistrplus)
  library(tidyverse)
  # 创建模拟数据
  set.seed(123)  # 确保结果可重现
  for(i in 1:nrow(ws)){
    cl<-c(runif(ws$freq[i], ws$lower[i], ws$upper[i]))
    if(exists("ws_cl")){
      ws_cl<-c(ws_cl,cl)
    }else{
      ws_cl<-cl
    }
  }
  for(i in 1:nrow(ws)){
    left0<-rep(ws$lower[i],ws$freq[i])
    if(exists("left")){
      left<-c(left,left0)
    }else{
      left<-left0
    }
  }
  for(i in 1:nrow(ws)){
    right0<-rep(ws$upper[i],ws$freq[i])
    if(exists("right")){
      right<-c(right,right0)
    }else{
      right<-right0
    }
  }
  ws1<-cbind(ws_cl,left,right)
  ws2<-as.data.frame(ws1)
  # Create data frame with just "left" and "right" columns, one row per respondent, ready for fitdistcens
  ws2_cl <- dplyr::select(ws2, left, right) %>%
    as.data.frame()
  # 归一化数据 (0-1范围)
  ws_df_scaled <- ws2_cl/factor
  # 拟合分布
  fits <- list(
    gamma = fitdistcens(ws_df_scaled, "gamma"),
    weibull = fitdistcens(ws_df_scaled, "weibull"),
    normal = fitdistcens(ws_df_scaled, "norm"),
    beta = tryCatch(fitdistcens(ws_df_scaled, "beta"), error = function(e) NULL)
  )

  # 移除拟合失败的分布
  fits <- fits[!sapply(fits, is.null)]

  # 比较AIC并选择最佳模型
  aic_values <- sapply(fits, function(f) f$aic)
  best_fit_name <- names(which.min(aic_values))
  best_fit <- fits[[best_fit_name]]

  # 计算分位数
  quantiles <- c(0.025, 0.10, 0.90, 0.975)
  q_values <- quantile(best_fit, quantiles)

  # 还原到原始尺度
  scaled_quantiles <- q_values$quantiles * factor
  names(scaled_quantiles) <- c("tolerance_lower", "preference_lower",
                               "preference_upper", "tolerance_upper")

  #
  library(ggplot2)
  p<-ggplot(ws2) +
    geom_density(aes(x = ws_cl/factor)) +
    stat_function(fun = dgamma, args = fits$gamma$estimate, colour = "red")+
    stat_function(fun = dbeta, args = fits$beta$estimate, colour = "blue")+
    stat_function(fun = dweibull, args = fits$weibull$estimate, colour = "green")+
    stat_function(fun = dnorm, args = fits$normal$estimate, colour = "purple")+
    xlab("Ocean temperature")+
    ylab("Density")
  # 返回结果
  list(
    best_model = best_fit_name,
    aic_values = aic_values,
    preference_range = scaled_quantiles[c("preference_lower", "preference_upper")],
    tolerance_range = scaled_quantiles[c("tolerance_lower", "tolerance_upper")],
    all_fits = fits,
    quantiles = scaled_quantiles,
    plots = p
  )
}
