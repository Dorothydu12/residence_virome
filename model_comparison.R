
calculate_errors <- function(actual, forecast) {
  mae <- mean(abs(actual - forecast))
  mse <- mean((actual - forecast)^2)
  rmse <- sqrt(mse)
  r2 <- 1 - sum((actual - forecast)^2) / sum((actual - mean(actual))^2)
  return(list(MAE = mae, MSE = mse, RMSE = rmse, R2 = r2))
}

# LM with Log Transformation Model - 矩阵/array专用版
run_lm_log <- function(signal_train, signal_test, env_train, env_test, test_size) {
  env_train <- as.matrix(env_train)
  env_test <- as.matrix(env_test)
  if (ncol(env_train) != ncol(env_test)) {
    stop("env_train 和 env_test 的列数不匹配")
  }
  log_signal_train <- log(signal_train + 1)
  X_train <- cbind(1, env_train)
  X_test <- cbind(1, env_test)
  lm_fit <- lm.fit(x = X_train, y = log_signal_train)
  forecasted_log_signal <- X_test %*% lm_fit$coefficients
  forecasted_signal_lm_log <- exp(forecasted_log_signal) - 1
  if (is.matrix(forecasted_signal_lm_log)) {
    forecasted_signal_lm_log <- as.vector(forecasted_signal_lm_log)
  }
  errors_lm_log <- calculate_errors(signal_test, forecasted_signal_lm_log)
  return(list(forecast = forecasted_signal_lm_log, errors = errors_lm_log))
}

# ARIMA Model
run_arima <- function(signal_train, signal_test, env_train, env_test, test_size) {
  arima_model <- auto.arima(signal_train, xreg = env_train)
  forecasted_signal_arima <- forecast(arima_model, xreg = env_test, h = test_size)
  predicted_values <- as.numeric(forecasted_signal_arima$mean)
  errors_arima <- calculate_errors(signal_test, predicted_values)
  return(list(forecast = forecasted_signal_arima, errors = errors_arima))
}

# NNETAR Model
#run_nnetar <- function(signal_train, signal_test, env_train, env_test, test_size) {
  nnetar_model <- nnetar(ts(signal_train), xreg = env_train, seasonal=FALSE)
  predicted_signal_nnetar <- predict(nnetar_model, h = test_size, xreg = env_test)
  errors_nnetar <- calculate_errors(signal_test, predicted_signal_nnetar$mean)
  return(list(predict = predicted_signal_nnetar, errors = errors_nnetar))
}
run_nnetar <- function(signal_train, signal_test, env_train, env_test, test_size) {
  
  # 训练模型
  nnetar_model <- forecast::nnetar(
    y = signal_train, 
    xreg = env_train, 
    seasonal = FALSE
  )
  
  # 预测 - predict.nnetar 可能返回列表或向量，需要处理
  predicted_signal_nnetar <- predict(nnetar_model, h = test_size, xreg = env_test)
  
  # 【关键】确保提取的是数值向量，不是列表
  if (is.list(predicted_signal_nnetar)) {
    predicted_values <- as.numeric(predicted_signal_nnetar$mean)  # 或 pred_result$fitted
  } else {
    predicted_values <- as.numeric(predicted_signal_nnetar)
  }
  
  # 确保长度正确
  predicted_values <- predicted_values[1:test_size]
  actual_values <- as.numeric(signal_test)[1:test_size]
  
  # 计算误差
  errors_nnetar <- calculate_errors(actual_values, predicted_values)
  
  # 返回统一结构
  return(list(predict = predicted_signal_nnetar, errors = errors_nnetar))
}

vOTU <- read.csv("/Users/shicong/Desktop/Residence-3#/temproal_prediction/R1/Doorknob_vOTU.csv",header = T,row.names = 1)
vOTU$Day_1 <- as.Date(vOTU$Day_1, format = "%d/%m/%Y")
vOTU_numeric <- vOTU %>%
  dplyr::select(-Day_1) %>%
  mutate(across(everything(), ~ scales::rescale(as.numeric(.), to = c(0, 1))))
data_matrix_vOTU <- as.matrix(vOTU_numeric)
svd_result_vOTU <- svd(data_matrix_vOTU)

explained_variance <- svd_result_vOTU$d^2 / sum(svd_result_vOTU$d^2)
cumulative_variance <- cumsum(explained_variance)
threshold <- 0.80
optimal_k_cumulative <- which(cumulative_variance >= threshold)[1]
cat("Virus: Optimal k using cumulative variance threshold (80%):", optimal_k_cumulative, "\n")
#write.csv(explained_variance,"location/R2/virus_handrail_park_explained_variance.csv")
U1 <- svd_result_vOTU$u
D1 <- diag(svd_result_vOTU$d)
V1 <- svd_result_vOTU$v
k = 2
U1_k <- U1[, 1:k]
###training 80%sample

train_size = floor(nrow(U1)*0.8)
test_size = nrow(U1)-train_size
U1_k_train <- U1_k[1:train_size, ]
U1_k_test <- U1_k[-(1:train_size), ]
train_size
test_size
# Prophet Model
run_prophet <- function(signal_train, signal_test, dates, env_train = NULL, env_test = NULL) {
  
  # 日期分割
  date_train <- dates[1:13]
  date_test  <- dates[14:17]
  
  # 构建数据框
  df_train <- data.frame(ds = date_train, y = as.vector(signal_train))
  df_test  <- data.frame(ds = date_test)
  
  # 添加外生变量
  if (!is.null(env_train)) {
    env_train <- as.data.frame(env_train)
    colnames(env_train) <- paste0("xreg_", 1:ncol(env_train))
    df_train <- cbind(df_train, env_train)
    
    env_test <- as.data.frame(env_test)
    colnames(env_test) <- colnames(env_train)
    df_test <- cbind(df_test, env_test)
  }
  
  # ===== 关键修复：小样本参数调整 =====
  m <- prophet(
    # 趋势设置（针对10个观测值优化）
    growth = 'linear',           # 线性趋势（小样本更稳定）
    n.changepoints = 3,          # 显式设置变点数量（< 10）
    changepoint.range = 0.8,     # 变点只在数据前80%范围内
    
    # 关闭所有季节性（数据量不足）
    yearly.seasonality = FALSE,
    weekly.seasonality = FALSE,
    daily.seasonality = FALSE,
    
    # 减少先验强度（避免过拟合）
    changepoint.prior.scale = 0.05,  # 趋势变化更保守
    interval.width = 0.95            # 预测区间
  )
  
  # 添加回归量
  if (!is.null(env_train)) {
    for (col in colnames(env_train)) {
      m <- add_regressor(m, name = col, prior.scale = 0.1)  # 弱先验
    }
  }
  
  # 拟合与预测
  m <- fit.prophet(m, df_train)
  fc <- predict(m, df_test)
  
  return(list(
    predict = fc$yhat,
    errors = calculate_errors(signal_test, fc$yhat),
    model = m
  ))
}

###rMAG
rMAG <-read.csv("/Users/shicong/Desktop/Residence-3#/temproal_prediction/R1/Doorknob_rMAG.csv",header = T,row.names = 1)
rMAG_numeric <- rMAG %>%
  mutate(across(everything(), ~ scales::rescale(as.numeric(.), to = c(0, 1))))
data_matrix_rMAG <- as.matrix(rMAG_numeric)
svd_result_rMAG <- svd(data_matrix_rMAG)

explained_variance1 <- svd_result_rMAG$d^2 / sum(svd_result_rMAG$d^2)
cumulative_variance1 <- cumsum(explained_variance1)
threshold <- 0.80
optimal_j_cumulative <- which(cumulative_variance1 >= threshold)[1]
cat("Bacteria::Optimal j using cumulative variance threshold (80%):", optimal_j_cumulative, "\n")

U <- svd_result_rMAG$u
D <- diag(svd_result_rMAG$d)
V <- svd_result_rMAG$v
j = 2
U_j <- U[, 1:j]

U_j_train <- U_j[1:train_size, ]
U_j_test <- U_j[-(1:train_size), ]


metrics_lm_log <- list()
metrics_arima <- list()
metrics_prophet <- list()
metrics_nnetar <- list()
plot_data <- list()
# Loop over each signal to collect plot data
for (i in 1:j) {
  signal_train <- U_j_train[, i]
  signal_test <- U_j_test[, i]
  signal_train <- ts(signal_train)
  signal_test <- ts(signal_test)
  lm_log_result <- run_lm_log(signal_train, signal_test, U1_k_train, U1_k_test, test_size)
  arima_result <- run_arima(signal_train, signal_test, U1_k_train, U1_k_test, test_size)
  prophet_result <- run_prophet(signal_train, signal_test, vOTU$Day_1,U1_k_train, U1_k_test)
  nnetar_result <- run_nnetar(signal_train, signal_test, U1_k_train, U1_k_test,test_size)
  # Store the metrics
  metrics_prophet[[i]] <- data.frame(Signal = paste("Signal", i), Model = "Prophet", MAE = prophet_result$errors$MAE, RMSE = prophet_result$errors$RMSE, R2 = prophet_result$errors$R2)
  metrics_arima[[i]] <- data.frame(Signal = paste("Signal", i), Model = "ARIMA", MAE = arima_result$errors$MAE, RMSE = arima_result$errors$RMSE, R2 = arima_result$errors$R2)
  metrics_nnetar[[i]] <- data.frame(Signal = paste("Signal", i), Model = "NNETAR", MAE = nnetar_result$errors$MAE, RMSE = nnetar_result$errors$RMSE, R2 = nnetar_result$errors$R2)
  metrics_lm_log[[i]] <- data.frame(Signal = paste("Signal", i), Model = "lm_log", MAE = lm_log_result$errors$MAE,  RMSE = lm_log_result$errors$RMSE, R2 = lm_log_result$errors$R2)
  
  # Combine data for plotting
  plot_data[[i]] <- data.frame(
    Time = c(1:train_size, (train_size + 1):(train_size + test_size)),
    Value = c(signal_train, signal_test),
    Model = "Real",
    Signal = paste("Signal", i)
  )
  plot_data[[i]] <- rbind(
    plot_data[[i]],
    data.frame(Time = (train_size + 1):(train_size + test_size), 
               Value = prophet_result$predict, 
               Model = "Prophet", 
               Signal = paste("Signal", i)),
    data.frame(Time = (train_size + 1):(train_size + test_size), 
               Value = as.numeric(arima_result$forecast$mean),  # convert ts to numeric
               Model = "ARIMA", 
               Signal = paste("Signal", i)),
    data.frame(Time = (train_size + 1):(train_size + test_size), 
               Value = as.numeric(nnetar_result$predict$mean),  # extract mean from forecast object
               Model = "NNETAR", 
               Signal = paste("Signal", i)),
    data.frame(Time = (train_size + 1):(train_size + test_size), 
               Value = lm_log_result$forecast, 
               Model = "lm_log", 
               Signal = paste("Signal", i))
  )
  #plot_data[[i]] <- rbind(
    #plot_data[[i]],
    #data.frame(Time = (train_size + 1):(train_size + test_size), Value = prophet_result$predict, Model = "Prophet", Signal = paste("Signal", i)),
    #data.frame(Time = (train_size + 1):(train_size + test_size), Value = arima_result$forecast$mean, Model = "ARIMA", Signal = paste("Signal", i)),
    #data.frame(Time = (train_size + 1):(train_size + test_size), Value = nnetar_result$predict, Model = "NNETAR", Signal = paste("Signal", i)),
    #data.frame(Time = (train_size + 1):(train_size + test_size), Value = lm_log_result$forecast, Model = "lm_log", Signal = paste("Signal", i))
  #)
}

# Combine all plot data
plot_data_combined <- do.call(rbind, plot_data)
write.csv(plot_data_combined,"../revision/model_comparasion/plot_data_combined_Doorknob_R1.csv")

metrics_lm_log_results <- do.call(rbind, metrics_lm_log)
metrics_arima_results <- do.call(rbind, metrics_arima)
metrics_prophet_results <- do.call(rbind, metrics_prophet)
metrics_nnetar_results <- do.call(rbind, metrics_nnetar)

comparison_metrics <- rbind(metrics_nnetar_results,metrics_lm_log_results,metrics_arima_results, metrics_prophet_results)
# Display the comparison metrics
print(comparison_metrics)
write.csv(comparison_metrics,"../revision/model_comparasion/parameter_Doorknob_R1.csv")

# Plot the data for all signals with actual and predicted values from ARIMA, Prophet, and NNETAR
ggplot(plot_data_combined, aes(x = Time, y = Value, color = Model)) +
  geom_line() +
  facet_wrap(~ Signal, scales = "free_y", ncol = 2) +  # Adjust columns based on how you want to arrange the plots
  labs(title = "subway exit handrail in Residece 4",
       x = "Time", y = "Value") +
  theme_classic() +
  scale_color_manual(values = c("ARIMA" = "blue", "Prophet" = "green", "NNETAR" = "red", "Actual" = "black","lm_log"= "purple"))+scale_x_continuous(breaks = seq(0, 19, by = 1))



