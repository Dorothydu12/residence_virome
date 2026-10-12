library(svd)
library(forecast)
####
calculate_errors <- function(actual, forecast) {
  mae <- mean(abs(actual - forecast))
  mse <- mean((actual - forecast)^2)
  rmse <- sqrt(mse)
  r2 <- 1 - sum((actual - forecast)^2) / sum((actual - mean(actual))^2)
  return(list(MAE = mae, MSE = mse, RMSE = rmse, R2 = r2))
}

run_nnetar <- function(signal_train, signal_test, env_train, env_test, test_size) {
  nnetar_model <- nnetar(ts(signal_train), xreg = env_train, seasonal=FALSE)
  predicted_signal_nnetar <- predict(nnetar_model, h = test_size, xreg = env_test)
  errors_nnetar <- calculate_errors(signal_test, predicted_signal_nnetar$mean)
  
  return(list(predict = predicted_signal_nnetar, errors = errors_nnetar))
}



vOTU <- read.csv("/Users/shicong/Desktop/Residence-3#/DL_model/individual_specific/R3_Leftpalm_vOTU.csv",header = T,row.names = 1)
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

ggplot(data.frame(Component = 1:length(explained_variance),
                  CumulativeVariance = cumulative_variance),
       aes(x = Component, y = CumulativeVariance)) +
  geom_line(size = 1) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "red", size = 1) +
  theme_bw()

ggsave("location/R4/virus_top_5_signals_handrail_subway.pdf", height = 5.04, width = 7.35)


U1 <- svd_result_vOTU$u
D1 <- diag(svd_result_vOTU$d)
V1 <- svd_result_vOTU$v
k = 3
U1_k <- U1[, 1:k]
###training 80%sample

train_size = floor(nrow(U1)*0.8)
test_size = nrow(U1)-train_size
U1_k_train <- U1_k[1:train_size, ]
U1_k_test <- U1_k[-(1:train_size), ]



###rMAG
rMAG <-read.csv("/Users/shicong/Desktop/Residence-3#/DL_model/individual_specific/R4_Leftpalm_rMAG.csv",header = T,row.names = 1)
rMAG_numeric <- rMAG %>%
  mutate(across(everything(), ~ scales::rescale(as.numeric(.), to = c(0, 1))))
data_matrix_rMAG <- as.matrix(rMAG_numeric)
svd_result_rMAG <- svd(data_matrix_rMAG)

explained_variance1 <- svd_result_rMAG$d^2 / sum(svd_result_rMAG$d^2)
cumulative_variance1 <- cumsum(explained_variance1)
threshold <- 0.80
optimal_j_cumulative <- which(cumulative_variance1 >= threshold)[1]
cat("Bacteria::Optimal j using cumulative variance threshold (80%):", optimal_j_cumulative, "\n")
write.csv(explained_variance1,"location/R4/bacteria_handrail_subway_explained_variance.csv")

ggplot(data.frame(Component = 1:length(explained_variance),
                  CumulativeVariance = cumulative_variance),
       aes(x = Component, y = CumulativeVariance)) +
  geom_line(size = 1) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "red", size = 1) +
  theme_bw()

ggsave("location/R4/bacteria_top_2_signals_handrail_subway.pdf", height = 5.04, width = 7.35)

U <- svd_result_rMAG$u
D <- diag(svd_result_rMAG$d)
V <- svd_result_rMAG$v
j = 3
U_j <- U[, 1:j]

U_j_train <- U_j[1:train_size, ]
U_j_test <- U_j[-(1:train_size), ]

metrics_nnetar <- list()
plot_data <- list()

for (i in 1:j) {
  signal_train <- U_j_train[, i]
  signal_test <- U_j_test[, i]
  nnetar_result <- run_nnetar(signal_train,signal_test, U1_k_train, U1_k_test, test_size)

  
  metrics_nnetar[[i]] <- data.frame(Signal = paste("Signal", i), Model = "NNETAR", MAE = nnetar_result$errors$MAE, RMSE = nnetar_result$errors$RMSE, R2 = nnetar_result$errors$R2)

  
  plot_data[[i]] <- data.frame(
    Time = c(1:train_size, (train_size + 1):(train_size + test_size)),
    Value = c(signal_train, signal_test),
    Model = "Real",
    Signal = paste("Signal", i)
  )
  
  plot_data[[i]] <- rbind(
    plot_data[[i]],
    data.frame(Time = (train_size + 1):(train_size + test_size),Value = nnetar_result$predict$mean, Model = "NNETAR", Signal = paste("Signal", i))
  )
}

plot_data_combined <- do.call(rbind, plot_data)
metrics_combined <- rbind(do.call(rbind, metrics_nnetar))
print(metrics_combined)

write.csv(plot_data_combined, "/Users/shicong/Desktop/Residence-3#/DL_model/individual_specific/obs_pred_r3_r4.csv", row.names = FALSE)
write.csv(metrics_combined, "/Users/shicong/Desktop/Residence-3#/DL_model/individual_specific/model_parameter_r3_r4.csv", row.names = FALSE)

plot_data_combined$Real_Color <- ifelse(plot_data_combined$Model == "Real", "black","red")

ggplot(plot_data_combined, aes(x = Time, y = Value, color = ifelse(Model == "Real", Real_Color, Model), linetype = Model)) +
  geom_line(size = 1, alpha = 0.9) +
  scale_color_manual(values = c("black" = "black", "red" = "red", "NNETAR" = "red")) +
  scale_linetype_manual(values = c("Real" = "solid", "ARIMA" = "solid", "NNETAR" = "solid", "Prophet" = "solid", "Neural Network" = "solid")) +
  facet_wrap(~ Signal, scales = "free_y", ncol = 3) +theme_bw() +theme(legend.position = "none")

ggsave("/Users/shicong/Desktop/Residence-3#/DL_model/individual_specific/r3_r4.pdf", height = 3.66, width = 15.8)

explained_variance1

##plot
df=read.csv("/Users/shicong/Desktop/Residence-3#/DL_model/location/obs_pred1.csv",header = T)
df$site=factor(df$site,levels=c("left palm","right palm","door knob","bed headboard","park/campus handrail","subway_exit handrail"))
ggplot(df,aes(y=Real,x=NNETAR))+geom_point(size=2,alpha=0.5)+geom_smooth(method = lm,se=FALSE)+geom_abline(intercept = 0, slope = 1,colour = "red", linetype = "dashed")+facet_wrap(~site)+ stat_regline_equation(aes(label = ..rr.label..))+theme_bw()+theme(legend.title= element_blank())+ylab("Real value") +xlab("Predicted value")

