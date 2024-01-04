library(lubridate)
library(forecast)
library(ggplot2)

##### Basic exploration of data 
data <- read.delim("drug.txt", header = TRUE, sep = ",")
head(data, 5)
begin <- ymd(as.Date(data$date[1])) ; begin # Data begins at 1997-07-01
end <- ymd(as.Date(data$date[nrow(data)])) ; end # Data ends at 2008-06-01
drug <- ts(data$value, start = c(year(begin), month(begin)), 
                       end = c(year(end), month(end)), 
                       frequency = 12) 
plot(drug, xlab = "Year", ylab = "Drug Sales") # Upward trend and seasonal pattern in time series 
par(mfrow = c(1,2))
acf(data$value) # ACF dies down slowly and thus, the data is non-stationary. 
pacf(data$value)   
par(mfrow = c(1,1))
stl_result <- stl(drug, s.window = "periodic")
autoplot(stl_result) # Increasing trend with seasonality in time series. The remainder is close to white noise. 

##### Split data to train test set 
train_size <- round(0.8 * length(drug)) ; train_size #163 training observations 
train <- head(drug, train_size)
test <- tail(drug, round(length(drug) - train_size)) ; length(test) # 41 test observations 

# Making training data stationary 
BoxCox.lambda(train) # 0.26. Use lambda = 0. 
train_1 <- BoxCox(train, 0) # Apply Box-cox transformation with lambda = 0 to stabilize variance
train_2 <- diff(train_1, differences = 1) # One time differencing to remove trend 
autoplot(train_2)
train_3 <- diff(train_2, lag = 12) # Differencing with lag 12 to remove seasonality 
autoplot(train_3) # Time series is now pure random noise fluctuation 
# For SARIMA model, use d = 1 and D = 1. 
par(mfrow=c(1,2))
acf(train_3, lag = 50) ; pacf(train_3, lag = 50) # p = 5 and q = 4 
par(mfrow = c(1, 1))

##### Model building for SARIMA and forecasting 
arima1 <- arima(train_1, 
                order = c(5, 1, 4), 
                seasonal = list(order = c(0,1,0), period = 12))
checkresiduals(arima1)
tsdiag(arima1) # The model is adequate 
pred.arima1 <- InvBoxCox(forecast(arima1, h = length(test))$mean, lambda = 0)
autoplot(pred.arima1)
accuracy(pred.arima1, test) # RMSE = 2.919 and MAPE = 13.07 

autoarima <- auto.arima(train_1)
autoarima # Proposed model is ARIMA(0,1,1)(2,1,1)[12]
checkresiduals(autoarima) # Residuals resemble white noise 
tsdiag(autoarima) # p-value larger than 0.005. The model is adequate. 
pred.autoarima <- InvBoxCox(forecast(autoarima, h = length(test))$mean, lambda = 0) 
accuracy(pred.autoarima, test) # RMSE = 1.66 and MAPE = 7.61 

##### Forecasting with Holt-Winter model 
hw_additive <- hw(train_1, seasonal= "additive")
hw_multi <- hw(train_1, seasonal = "multiplicative")
checkresiduals(hw_additive)
checkresiduals(hw_multi) # Residual for both resemble white noise 

pred.hw_add <- InvBoxCox(forecast(hw_additive$mean, h = length(test))$mean, lambda = 0) 
accuracy(pred.hw_add, test) # RMSE = 1.86 and MAPE = 7.47 
pred.hw_multi <- InvBoxCox(forecast(hw_multi$mean, h = length(test))$mean, lambda = 0)
accuracy(pred.hw_multi, test) # RMSE = 3.25 and MAPE = 12.61 

#### Based on RMSE and MAPE, the best model is ARIMA(0,1,1)(2,1,1)[12]. 
pred.autoarima.forecast <- InvBoxCox(forecast(autoarima, h = length(test) + 20)$mean, lambda = 0) 
autoplot(pred.autoarima.forecast, ylab = "Forecasted Drug Sales")