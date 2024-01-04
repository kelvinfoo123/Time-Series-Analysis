library(forecast)

##### Load data 
data <- scan("wwwusage.txt", skip = 1) # The first entry is header so we skip it. 
head(data, 5)
min(data) ; max(data) ; mean(data) # Min = 83 
                                   # Max = 228 
                                   # Mean = 137.08 

##### Time series plot 
plot(data, xlab = "Time", ylab = "Number of users")
lines(data) # There is a clear upward trend at the end of the time series. The time series is not stationary. 

##### First order differencing 
## We perform first order differencing on the original data to eliminate the trend and see if the time series is more stationary. 
d1 <- diff(data)
plot(d1, xlab = "Time", ylab = "Values after first order differencing")
lines(d1) # The time series is more stationary after first order differencing. 

## ACF and PACF 
acf(d1, lag.max = 50) # Cuts off after lag 24 
pacf(d1, lag.max = 50) # Cuts off after lag 3 
# An appropriate model might be AR(3) looking at the ACF and PACF. 

## Yule Walker Estimation 
yw1 <- ar.yw(d1) 
yw1 # Yule walker suggested AR(3) model: X_t = Z_t + 1.11 X_{t-1} - 0.60 X_{t-2} + 0.30 X_{t-3}

## Fit ARIMA(3, 1, 0) model 
arima1 <- arima(data, order = c(3, 1, 0))
arima1
tsdiag(arima1) # The fitted ARIMA(3, 1, 0) model is adequate because the residual terms look random, making it white noise, the ACF cuts off after lag 0 and the p-value indicates significance. 
AIC(arima1) # 511.994 

##### Second order differencing 
d2 <- diff(d1)
plot(d2, xlab = "Time", ylab = "Values after second order differencing")
lines(d2)

## ACF and PACF 
acf(d2, lag.max = 50) # Cuts off after lag 27 
pacf(d2, lag.max = 50) # Cuts off after lag 2 
# An appropriate model might be AR(2) after looking at the ACF and PACF. 

## Yule Walker Estimation
yw2 <- ar.yw(d2)
yw2 # Yule Walker suggested AR(2) mode: X_t = Z_t + 0.25 X_{t-1} - 0.43 X_{t-2}

arima2 <- arima(data, order = c(2, 2, 0))
arima2
tsdiag(arima2) # The fitted ARIMA(2, 2, 0) model is adequate because the residual terms look random, making it white noise, the ACF cuts off after lag 0 and the p-value indicates significance. 
AIC(arima2) # 511.465 

##### Auto ARIMA 
autoarima <- auto.arima(data) # autoarima suggested ARMA(1, 1, 1) model 
autoarima 
arima3 <- arima(data, order = c(1, 1, 1))
arima3 
tsdiag(arima3) # The fitted ARIMA(1, 1, 1) model is adequate because the residual terms look random, making it white noise, the ACF cuts off after lag 0 and the p-value indicates significance. 
AIC(arima3) # 514.300

# The best model in terms of AIC is ARIMA(2, 2, 0). 