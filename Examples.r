# Some examples of estimating quantiles 
source('CaViaR.R')

# Stationary AR(1) with gaussian innovations
set.seed(123)
y <- arima.sim(model = list(ar = 0.9),
               n = 1000,
               sd = 1)
plot(y, type = "l")

# 'SAV' model
y_quantiles <- CAViaR(y,
                      alpha = 0.05,
                      model = 'SAV')
y_quantiles
plot(y, type = "l")
lines(y_quantiles$quantile, col = "red")






