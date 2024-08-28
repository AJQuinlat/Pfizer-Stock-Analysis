# installing and importing needed libraries
library(quantmod)
library(forecast)
library(tseries)
library(nortest)
library(lmtest)
library(fDMA)
library(moments)
library(fGarch)
library(rugarch)
library(TSA)
library(FinTS)
library(psych)

# Getting Pfizer closing prices
price_pfizer <- getSymbols("PFE", src = "yahoo",  auto.assign=FALSE)
head(price_pfizer)
closing_pr <- Cl(price_pfizer)
plot(closing_pr)

# Converting Pfizer closing prices to months and inspecting trend, cycle, and irregularity
closing_pr_m <- Cl(to.monthly(price_pfizer))
dc <- decompose(as.ts(closing_pr_m))
plot(dc)

# Converting Pfizer closing prices to log returns and differencing to lag 1
par(mfrow=c(1,1))
logreturn <- diff(log(closing_pr), lag=1)
plot.ts(logreturn, xlab="Time", ylab="Log Gross Return")

# Inspecting summary statistics
describe(logreturn)

# Inspecting shape of log returns
hist(logreturn, prob=TRUE)
lines(density(na.omit(logreturn)))

# Testing for normality of log returns
lillie.test(logreturn)

# Omitting missing values
pfizer <- na.omit(logreturn)

# Testing for Stationarity
adf.test(closing_pr)
# Conclusion: Closing prices are not stationary
adf.test(pfizer)
# Conclusion: Log returns are stationary

# Inspecting the ACF and PACF of the monthly log returns
par(mfrow=c(2,1))
pfizer_acf <- acf(pfizer)
pfizer_pacf <- pacf(pfizer)
par(mfrow=c(1,1))

# Using auto.arima() function to determine best ARMA model
arma <- auto.arima(pfizer, trace=TRUE)
summary(arma)
coeftest(arma)

# Testing residuals for model adequacy
res <- residuals(arma)
res
plot(res)

# Testing for zero mean of residuals
t.test(res, mu=0)

# Testing for ARCH effects
archtest(as.ts(res),lag=NULL)

# Testing for normality of residuals
jarque.bera.test(res)
hist(res)
kurtosis(res)
skewness(res)

# Testing for independence of residuals
Box.test(res, type="Ljung")


# _________________________________________________________________

par(mfrow=c(3,1))
plot(pfizer, type="l", ylab='value of Pfizer log return',
     main='Plot of 2007-2024 monthly Pfizer log return')
plot(abs(pfizer), type="l", ylab='value of Pfizer log return',
     main='Plot of 2007-2024 monthly Pfizer absolute log return')
plot(pfizer^2, type="l", ylab='value of Pfizer log return',
     main='Plot of 2007-2024 monthly Pfizer squared log return')

acf(pfizer)
acf(abs(pfizer))
acf(pfizer^2)
par(mfrow=c(1,1))

# Getting the conditional mean model
model <- auto.arima(pfizer, trace=TRUE)
summary(model)

# Testing the adequacy of the model
a <- residuals(model)
jarque.bera.test(a) #Not Normally Distributed
hist(a)
kurtosis(a)
skewness(a)

# Testing if series has ARCH effects
ArchTest(pfizer, lags=1, demean=TRUE)
# p-value<0.0001 (There is ARCH effects)

# Inspecting the best GARCH model to use using EACF
eacf(abs(pfizer))

# ARMA(0,4) GARCH(1,3) ________________________________________________
sgarch13 <- ugarchspec(variance.model=list(model="sGARCH",
                                           garchOrder=c(1,3)), mean.model=list(armaOrder=c(0,4)),
                       distribution.model="std")
garchFit1 <- ugarchfit(spec=sgarch11, data=pfizer)
coef(garchFit)

stdres <- residuals(garchFit, standardize=T)
stdres2 <- stdres^2

# Testing for Normality in standardized residuals
jarque.bera.test(stdres)
hist(stdres)
kurtosis(stdres)
skewness(stdres)
# Not normal, but approximately normal

# Testing for Independence in standardized squared residuals
Box.test(stdres2, type="Ljung")
# p-value = 0.951 (The square standardized residuals are independent)

# Testing for ARCH effects in standardized squared residuals
archTest <- ArchTest(stdres2, lags = 1, demean = TRUE)
archTest
# p-value=0.9311 (There is no ARCH effects)


# 10-day Ahead Forecast
rff <- ugarchfit(spec=sgarch11, data=pfizer,out.sample=10)
rf <- ugarchforecast(rff, n.ahead=10, n.roll=0)
rf
plot(rf)
