library(ellipse)
library(tseries)
library(dplyr)
library(forecast)
require(fUnitRoots)

path <- "C:/Users/Valentin/Documents/Travail/ENSAE/valeurs_mensuelles.csv"
setwd(path) 
getwd()

# Loading of data
data <- as.data.frame(read.csv(path, sep = ";", header = TRUE, stringsAsFactors = FALSE))
data <- slice(data, -(1:3))
data <- rename(data, date := Libellé)
data <- rename(data, 'Indice brut de la production industrielle aéoronautique et spatiale' := !! colnames(data)[2] )
data <- data[, 1:2]
data <- arrange(data, date)
  
## Part I ##

# Question 1

# first visualisation
time_series <- ts(as.numeric(data[,2]), start=c(1990, 1), frequency=12)
plot(time_series, xlab="Date", ylab="space and aeronautics construction", main = "space and aeronautics construction")
monthplot(time_series)

# before Covid19
data <- data %>% filter(as.Date(paste0(date, "-01")) < as.Date("2020-03-01"))
time_series <- ts(as.numeric(data[,2]), start=c(1990, 1), frequency=12)
plot(time_series, xlab="Date", ylab="space and aeronautics construction", main = "space and aeronautics construction")
monthplot(time_series)
lag.plot(time_series, lags=12, layout=c(3,4), do.lines=FALSE)
decomp <- decompose(time_series)
plot(decomp)

# acf and pacf
acf(time_series)
pacf(time_series)
n <- length(time_series)

# linear model
summary(lm(time_series~seq(1,n)))

# adf and kpss test
adfTest_valid <- function(series,kmax,type){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k, " lags: residuals OK? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T) == 0) {
      noautocorr <- 1; cat("OK \n")}
    else cat("nope \n")
    k <- k + 1
  }
  return(adf)
}
adf <- adfTest_valid(time_series,24,"ct")
adf
kpss.test(time_series, null="Trend")

# Question 2

# differantiating the time series
diff_time_series <- diff(time_series, 1)
plot(diff_time_series)

# linear model
summary(lm(diff_time_series ~ seq(1, length(diff_time_series))))

# adf and kpss test
adf <- adfTest_valid(diff_time_series,24, type="nc")
adf
kpss.test(diff_time_series)

# Question 3

plot(cbind(time_series, diff_time_series))

## Part II ##

# Question 4

# acf and pacf
par(mfrow=c(1,2))
acf(diff_time_series);pacf(diff_time_series)

q_max <- 1
p_max <- 4

#LB tests for orders 1 to 24
arima401 <- arima(diff_time_series,c(4,0,1)) 
Box.test(arima401$residuals, lag=6, type="Ljung-Box", fitdf=5)

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
Qtests(arima401$residuals, 24, 5)

signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}
signif(arima401)

#test of all the possible models
arimafit <- function(estim){
  
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals,24,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),4)
  cat("coefficients nullity tests :\n")
  print(adjust)
  cat("\n tests of autocorrelation of the residuals : \n")
  print(pvals)
}

estim <- arima(diff_time_series,c(1,0,0)); arimafit(estim)
estim <- arima(diff_time_series,c(2,0,0)); arimafit(estim)
estim <- arima(diff_time_series,c(3,0,0)); arimafit(estim)
estim <- arima(diff_time_series,c(4,0,0)); arimafit(estim)
estim <- arima(diff_time_series,c(0,0,1)); arimafit(estim)
estim <- arima(diff_time_series,c(1,0,1)); arimafit(estim)
estim <- arima(diff_time_series,c(2,0,1)); arimafit(estim)
estim <- arima(diff_time_series,c(3,0,1)); arimafit(estim)

# AIC and BIC
mat <- matrix (NA, nrow=p_max+1, ncol=q_max+1)
rownames(mat) <- paste0("p=",0:p_max) 
colnames(mat) <- paste0("q=",0:q_max) 
AICs <- mat #
BICs <- mat
pqs <- expand.grid(0:p_max, 0:q_max)
for (row in 1:dim(pqs)[1]){
  p <- pqs[row, 1] 
  q <- pqs[row, 2] 
  estim <- try(arima(diff_time_series, c(p, 0, q), include.mean = F))
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim) 
}

AICs
AICs==min(AICs)
BICs 
BICs==min(BICs)

# final model
arma40 <- arima(diff_time_series, c(4, 0, 0), include.mean=F)
arma40

#adjusted R2
adj_r2 <- function(model){
  ss_res <- sum(model$residuals^2)
  p <- model$arma[1] 
  q <- model$arma[2]
  ss_tot <- sum(diff_time_series[-c(1:max(p, q))]^2) 
  n <- model$nobs-max(p, q) 
  adj_r2 <- 1-(ss_res/(n-p-q-1)) / (ss_tot/(n-1))
  return (adj_r2)
}
adj_r2(arma40)

# Question 5

arima410 <- arima(time_series, c(4, 1, 0), include.mean=F)
arima410

plot(time_series, xlab="Date" , ylab="Indice", main = "Observed vs. Predicted" )
lines(fitted(arima410), col = "red")
plot(diff_time_series, xlab="Date", ylab="Indice", main="Observed vs. Predicted" )
lines(fitted(arma40), col = "red")

## Part III ##

# Question 7

# discussion on the gaussian hypothesis
tsdiag(arma40)
jarque.bera.test(arma40$residuals)
qqnorm(arma40$residuals)
qqline(arma40$residuals, col = "red")
plot(density(arma40$residuals), xlim=c(-10,10), main="Density of residuals")
mu <- mean(arma40$residuals)
sigma <- sd(arma40$residuals)
x <- seq(-10,10)
y <- dnorm(x,mu,sigma)
lines(x, y, lwd=0.5, col="blue")

arma40$coef
phi_1 <- as.numeric(arma40$coef[1])
phi_2 <- as.numeric(arma40$coef[2])
phi_3 <- as.numeric(arma40$coef[3])
phi_4 <- as.numeric(arma40$coef[4])
sigma2 <- as.numeric(arma40$sigma)
phi_1
phi_2
phi_3
phi_4
sigma2

# checking of the roots
ar_coefs <- c(phi_1, phi_2, phi_3, phi_4)
ar_roots <- polyroot(c(1, -ar_coefs))
abs(ar_roots)
all(abs(ar_roots) > 1)

# Question 8

# prediction
XT1 = predict(arma40, n.ahead=2)$pred[1]
XT2 = predict(arma40, n.ahead=2)$pred[2]
XT1
XT2

fore = forecast(arma40, h=5, level=95)
par(mfrow=c(1,1))
plot(fore, xlim=c(2016,2021), col=1, fcol=2, shaded=TRUE, xlab="Time" , ylab="Value", main="Forecast for the time series")

# Bivariate confidence region
mean <- c(XT1, XT2)
sigma2 <- arma40$sigma2
phi_1 <- arma40$coef[1]
cov_matrix <- matrix(c(sigma2, sigma2*phi_1, sigma2*phi_1, sigma2*(1+phi_1*phi_1)), nrow=2)

alpha <- 0.05
chi2_val <- qchisq(1 - alpha, df = 2)

ellipse_points <- ellipse(cov_matrix, centre = mean, level = 1 - alpha)
plot(ellipse_points, type = 'l', xlab = 'Forecast for X_{T+1}', ylab = 'Forecast for X_{T+2}', main = '95 % bivariate confidence region', col = 'red')
points(mean[1], mean[2], pch = 19)
abline(h=XT2,v=XT1, col="blue")
abline(h=0,v=0)
grid()