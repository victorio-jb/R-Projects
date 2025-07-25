library(zoo)
library(forecast) #Needed for forecasting
library(tseries) #Needed for Jarque-Bera Test, ADF Test, Acf
library(e1071) #Needed to check for skewness
library(FinTS) #For the ARCH LM Test
library(urca)

OFWremittancedata <- read.csv('OFW Remittances Dataset.csv',h=T)$remmittance
OFWremittance <- ts(OFWremittancedata,start=c(2000,1),end=c(2024,8),frequency=12) #Creating a time series

OFWremittance1 <-ts(OFWremittancedata,start=c(2000,1),end=c(2013,12),frequency=12)
OFWremittance2 <- ts(OFWremittancedata[-c(1:length(OFWremittance1))],start=c(2014,1),end=c(2019,12),frequency=12)

################### PRELIMINARY STAGE ###########################
plot(OFWremittance2,col='blue',main='',ylab='', xlab='')
adf.test(OFWremittance2) #Stationary according to the ADF Test
summary(ur.ers(OFWremittance2)) #Not Stationary according to the DF-GLS ERS Test

Acf(as.numeric(OFWremittance2),main='ACF of OFWremittance series') #Does not decay to zero fast, implying non-stationarity
Pacf(as.numeric(OFWremittance2),main='PACF of OFWremittance series')

#Transformations for Stationarity
yt<-log(OFWremittance2)
#Seasonal Differencig + 1st Differencing of Logarithm Transformed Series
dlog.OFWremittance <- ts(diff(diff(yt,12)),start=c(2015,2),end=c(2019,12),frequency=12)
plot(dlog.OFWremittance^2,col='blue',main='',ylab='Seasonal First Diff of Log OFW Remittance', xlab='')
#Plot Shows no more evidece of trend nor seasonality
adf.test(dlog.OFWremittance)
summary(ur.ers(dlog.OFWremittance))
#Series is now stationary according to the ADF and ERS Test

############## IDENTIFICATION AND ESTIMATION ##########################

#AR Terms (1,2,10)
Pacf(as.numeric(dlog.OFWremittance),main='')

#MA Terms(1,9,10,11) include SMA Term
Acf(as.numeric(dlog.OFWremittance),main='')

#AR 1 is not significant, significant spike in PACF
xx <- arima(yt,order=c(10,1,9),seasonal=list(order=c(0,1,1)),fixed=c(NA,NA,0,0,0,0,0,0,0,NA,NA,0,0,0,0,0,0,NA,NA,NA))
summary(xx)

#All coefficients are significant
xx <- arima(yt,order=c(2,1,9),seasonal=list(order=c(0,1,1)),fixed=c(NA,NA,NA,0,0,0,0,0,0,NA,NA,NA))
summary(xx)

####################### RESIDUAL DIAGNOSTICS ########################
plot(xx$residuals)
#No significant lags for both ACF and PACF
Acf(xx$residuals,main='')
Pacf(xx$residuals,main='')
#Process is both stationary and invertible
autoplot(xx)
#Residuals are Stationart
adf.test(xx$residuals) 
summary(ur.ers(xx$residuals))

#Diagnostic Checking
#Do not Reject
Box.test(xx$residuals, lag = 24, type = "Ljung-Box", fitdf = 0) #Ljung-Box Test
#Do not Reject
ArchTest (xx$residuals)  #Testing constancy of variance
#Do not Reject
jarque.bera.test(xx$residuals) #Testing normality

#Other Visualizations for Normality
hist(xx$residuals)
skewness(xx$residuals)


################### POINT FORECASTS ############################
  
#Forecasts for the series 
expforecast <- exp(forecast(xx,60)$mean)
cbind(expforecast,exp(forecast(xx,60)$lower[,2]),exp(forecast(xx,60)$upper[,2]))

#Plotting Forecasts against actual values
OFWremittance3 <- ts(OFWremittancedata[-c(1:(length(OFWremittance1)+length(OFWremittance2)))],start=c(2020,1),end=c(2024,8),frequency=12)
plot(OFWremittance3,col='black',main='',ylab='', xlab="",xlim=c(2020,2025),ylim=c(500,5000))
lines(expforecast, col='violet')
lines(exp(forecast(xx,60)$lower[,2]),col='red') #lower 95%
lines(exp(forecast(xx,60)$upper[,2]),col='red') #upper 95%
legend('bottomright',title='Series',pt.cex=2,cex=0.8, c("Actual Values of OFW Remittances from 2020-2024", "Forecasts of OFW Remittances from 2020-2024"
                                                        ,"95% Prediction Intervals"), col=c("black", "violet","red"), lty=c(1,1),y.intersp=1 )


