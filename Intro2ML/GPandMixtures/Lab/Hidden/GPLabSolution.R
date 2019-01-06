################################################
# Problem 1 - Gaussian Process Regression
################################################

# 1 a)

# Setting up the squared exponential kernel
SquaredExpKernel <- function(x1,x2,sigmaF=1,l=1){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}

# Defining the PosteriorGP function that computes the posterior mean and variance
PosteriorGP <- function(x, y, K=SquaredExpKernel, hyperParam, sigmaNoise, xStar=x){
  
  k_x_x <- K(x,x,hyperParam[1],hyperParam[2])
  k_x_xStar <- K(x,xStar,hyperParam[1],hyperParam[2])
  k_xStar_x <- K(xStar,x,hyperParam[1],hyperParam[2])
  k_xStar_xStar <- K(xStar,xStar,hyperParam[1],hyperParam[2])
  nx <- ncol(k_x_x)

  # Calculate the mean and covariance functions
  fPostMean <- k_xStar_x%*%solve(k_x_x + sigmaNoise^2*diag(1,nx))%*%y
  fPostCov <- k_xStar_xStar - k_xStar_x%*%solve(k_x_x + sigmaNoise^2*diag(1,nx))%*%k_x_xStar

  return(list(fPostMean = fPostMean, fPostCov = fPostCov))
}


# 1 b)

# Data
x <- c(-1.0,-0.6,-0.2,0.4,0.8)
y <- c(0.768,-0.044,-0.940,0.719,-0.664)

# Only first observation
xStar <- seq(-1,1,by =0.01)
GPPost <- PosteriorGP(x[4], y[4], K=SquaredExpKernel, hyperParam=c(1,0.3), sigmaNoise = 0.1, xStar)
plot(x,y, ylim = c(-3,3))
lines(xStar,GPPost$fPostMean,type="l")
lines(xStar,GPPost$fPostMean - 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")
lines(xStar,GPPost$fPostMean + 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")


# 1 c)

# Only two first observations
xStar <- seq(-1,1,by =.1)
GPPost <- PosteriorGP(x[c(2,4)], y[c(2,4)], K=SquaredExpKernel, hyperParam=c(1,0.3), sigmaNoise = 0.1, xStar)
plot(x,y, ylim = c(-3,3))
lines(xStar,GPPost$fPostMean,type="l")
lines(xStar,GPPost$fPostMean - 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")
lines(xStar,GPPost$fPostMean + 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")

# 1 d)
# All 5 observations
xStar <- seq(-1,1,by =.1)
GPPost <- PosteriorGP(x, y, K=SquaredExpKernel, hyperParam=c(1,0.3), sigmaNoise = 0.1, xStar)
plot(x,y, ylim = c(-3,3))
lines(xStar,GPPost$fPostMean,type="l")
lines(xStar,GPPost$fPostMean - 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")
lines(xStar,GPPost$fPostMean + 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")


# 1 e)

# Repeat with lengthScale = 1

# Only first observation
xStar <- seq(-1,1,by =.1)
GPPost <- PosteriorGP(x[4], y[4], K=SquaredExpKernel, hyperParam=c(1,1), sigmaNoise = 0.1, xStar)
plot(x,y, ylim = c(-3,3))
lines(xStar,GPPost$fPostMean,type="l")
lines(xStar,GPPost$fPostMean - 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")
lines(xStar,GPPost$fPostMean + 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")

# Only two first observations
xStar <- seq(-1,1,by =.1)
GPPost <- PosteriorGP(x[c(2,4)], y[c(2,4)], K=SquaredExpKernel, hyperParam=c(1,1), sigmaNoise = 0.1, xStar)
plot(x,y, ylim = c(-3,3))
lines(xStar,GPPost$fPostMean,type="l")
lines(xStar,GPPost$fPostMean - 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")
lines(xStar,GPPost$fPostMean + 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")

# All 5 observations
xStar <- seq(-1,1,by =.1)
GPPost <- PosteriorGP(x, y, K=SquaredExpKernel, hyperParam=c(1,1), sigmaNoise = 0.1, xStar)
plot(x,y, ylim = c(-3,3))
lines(xStar,GPPost$fPostMean,type="l")
lines(xStar,GPPost$fPostMean - 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")
lines(xStar,GPPost$fPostMean + 1.96*sqrt(diag(GPPost$fPostCov)),type="l", col="red")

# 2) Japan temp data
# read data
japanTemp <- read.table('~/Dropbox/Teaching/IntroToML/GPandMixtures/Lab/JapanTemp.dat', header = TRUE)
y <- japanTemp$temp
x <- japanTemp$time
xStar <- x
lmModel <- lm(y ~ x)
GPPost <- PosteriorGP(x, y, K=SquaredExpKernel, hyperParam=c(100,0.2), sigmaNoise = sd(lmModel$residuals), xStar)
plot(x,y, cex=0.4)
lines(xStar,GPPost$fPostMean, type="l", col="blue", lwd = 3)
lines(xStar,GPPost$fPostMean - 1.96*sqrt(diag(GPPost$fPostCov)),type="l", lwd = 2, col="red")
lines(xStar,GPPost$fPostMean + 1.96*sqrt(diag(GPPost$fPostCov)),type="l", lwd = 2, col="red")


################################################
# Problem 2 - EM for mixture models 
################################################

# 2 a)

mixtureEM <- function(data, K, initMu, initSigma, initPi, tol){
  
  # Preliminaries
  count <- 0
  nTotal <- length(data)
  N <- rep(0,K)
  gamma = matrix(0,nTotal,K)
  Mu = initMu
  Sigma = initSigma
  Pi = initPi
  
  LogLOld <- 10^10
  LogLDiff <- 10^10
  while (LogLDiff > tol){
    count <- count + 1
    # E-step
  
    for (k in 1:K){
      gamma[,k] = Pi[k]*dnorm(data, mean=Mu[k], sd = sqrt(Sigma[k]))
    }
    sumGamma <- rowSums(gamma)
    for (k in 1:K){
      gamma[,k] = gamma[,k]/sumGamma
    }
 
    # M-step
    for (k in 1:K){
      N[k] <- sum(gamma[,k])
      Mu[k] = (1/N[k])*sum(gamma[,k]*data)
      Sigma[k] = sqrt((1/N[k])*sum(gamma[,k]*(data-Mu[k])^2))
      Pi[k] = N[k]/nTotal
    }
   
    # Log-Likelihood computation
    for (k in 1:K){
      gamma[,k] = Pi[k]*dnorm(data, mean=Mu[k], sd = sqrt(Sigma[k]))
    }
    LogL = sum(log(rowSums(gamma)))
    print(c(LogL,LogLOld))
    LogLDiff = abs(LogL - LogLOld)
    LogLOld = LogL
    
  }
  return(list(Mu = Mu,Sigma = Sigma, Pi = Pi, LogL = LogL, nIter = count))
}

# Generate some data for testing
data <- c(rnorm(50,1,1),rnorm(50,5,0.5))
initMu = c(-1,1)
initSigma = c(1,1)
initPi = c(0.5,0.5)

# Run the EM
K <- 2
EMfit <- mixtureEM(data, K, initMu, initSigma, initPi, tol = 0.0000001)
EMfit

# Plot the histogram and density fit
dataGrid <- seq(min(data)-0.2*IQR(data),max(data)+0.2*IQR(data), length = 1000)
densEst <- rep(0,length(dataGrid))
for (k in 1:K){
  densEst <- densEst + EMfit$Pi[k]*dnorm(dataGrid, mean=EMfit$Mu[k], sd = sqrt(EMfit$Sigma[k]))
}
hist(data, 30, freq = FALSE, xlim = c(min(data)-0.2*IQR(data),max(data)+0.2*IQR(data)))
lines(dataGrid,densEst, col="red", lwd = 2)

# Checking the correctness of the code by comparing to widely used package
#install.packages("mixtools")
library(mixtools)
EMfitMixTools <- normalmixEM(data, lambda = initPi[1], mu = initMu, sigma = initSigma)

EMfitMixTools$lambda
EMfit$Pi

EMfitMixTools$mu
EMfit$Mu

EMfitMixTools$sigma
EMfit$Sigma

EMfit

# Plotting also the fit from MixTools
densEst <- rep(0,length(dataGrid))
for (k in 1:K){
  densEst <- densEst + EMfitMixTools$lambda[k]*dnorm(dataGrid, mean=EMfitMixTools$mu[k], sd = sqrt(EMfitMixTools$sigma[k]))
}
lines(dataGrid,densEst, col="blue", lwd = 2)

# 2 b)
fish <- read.table('~/Dropbox/Teaching/IntroToML/GPandMixtures/Lab/Fish.dat')
data <- as.matrix(fish)

# Run the EM
# Plot the histogram and density fit
dataGrid <- seq(min(data)-0.2*IQR(data),max(data)+0.2*IQR(data), length = 1000)
par(mfrow = c(2,2))
for (K in 2:4){
  densEst <- rep(0,length(dataGrid))
  initMu <- quantile(data,seq(1/(K+1),1-1/(K+1), length=K)) + 10*rnorm(K)
  EMfit <- mixtureEM(data, K, initMu = initMu, initSigma  = rep(sd(data),K), initPi = rep(1/K,K), tol = 0.0000001)
  for (k in 1:K){
    densEst <- densEst + EMfit$Pi[k]*dnorm(dataGrid, mean=EMfit$Mu[k], sd = sqrt(EMfit$Sigma[k]))
  }
  hist(data, 30, freq = FALSE, ylim = c(0,0.1), xlim = c(min(data)-0.2*IQR(data),max(data)+0.2*IQR(data)))
  lines(dataGrid,densEst, col="red", lwd = 2)
}


