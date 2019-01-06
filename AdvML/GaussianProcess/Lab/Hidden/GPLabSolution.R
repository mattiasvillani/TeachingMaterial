# Solutions to Gaussian Process lab - Advanced ML
# By Mattias Villani


###############################################
###            Tullinge temp                ### 
###############################################

PeriodicSine <- function(sigmaf = 1, ell = 1, d = 1) 
{
  rval <- function(x, y = NULL) {
    r = sqrt(crossprod(x-y));
    return(sigmaf^2*exp(-2*(sin(pi*r/d)^2)/ell^2))
  }
  class(rval) <- "kernel"
  return(rval)
} 

LocallyPeriodicSine <- function(sigmaf = 1, ellPeriod = 1, ellGlobal = 1, d = 1) 
{
  rval <- function(x, y = NULL) {
    r = sqrt(crossprod(x-y));
    return(sigmaf^2*exp(-2*(sin(pi*r/d)^2)/ellPeriod^2)*exp(-0.5*r^2/ellGlobal^2))
  }
  class(rval) <- "kernel"
  return(rval)
} 

# Plot the correlation function
PeriodicFunc <- LocallyPeriodicSine(sigmaf = 1, ellPeriod = 1, ellGlobal = 10, d = 4)
d = 4;
xGrid <- seq(0,4*d, by = 1/(40*d))
kVals <- matrix(NA,length(xGrid))
for (i in 1:length(xGrid)){
  kVals[i] <- PeriodicFunc(0,xGrid[i])
}
plot(xGrid,kVals,type='l', xlab='r', ylab='k(r)', main = 'Correlation function')



tempData <- read.csv('https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv', header=TRUE, sep=';')
temp <- tempData$temp
plot(temp, type="l")
time = 1:length(temp)
day = rep(1:365,6)

# Extract every k:th observation
subset <- seq(1, length(temp), by = 5)
temp <- temp[subset]
time = time[subset]
day <- day[subset]
plot(time,temp, type="l")

polyFit <- lm(temp ~  time + I(time^2)  + I(time^2))
sigmaNoise = sd(polyFit$residuals)

# GP on time - SE kernel
sigmaf <- 20
ell <- 0.2
GPfit <- gausspr(time, temp, var = sigmaNoise^2) 
GPfit <- gausspr(time, temp, kernel = Matern32, kpar = list(sigmaf = sigmaf, ell=ell), var = sigmaNoise^2) 
meanPred <- predict(GPfit, time)
lines(time, meanPred, col="blue", lwd = 2)


# GP on day - SEkernel
sigmaf <- 20
ell <- 0.2*6 # kernlab standardized data so the distance between two days is (6 times) larger here on the standardized scale than when we used time above.
GPfit <- gausspr(day, temp, kernel = Matern32, kpar = list(sigmaf = sigmaf, ell=ell), var = sigmaNoise^2) 
meanPred <- predict(GPfit, day)
lines(time, meanPred, col="red", lwd = 2)

# GP on time - periodic kernel
sigmaf <- 20
ell <- 0.2
d <- 365/sd(time)
GPfit <- gausspr(time, temp, kernel = PeriodicSine, kpar = list(sigmaf = sigmaf, ell=ell, d = d), var = sigmaNoise^2) 
meanPred <- predict(GPfit, time)
lines(time, meanPred, col="purple", lwd = 2)

# GP on time - locally periodic kernel
sigmaf <- 20
ellPeriod <- 1
ellGlobal <- 10
d <- 365/sd(time)
GPfit <- gausspr(time, temp, kernel = LocallyPeriodicSine, kpar = list(sigmaf = sigmaf, ellPeriod=ellPeriod, ellGlobal = ellGlobal, d = d), var = sigmaNoise^2) 
meanPred <- predict(GPfit, time)
lines(time, meanPred, col="green", lwd = 2)


# Evalating the log marginal likelihood over a grid of values
sigmafGrid <- seq(1, 20, by = 1)
ellGrid <- seq(0.1, 3, by = 0.1)
lml <- matrix(NA, length(ellGrid), length(sigmafGrid))
for (j in 1:length(sigmafGrid)){
  for (i in 1:length(ellGrid)){
    print(c(j,i))
    lml[i,j] <- MargLikeGPReg(X = time, y = temp, kernel = Matern32, kpar = list(sigmaf = sigmafGrid[j], ell = ellGrid[i]), sigmaNoise)
  }
}

# Find the hyperparametes with largest log marginal likelihood over the grid
optEll <- ellGrid[which(lml == max(lml), arr.ind =TRUE)[1]]
optSigmaf <- sigmafGrid[which(lml == max(lml), arr.ind =TRUE)[2]]
message(c("The optimal hyperparameters are: ell = ",optEll," and ", "sigmaf = ",optSigmaf))







##############################################
#####       Banknote fraud           ########
##############################################
data <- read.csv('https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/banknoteFraud.csv', header=FALSE, sep=',')
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])
set.seed(111); SelectTraining <- sample(1:dim(data)[1], size = 1000, replace = FALSE) # Randomly select training samples.
y <- data[,5]
X <- as.matrix(data[,1:4]) # Adding a column of ones for the intercept
nPara <- dim(X)[2]       # Number of covariates incl intercept
yTrain <- y[SelectTraining]
yTest <- y[-SelectTraining]
XTrain <- X[SelectTraining,]
XTest <- X[-SelectTraining,]

# First with only two covariates so that we can plot
selVars <- c(1,2)
GPfitFraud <- gausspr(x = XTrain[,selVars], y = yTrain, kernel = "rbfdot", kpar = 'automatic')
GPfitFraud

# predict on the training set
predTrain <- predict(GPfitFraud,XTrain[,selVars])
table(predTrain, yTrain) # confusion matrix for training sample
accuracyTrain <-sum(predTrain==yTrain)/length(yTrain)
predTest <- predict(GPfitFraud,XTest[,selVars])
table(predTest, yTest) # confusion matrix for test sample
accuracyTest <-sum(predTest==yTest)/length(yTest)
print(c(accuracyTrain,accuracyTest))

# class probabilities 
Xplot <- XTrain[,selVars[1:2]]
yplot <- yTrain

# Make grid points to plot over
x1 <- seq(min(Xplot[,1]),max(Xplot[,1]),length=100)
x2 <- seq(min(Xplot[,2]),max(Xplot[,2]),length=100)
gridPoints <- meshgrid(x1, x2)
gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))
gridPoints <- data.frame(gridPoints)
probPreds <- predict(GPfitFraud, gridPoints, type="probabilities")

# Plotting for Prob(fraud)
contour(x1, x2, matrix(probPreds[,1],100), 20, xlab = "x1", ylab = "x2", main = 'Prob(fraud)')

# Adding data points
points(Xplot[yplot==0,1],Xplot[yplot==0,2], col="red", cex=6, pch=".")
points(Xplot[yplot==1,1],Xplot[yplot==1,2], col="blue", cex=6, pch=".")

# Then with all four covariates
selVars <- c(1,2,3,4)
GPfitFraud <- gausspr(x = XTrain[,selVars], y = yTrain, kernel = "rbfdot", kpar = 'automatic')
GPfitFraud

# predict on the training set
predTrain <- predict(GPfitFraud,XTrain[,selVars])
table(predTrain, yTrain) # confusion matrix for training sample
accuracyTrain <-sum(predTrain==yTrain)/length(yTrain)
predTest <- predict(GPfitFraud,XTest[,selVars])
table(predTest, yTest) # confusion matrix for test sample
accuracyTest <-sum(predTest==yTest)/length(yTest)
print(c(accuracyTrain,accuracyTest))




## THE MARGINAL LIKELIHOOD STUFF BELOW GIVES ME STRANGE RESULTS. IGNORE FOR NOW ##

## Defining the log marginal likelihood function
MargLikeGPReg <- function(X, y, kernel, kpar, sigmaNoise){
  kernelFunc <- do.call(kernel, kpar)
  K <- kernelMatrix(kernel = kernelFunc, x = X) # K(X,X)
  n <- length(y)
  
  L = chol(K+sigmaNoise^2*diag(n));
  alpha_ <- solve(t(L),solve(L,y))
  return( -0.5*crossprod(y,alpha_) -sum(log(diag(L))) )
}

# Evalating the log marginal likelihood over a grid of values
sigmafGrid <- seq(0.1,1,by = 0.2)
ellGrid <- seq(1,10,by = 1)
lml <- matrix(NA, length(ellGrid), length(sigmafGrid))
for (j in 1:length(sigmafGrid)){
  for (i in 1:length(ellGrid)){
    lml[i,j] <- MargLikeGPReg(X = Distance, y = LogRatio, kernel = Matern32, kpar = list(sigmaf = sigmafGrid[j], ell = ellGrid[i]), sigmaNoise)
  }
}
# Find the hyperparametes with largest log marginal likelihood over the grid
optEll <- ellGrid[which(lml == max(lml), arr.ind =TRUE)[1]]
optSigmaf <- sigmafGrid[which(lml == max(lml), arr.ind =TRUE)[2]]
message(c("The optimal hyperparameters are: ell = ",optEll," and ", "sigmaf = ",optSigmaf))

# Plotting the fit based on the optimal hyperparameters
GPfit <- gausspr(Distance, LogRatio, kernel = Matern32, kpar = list(sigmaf = optSigmaf, ell=optEll), var = sigmaNoise^2) 
meanPred <- predict(GPfit, Distance)
lines(Distance, meanPred, col="red", lwd = 3)

