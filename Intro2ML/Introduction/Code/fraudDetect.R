# Logistic regression for banknote fraud data from http://archive.ics.uci.edu/ml/datasets/banknote+authentication#
# Author: Mattias Villani, Statistics and Machine Learning, Linkoping University, Sweden. e-mail: mattias.villani@liu.se

# Just some figure settings
#install.packages("RColorBrewer")
library("RColorBrewer")
plotColors = brewer.pal(12, "Paired")
pointColor = plotColors[5] # Color for single dots
lwdDef = 8                 # Default line thickness
lwdThin = 6
lwdThinner = 3
pointSizeDef = 4
cexLabDef = 1.5            # Default scaling of font size labels
cexAxisDef = 1.5           # Default scaling of tick labels

# install.packages("mvtnorm") # Loading a package that contains the multivariate normal pdf
library("mvtnorm") # This command reads the mvtnorm package into R's memory. NOW we can use dmvnorm function.

tau = 10
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- tau^2*diag(nPara);

# Defining the log-likelihood function
LogPost <- function(betaVect, y, X, mu, Sigma){
  nPara = length(mu)
  linPred = X%*%betaVect # matrix product
  logLik <- sum( linPred*y -log(1 + exp(linPred)))
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE)
  return(logLik + logPrior)
}


# Reading in fraud data from file
data <- read.csv('~/Dropbox/Teaching/ProbStatUProg/Data/banknoteFraud.csv', header = FALSE)
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
SelectTraining <- sample(1:dim(data)[1], size = 1000, replace = FALSE)
y <- data[,5]
X <- as.matrix(cbind(1,data[,1:4])) # Adding a column of ones for the intercept
nPara <- dim(X)[2]       # Number of covariates incl intercept
yTrain <- y[SelectTraining]
yTest <- y[-SelectTraining]
XTrain <- X[SelectTraining,]
XTest <- X[-SelectTraining,]

# Optimize to the find the ML estimates. control = list(fnscale=-1) puts a minus sign in front of LogLik
# We need this since optim minimizes, not maximizes.
initVal <- as.vector(solve(crossprod(XTrain,XTrain))%*%t(XTrain)%*%yTrain); # Initial values by OLS
optimResults <- optim(initVal, LogPost, gr = NULL, yTrain, XTrain, 
                      mu, Sigma, method=c("BFGS"), control=list(fnscale=-1), hessian=TRUE)
modePost <- optimResults$par
CovPost <- -solve(optimResults$hessian)
approxPostStd <- sqrt(diag(CovPost)) # Computing approximate standard deviations.


# Predition function
PredictLogistic <- function(threshold = 0.5, XTest, yTest, betaHat){
  linFunc = XTest%*%betaHat # matrix product
  thetaVect = exp(linFunc)/(1+exp(linFunc))
  results = list()
  results$probs <- thetaVect 
  results$preds <- ifelse(thetaVect>threshold,1,0)
  results$confusionMatrix <- table(results$preds,yTest)
  results$accuracy = sum(diag(results$confusionMatrix))/dim(XTest)[1] # What proportion of notes were correctly classified?
  results$precision = results$confusionMatrix[2,2]/sum(results$confusionMatrix[2,]) # Out of those selected (marked as fraud) what proportion were right? 
  results$recall = results$confusionMatrix[2,2]/sum(results$confusionMatrix[,2])# What proportion of frauds were detected? Sensitivity. True positive rate.
  results$FPR = results$confusionMatrix[2,1]/sum(results$confusionMatrix[,1]) # False Positive Rate
  return(results)
}

# Predition function
PredictLogisticBayes <- function(threshold = 0.5, XTest, yTest, betaSample){
  nSim = dim(betaSample)[1]
  thetaVect = matrix(NA,dim(XTest)[1],nSim)
  for (i in 1:nSim){
    linFunc = XTest%*%betaSample[i,] # matrix product
    thetaVect[,i] = exp(linFunc)/(1+exp(linFunc))
  }
  return(thetaVect)
}
 
# Predicting the test set and evaluating the results with threshold = 0.5
results <- PredictLogistic(threshold = 0.5, XTest, yTest, modePost)


par(cex.lab=cexLabDef, cex.axis = cexAxisDef)
#par(mfrow = c(2,3))
plot(results$probs, axes=FALSE, pch = 16,  col = plotColors[2], xlab = "Obs number", ylab = "Prediction probability")
axis(side = 1, at = c(0,100,200,300,400))
axis(side = 2, at = seq(0, 1, by = 0.2))
falsePosIndex = which(results$preds!=yTest)
points(falsePosIndex,results$probs[falsePosIndex], col = plotColors[7], pch = 16)


# Sample from posterior
nSim = 1000
betaSample = rmvnorm(nSim, mean = modePost, sigma = CovPost)
predSample = PredictLogisticBayes(threshold = 0.5, XTest, yTest, betaSample)
predUpper = apply(predSample, 1, quantile,0.975)
predLower = apply(predSample, 1, quantile,0.025)

arrows(falsePosIndex, predLower[falsePosIndex], 
       falsePosIndex, predUpper[falsePosIndex], 
       length=0.05, angle=90, code=3, col = plotColors[7])

ROC <-function(PredictFunction, thresholds, ...){
  # PredictFunction is a function that returns a list 'results' with elements results$precision and results$recall
  # First element of PredictFunction must be threshold, the other arguments can be anything needed for that function (hence the  ... argument)
  nThres <- length(thresholds)
  precision <- rep(NA,nThres)
  recall <- rep(NA,nThres)
  count <- 0 
  for (threshold in thresholds){
    count = count + 1
    results <- PredictFunction(threshold = threshold, ...)
    precision[count] <- results$precision
    recall[count] <- results$recall # This is also True Positive Rate (TPR)
    FPR <- results$confusionMatrix[2,1]/sum(results$confusionMatrix[,1]) # False Positive Rate
  }
  return(cbind(FPR,recall))
}

ROCResults <- ROC(PredictLogistic, thresholds = seq(0.00000001,0.1,length = 10000), XTest, modePost)
plot(ROCResults[,1],ROCResults[,2])
