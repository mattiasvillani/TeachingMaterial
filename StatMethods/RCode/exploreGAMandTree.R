# install.packages("gam") # downloading the gam package from the internet
library(gam) # Loading the gam package into memory

# Loading data from file and separating into training and testing data
spamData<-read.table("SpamReduced.dat",header=TRUE)  # Spam data from Hastie et al.
spam <- spamData$spam;
nTotalObs <- length(spam)
trainingData <-  sample(1:nTotalObs,3000); # Choosing a random subset of the data to train/estimate the model
testData <- setdiff(1:nTotalObs,trainingData) # These are the observations not in the training data.


# Estimating a GAM model with the back-fitting algorithm
gamObj <- gam(spam ~ s(CapRunTotal,df = 4) + s(our,df = 4) + s(over,df = 4)  + s(remove,df = 4) + s(internet,df = 4) + s(free,df = 4) + s(hpl,df = 4) + s(X.,df = 4)+ s(X..1,df = 4) + s(CapRunMax,df = 4) + s(CapRunTotal,df = 4) + s(hp,df = 4) + s(george,df = 4) + s(X1999,df = 4) + s(re,df = 4) + s(edu,df = 4), family = binomial(link="logit"), data = spamData, subset = trainingData)

summary(gamObj)
                                     
# Predictions GAM
predProbsGAM <- predict(gamObj, newdata = spamData[testData,], type="response")

L01 <- 10; # Loss of predicting spam when the email is good
L10 <- 1;  # Loss of predicting good email when the email is spam

preds <- (predProbsGAM > L01/(L01 + L10));
spamTestData <- spamData[testData,1]; # This is the observations on the response in the test data.
p11 = sum((preds*spamTestData)/sum(spamTestData))  # Pr(signal spam | email is spam)
p10 = sum((preds*(1-spamTestData))/sum(1-spamTestData))  # Pr(signal spam | email is good)
p01 = sum(((1-preds)*spamTestData)/sum(spamTestData))    # Pr(signal good email | email is spam)
p00 = sum(((1-preds)*(1-spamTestData))/sum(1-spamTestData))    # Pr(signal good email | email is good)

# Growing a tree
install.packages("tree") # downloading the tree package from the internet
library(tree) # Loading the tree package into memory
treeObj <- tree(spam ~ CapRunTotal + our + over  + remove + internet + free + hpl + X. + X..1 + CapRunMax + CapRunTotal + hp + george + X1999 + re + edu, data = spamData, subset = trainingData)

summary(treeObj)
                                     
# Predictions
predProbsTree <- predict(treeObj, newdata = spamData[testData,])

predsTree <- (predProbsTree > L01/(L01 + L10));
spamTestData <- spamData[testData,1]; # This is the observations on the response in the test data.
p11Tree = sum((predsTree*spamTestData)/sum(spamTestData))  # Pr(signal spam | email is spam)
p10Tree = sum((predsTree*(1-spamTestData))/sum(1-spamTestData))  # Pr(signal spam | email is good)
p01Tree = sum(((1-predsTree)*spamTestData)/sum(spamTestData))    # Pr(signal good email | email is spam)
p00Tree = sum(((1-predsTree)*(1-spamTestData))/sum(1-spamTestData))    # Pr(signal good email | email is good)
