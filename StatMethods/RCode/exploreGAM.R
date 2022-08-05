# install.packages("gam") # downloading the gam package from the internet
library(gam) # Loading the gam package into memory

# Loading data from file
spamData<-read.table("SpamReduced.dat",header=TRUE)  # Spam data from Hastie et al.

d <- 4 # This is the degrees of freedom for each covariate
gamObj <- gam(spam ~ s(CapRunTotal,df = d) + s(our,df = d) + s(over,df = d)  + s(remove,df = d) + s(internet,df = d) + s(free,df = d) + s(hpl,df = d) + s(X.,df = d)+ s(X..1,df = d) + s(CapRunMax,df = d) + s(CapRunTotal,df = d) + s(hp,df = d) + s(george,df = d) + s(X1999,df = d) + s(re,df = d) + s(edu,df = d), family = binomial, data = spamData)
d <- 4
gamObj <- gam(spam ~ s(CapRunTotal,df = d), family = binomial, data = spamData)

summary(gamObj)

plot(gamObj)

