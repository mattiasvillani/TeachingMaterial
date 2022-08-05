#################################################################################################
#################  Example of conjugate prior inference of the multinomial model ################ 
#################################################################################################
 

###########   Setting up data and prior  #################
y=c(36,87,77) # Data
alpha=c(1,1,1) # Dirichlet prior hyperparameters
NIter=10000 # Number of posterior draws

###########   Initializing storage matrices and performing posterior sampling from Dirichlet  #################
ThetaDraws=matrix(0,NIter,3) # Matrix where the posterior draws are stored
ThetaDraws[,1]=rgamma(NIter,y[1]+alpha[1],1); # Generating (intermediate) independent Gamma-draws
ThetaDraws[,2]=rgamma(NIter,y[2]+alpha[2],1);
ThetaDraws[,3]=rgamma(NIter,y[3]+alpha[3],1);
for (i in 1:NIter)
{
ThetaDraws[i,]=ThetaDraws[i,]/sum(ThetaDraws[i,]) # Diving every column of ThetaDraws by the sum of the elements in that column. Gives Dirichlet draws.
}


################ Computing Summary statistics from the posterior sample ###################
mean(ThetaDraws[,1])
mean(ThetaDraws[,2])
mean(ThetaDraws[,3])

sqrt(var(ThetaDraws[,1]))
sqrt(var(ThetaDraws[,2]))
sqrt(var(ThetaDraws[,3]))

sum(ThetaDraws[,2]>ThetaDraws[,3])/NIter # p(theta2>theta3|Data)

# Plots histograms of the posterior draws
plot.new() # Opens a new graphical window
par(mfrow=c(2,2)) # Splits the graphical window in four parts (2-by-2 structure)
hist(ThetaDraws[,1],25) # Plots the histogram of theta[,1] in the upper left subgraph
hist(ThetaDraws[,2],25)
hist(ThetaDraws[,3],25)

