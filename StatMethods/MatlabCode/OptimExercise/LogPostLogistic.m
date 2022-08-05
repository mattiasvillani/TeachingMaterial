function logPostValue = LogPostLogistic(betaVect,y,X,mu,Sigma,negative)

nPara = length(betaVect);
linPred = X*betaVect;
logLik = sum( linPred.*y -log(1 + exp(linPred)));
if isinf(logLik)
    logLik = -20000;
end
logPrior = log(mvnpdf(betaVect, mu, Sigma));
logPostValue = logLik + logPrior;
if negative
    logPostValue = -logPostValue; % Note Matlabs fminunc.m only minimizes, hence the minus sign here.
end