function logPostValue = LogPostProbit(betaVect,y,X,mu,Sigma,negative)

linPred = X*betaVect;

logLik = sum( y.*log(normcdf(linPred)) + (1-y).*log(1-normcdf(linPred))  );
if isinf(logLik)
    logLik = -2000;    
end
logPrior = log(mvnpdf(betaVect, mu, Sigma));
logPostValue = logLik + logPrior;
if negative
    logPostValue = -logPostValue; % Note Matlabs fminunc.m only minimizes, hence the minus sign here.
end
