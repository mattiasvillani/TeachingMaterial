function [betaSample,sigma2Sample,lambdaSample] = PostSampSmoothRegr(y,X,polyOrder,D,aSigma,bSigma,aSmooth,bSmooth,nIter)



% Preliminaries and storage
[n,p] = size(X); % n is then number of observations. p is the number of regressors (covariates)
nKnots = trace(D);
nOrig = p - nKnots; % Number of original variables
sigma2Sample = zeros(nIter,1);
betaSample = zeros(nIter,p);
lambdaSample = zeros(nIter,1);
lambda2 = 1; % Just some starting value
R = y'*y - y'*X*inv( X'*X + (lambda2^polyOrder)*D)*X'*y; % Residual sum of squares
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nIter
   
    % Updating the some quantities based on the newly simulated lambda
    H = (lambda2^polyOrder)*D; % This is the precision (inverse variance)
    betaTilde = (X'*X + H) \ X'*y; % Bayes estimate of regression coefficients
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        Draw sigma2 and beta jointly from p(beta,sigma2 | lambda, Data)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % First, sigma2 from p(sigma2 | lambda, Data)

    precision = gamrnd((n+2*aSigma-nOrig)/2,2/(R+2*bSigma));
    sigma2 = 1/precision;
    sigma2Sample(i) = sigma2;
   
    % Second, draw beta from p(beta | sigma2, lambda, Data)
    postCovBeta = sigma2*inv(X'*X + H);
    postMeanBeta = betaTilde;
    
    betaVect = mvnrnd(postMeanBeta,postCovBeta)';
    betaSample(i,:) = betaVect';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Draw lambda from p(lambda | beta,sigma2, data)                             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lambda2p = gamrnd(nKnots/2+aSmooth,1/(0.5*precision*betaVect'*D*betaVect + 1/bSmooth));
    lambda2 = (lambda2p)^(1/polyOrder);
    lambdaSample(i) = sqrt(lambda2);
    
end

    
    
    
    
