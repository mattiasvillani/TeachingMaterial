function [fittedVals,fittedValsStd,xGrid,knots,dfFit,Slambda,lambdaSample] = SplineRegBayes(y,x,lambda,p,nKnots,splineType,standardizeX,nGridPoints,nIter)

% Prelims

% Makes sure that x and y are column vectors, and sort them according to
% the x's
x = x(:);
y = y(:);
[x,sortIdx] = sort(x);
y = y(sortIdx);

% Normalize x
if standardizeX
    c0 = mean(x);
    c1 = std(x);
else
    c0=0;
    c1=1;
end
x = (x-c0)/c1;

n = size(x,1);
quantiles = linspace(0.5/nKnots,1-0.5/nKnots,nKnots);
knots = prctile(x,100*quantiles)';

% Estimating the spline
X = feval(splineType,x,knots,p); % This the covariate matrix with the spline covariates
D = eye(size(X,2));
D(1:(p+1),1:(p+1)) = zeros(p+1);


% Gibbs sampling
aSigma = 0.001;
bSigma = 0.001;
aSmooth = 0.001;
bSmooth = 100000;

% This code allows the user to specify the prior on lambda through the mean
% and standard deviation of lambda^(2*p)
% mSmooth = 10;
% sSmooth = 5;
% aSmooth = (mSmooth/sSmooth)^2;
% bSmooth = sSmooth^2/mSmooth;


[betaSample,sigma2Sample,lambdaSample] = PostSampSmoothRegr(y,X,p,D,aSigma,bSigma,aSmooth,bSmooth,nIter);


% Setting up the grid of x-values where we evaluate E(y|x)
if isnan(nGridPoints) % Evaluate the fit at the data points in the sample.
    xGrid = x;
else
    xGrid = linspace(min(x),max(x),nGridPoints)';
end
XGrid = feval(splineType,xGrid,knots,p); %

% Computing the fit over a grid of x-values
fittedVals = zeros(length(xGrid),nIter);
dfFit = zeros(nIter,1);
for i = 1:nIter
    fittedVals(:,i) = XGrid*betaSample(i,:)'; % Here we store the fits
    Slambda = X* ((X'*X + (lambdaSample(i)^(2*p))*D) \ X');
    dfFit(i) = trace(Slambda); % Degrees of freedom for the fit
end
fittedValsStd = std(fittedVals,1,2);
fittedVals = mean(fittedVals,2);

% Transform xGrid back to the original scale
xGrid = c0 + c1*xGrid;
knots = c0 + c1*knots;

end % end of function. Needed only when the m-file contains subfunctions, as below.

function X = PolyTrunc(x,knots,p)

n = length(x);
nKnots = length(knots);
X = [ones(n,1) zeros(n,p+nKnots)];
for j=1:p
    X(:,j+1) = x.^j;
end
for i = 1:nKnots
    idx = (x>knots(i)); % logical index for the data points where x>knot
    X(idx,i+p+1) = (x(idx)-knots(i)).^p;
end

end % PolyTrunc
