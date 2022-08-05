function [fittedVals,xGrid,knots,dfFit,Slambda] = SplineReg(y,x,lambda,p,nKnots,splineType,standardizeX,nGridPoints)

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
betaEst = (X'*X + (lambda^(2*p))*D) \ X'*y; % X \ y = inv(X'*X)*X'*y

% Setting up the grid of x-values where we evaluate E(y|x)
if isnan(nGridPoints) % Evaluate the fit at the data points in the sample.
    xGrid = x;
else
    xGrid = linspace(min(x),max(x),nGridPoints)';
end
XGrid = feval(splineType,xGrid,knots,p); %

% Computing the fit over a grid of x-values
fittedVals = XGrid*betaEst; % Here we store the fits

% Transform xGrid back to the original scale
xGrid = c0 + c1*xGrid;
knots = c0 + c1*knots;

Slambda = X* ((X'*X + (lambda^(2*p))*D) \ X');
dfFit = trace(Slambda); % Degrees of freedom for the fit

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
