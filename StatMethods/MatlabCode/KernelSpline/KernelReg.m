function [fittedVals,xGrid,kernelOnGrid,dfFit,Slambda] = KernelReg(y,x,b,p,kernelType,standardizeX,nGridPoints)

if nargin == 5
    standardizeX = 1;
end

% Prelims

% Makes sure that x and y are column vectors, and sort them according to the x's
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

% Setting up the grid of x-values where we evaluate E(y|x)
if isnan(nGridPoints) % Evaluate the fit at all data points.
    xGrid = x;
else
    xGrid = linspace(min(x),max(x),nGridPoints);
end
nGridPoints = length(xGrid);
fittedVals = zeros(nGridPoints,1); % Here we store the fits

% Estimating all the local regression fits.
Xlocal = [ones(n,1) zeros(n,p)]; % This is the local covariate matrix, later we modify it
W = eye(n,n); % This is the weight matrix, later we modify it
for i = 1:nGridPoints
    % Setting up the local regression covariates
    x0 = xGrid(i); % This is the x-value where we are currently making the local fit
    for j = 1:p
        Xlocal(:,j+1) = (x-x0).^j;
    end
    W = diag(feval(kernelType,x,x0,b));
    localBeta = (Xlocal'*W*Xlocal) \ (Xlocal'*W*y); % This is a fast and stable way to computing inv(Xlocal'*W*Xlocal)*(Xlocal'*W*y)
    fittedVals(i) = localBeta(1); % The fitted value is just the intercept in this model.
end % nGridPoints

% Computing the smoother matrix, S. Need to make a new pass through the
% x's, but this using the actual x's as grid.
Slambda = zeros(n,n); % Initializing the smoothing matrix.
W = eye(n,n); % This is the weight matrix, later we modify it
for i = 1:n
    % Setting up the local regression covariates
    x0 = x(i); % This is the x-value where we are currently making the local fit
    for j = 1:p
        Xlocal(:,j+1) = (x-x0).^j;
    end
    W = diag(feval(kernelType,x,x0,b));
    Slambda(i,:) = Xlocal(i,:)*((Xlocal'*W*Xlocal) \ (Xlocal'*W)); % This is the ith row of the smoother matrix.
end % loop over observations.
dfFit = trace(Slambda); % Degrees of freedom for the fit


% Computing the kernel function - ignore what's happening here, it is just
% stuff to get the graphical stuff good looking
middlePoint = round(nGridPoints/2);
x0 = xGrid(middlePoint); % Middle of the grid
kernelOnGrid = feval(kernelType,xGrid,x0,b);
nonZeroKernel = (kernelOnGrid~=0);
kernelOnGrid = kernelOnGrid/sum(kernelOnGrid);
kernelOnGrid = kernelOnGrid-kernelOnGrid(middlePoint); % kernelOnGrid is now zero at x0
kernelOnGrid = 0.3*kernelOnGrid.*abs(((min(y))./(min(kernelOnGrid))));
kernelOnGrid = kernelOnGrid + (fittedVals(middlePoint)+min(y))/2;
smallestNonZero = find(nonZeroKernel,1,'first');
largestNonZero = find(nonZeroKernel,1,'last');
kernelOnGrid(1:smallestNonZero-1) = nan;
kernelOnGrid(largestNonZero+1:end) = nan;

% Transform xGrid back to the original scale
xGrid = c0 + c1*xGrid; 

end % end of function. Needed only when the m-file contains subfunctions, as below.

function K = Epanechnikov(x,x0,b)

t = abs(x-x0)/b;
K = zeros(length(t),1);
idx = abs(t) <= 1; % These are the observations numbers for observations where t<=1
K(idx) = (3/4)*(1-t(idx).^2);
end


function K = Gaussian(x,x0,b)

t = abs(x-x0)/b;
K = (1/sqrt(2*pi*b^2))*normpdf(t);
end

function K = Uniform(x,x0,b)

t = abs(x-x0)/b;
K = zeros(length(t),1);
idx = abs(t) <= 1; % These are the observations numbers for observations where t<=1
K(idx) = 1;
end

