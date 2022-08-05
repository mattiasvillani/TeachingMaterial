%% USER INPUT
lambda = 0.1; % Shrinkage
splineType = 'PolyTrunc'; % Choices: 'PolyTrunc','ThinPlate' or 'Hardy' 
p = 2; % Order of the polynomial. Relevant only when splineType='polyTrunc'
nKnots = 24;
standardizeX = 1; % If standardizeX=1, covariates are standardize to have mean zero and unit variance.

% Simulate sine data
simulate=1;
if simulate
    n = 200; % Number of observations
    sigma = 0.5; % Standard deviations of the errors in the data generating model.
    x = -5 + rand(n,1)*10; % Simulate covariate values from Uniform[-5,5]
    x = sort(x);
    yTrue = sin(2*x); % This is the true mean curve
    y =  yTrue + randn(n,1)*sigma;
end

% Plotting the simulated data and the true mean curve
subplot(2,2,4)
plot(x,y,'ko','markersize',3.5)
hold on
plot(x,yTrue,'b')
xlabel('Range')
ylabel('LogRatio')

% Estimate a polynomial kernel regression model of order p
[fittedVals,xGrid,knots] = SplineReg(y,x,lambda,p,nKnots,splineType,standardizeX);
plot(xGrid,fittedVals,'r','linewidth',2) % Plotting E(y|x)
AddHashMarks(FigHandle,knots) % Adding the rug plot with knot locations
title(['Trunc. polynomial of order ', int2str(p),' with ',int2str(nKnots),' knots ', 'and', ' \lambda = ',num2str(lambda,3)])
legend('Data','Estimated E(y|x)')
axis tight
box off