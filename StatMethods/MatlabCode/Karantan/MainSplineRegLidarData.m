%% USER INPUT
lambda = 0.1; % Shrinkage
splineType = 'PolyTrunc'; % Choices: 'PolyTrunc','ThinPlate' or 'Hardy' 
p = 1; % Order of the polynomial. Relevant only when splineType='polyTrunc'
nKnots = 24;
standardizeX = 1; % If standardizeX=1, covariates are standardize to have mean zero and unit variance.

% Loading the Lidar data from file
load('CanadianWages.mat','y','X')
x = X;

% Plotting the data
FigHandle = figure;
plot(x,y,'k.')
hold on
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