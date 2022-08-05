function [lambdaOptCV, lambdaOptGCV, dfFitsOptCV, dfFitsOptGCV, CV, GCV] = CVLinearSmootherMany(y,x,model,p,nKnots,modelType,standardizeX,lambdaGrid)

n = length(y);
nLambdas = length(lambdaGrid);
CV = zeros(nLambdas,1);
GCV = zeros(nLambdas,1);
dfFits = zeros(nLambdas,1);
count = 0;
for lambda = lambdaGrid
    if strcmpi(model,'spline')
        [fittedVals,xGrid,knots,dfFit,S] = SplineReg(y,x,lambda,p,nKnots,modelType,standardizeX,nan);
    else
        [fittedVals,xGrid,kernelOnGrid,dfFit,S] = KernelReg(y,x,lambda,p,modelType,standardizeX,nan); % Note the bandwidth b = lambda
    end
    count = count + 1;
    residuals = y - fittedVals;
    CV(count) = mean((residuals./(1-diag(S))).^2);
    GCV(count) = mean((residuals./(1-trace(S)/n)).^2);
    dfFits(count) = dfFit;
end
lambdaOptCV = lambdaGrid(CV==min(CV));
lambdaOptGCV = lambdaGrid(GCV==min(GCV));
dfFitsOptCV = dfFits(CV==min(CV));
dfFitsOptGCV = dfFits(GCV==min(GCV));

