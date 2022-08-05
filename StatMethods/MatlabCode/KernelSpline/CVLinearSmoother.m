function [CV, GCV] = CVLinearSmoother(y,fittedVals,S)

residuals = y - fittedVals;
CV = mean((residuals./(1-diag(S))).^2);
GCV = mean((residuals./(1-trace(S)/n)).^2);