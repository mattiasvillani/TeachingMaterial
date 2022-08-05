function [changedVariable,Values] = ChangedVarSettings(lambda,nKnots,polyOrder,modelType)

if length(lambda)>1
    changedVariable = 'lambda';
    Values = lambda;
elseif length(nKnots)>1
    changedVariable = 'nKnots';
    Values = nKnots;
elseif length(polyOrder)>1
    changedVariable = 'polyOrder';
    Values = polyOrder;
elseif iscell(modelType)
    changedVariable = 'modelType';
    Values = modelType;
else
    disp('No model setting is varied')
    changedVariable = nan;
    Values = nan;
end