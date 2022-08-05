%% Kernel regression and spline regression with a single covariate.
% This is code I have used for teaching the course Statistical Methods at Stockholm University.
% This is in no way production-ready code. It was ment for TEACHING.
%
% Usage: The code estimates spline or kernel regression models for
% situation with a single covariate (again, the code is for TEACHING).
% A nice feature is that many of the inputs (e.g. nKnots) can be vectors
% and the code then loops over all values for that input. You can change 
% only one input at the time. Inputs that can be varied are: lambda, nKnots and polyOrder
%
%
% Author: Mattias Villani, Linköping University
%         mattias.villani@gmail.com
%

%%%%%% BEGIN USER INPUT %%%%%%
DataFileName = 'CanadianWages.mat'; % This is the name of the .mat file that holds the data
                                % y, X, XName (cell), yName (string)
model = 'Spline';               % Choices: 'Spline' or 'Kernel' (Bayesian estimation only for splines)
lambda = 'bayes';               % Shrinkage (Splines) or kernel bandwidth (kernel regression). 
                                % If lambda = 'cv', lambda is determined by (leave-one-out) cross-validation.
                                % If lambda = 'bayes', the Bayesian posterior distribution of lambda is computed.
splineType = 'PolyTrunc';       % Choices: 'PolyTrunc'
kernelType = 'Epanechnikov';    % Choices: 'Uniform','Gaussian' or 'Epanechnikov' 
polyOrder = [2];                          % Order of the polynomial. For splines this is only relevant when splineType='polyTrunc'
nKnots = [10];          % Number of knots
standardizeX = 1;               % If standardizeX = 1, covariates are standardized to have mean zero and unit variance.
nGridPoints = 100;              % The number of grid points (in X-space) where the fit is evaluated. nGridPoints = nan evaluates at all the points in the data.
plotFit = 1;                    % If plotFit=0, no plots are produced.
crossValidate = 1;              % If =1, then cross-validation of lambda
nIter = 5000;                   % The number of draws from the Bayesian posterior distribution
%%%%%   END USER INPUT  %%%%%%

% Choosing the 
switch lower(lambda(1))
    case 'c'
        CVOptimalLambda = 1;
        bayes = 0;
        lambda = nan;
    case 'b'
        CVOptimalLambda = 0;
        bayes = 1;
        lambda = nan;
    otherwise
        % Lambda is numerical.
        CVOptimalLambda = 0;
        bayes = 0;
end

% Loading the data from file DataFileName.mat
load(DataFileName,'y','X','XName','yName')
x = X(:,1); % Change here and the next line if some other column of X is the predictor
xName = XName{1};

% Setting up a grid of shrinkage factors for plotting the cross-validation
% of the lambda parameter
if standardizeX
    lambdaGrid = linspace(0.1,3,100);
else
    lambdaGrid = linspace(0.1,max(x),100);
end

if strcmpi(model,'Spline')
    modelType = splineType;
    FigHandle = figure('name',['Nonparametric spline regression for the ',DataFileName,' data']);
else
    modelType = kernelType;
    FigHandle = figure('name',['Nonparametric kernel regression for the ',DataFileName,' data']);
end
if crossValidate
    figHandleCV = figure('name','Cross-validation of lambda');
end
if bayes
    figHandlePostLambda = figure('name','Posterior distribution of lambda');
    figHandelPostDF = figure('name','Posterior distribution of degrees of freedom');
end

[changedVariable,Values] = ChangedVarSettings(lambda,nKnots,polyOrder,modelType); % Function that figures out which 
                                                                          % parameter you want to vary
nSettings = length(Values);
[nRows,nCols] = ConstructOptimalSubplot(nSettings); % Trivial function that finds out the best way to organize the plots.
if bayes
    dfFits = zeros(nSettings,nIter);
else
    dfFits = zeros(nSettings,1);
end

% Looping over the settings
for i = 1:nSettings
    
    % Changing a setting/parameter of the model
    if ~isnan(changedVariable)
        changeStr = [changedVariable,'=',num2str(Values(i))]; % Building up a string (letters) with the changed setting
        fprintf(1,changeStr) % Print the changed setting to the screen
        eval(changeStr); % Here we actaully change the setting
    end
    % Plotting the data
    if plotFit
        figure(FigHandle)
        subplot(nRows,nCols,i)
        plot(x,y,'k.')
        hold on
        xlabel(xName)
        ylabel(yName)
    end
    
    if CVOptimalLambda % CV-optimal lambda used
        lambdaOptCV = CVLinearSmootherMany(y,x,model,polyOrder,nKnots,modelType,standardizeX,lambdaGrid);
        lambda = lambdaOptCV;
    end
    
    if strcmpi(model,'Kernel')
        % Estimate a polynomial kernel regression model of order polyOrder
        [fittedVals,xGrid,kernelOnGrid,dfFit,Slambda] = KernelReg(y,x,lambda,polyOrder,modelType,standardizeX,nGridPoints); % Note the bandwidth b = lambda
    else
        % Estimate a spline regression model of order polyOrder
        if bayes
            [fittedVals,fittedValsStd,xGrid,knots,dfFit,Slambda,lambdaSample] = SplineRegBayes(y,x,lambda,polyOrder,nKnots,modelType,standardizeX,nGridPoints,nIter);
            dfFits(i,:) = dfFit';
            dfFit = mean(dfFit);
            lambda = mean(lambdaSample);          
        else % Non-Bayesian estimation
            [fittedVals,xGrid,knots,dfFit,Slambda] = SplineReg(y,x,lambda,polyOrder,nKnots,modelType,standardizeX,nGridPoints);
            fittedValsStd = nan(size(fittedVals)); % TODO: Standard errors for fit not yet implemented.
            dfFits(i) = dfFit;
        end
    end
    
    if plotFit
        plot(xGrid,fittedVals,'r','linewidth',2) % Plotting E(y|x)
        hold on
        plot(xGrid,fittedVals - 1.96*fittedValsStd,'linestyle','--','color','g','linewidth',2) % Plotting E(y|x)
        plot(xGrid,fittedVals + 1.96*fittedValsStd,'linestyle','--','color','g','linewidth',2) % Plotting E(y|x)

        % Adding marks for the knot locations
        if strcmpi(model,'Spline')
            AddHashMarks(FigHandle,knots) % Adding the rug plot with knot locations
            title(['PolyOrder = ', int2str(polyOrder),'. nKnots =  ',int2str(nKnots), '. \lambda = ',num2str(lambda,3), '. df = ',num2str(dfFit,2)])
        else % Kernel reg
            title(['PolyOrder = ', int2str(polyOrder),'. nKnots =  ',int2str(nKnots), '. \lambda = ',num2str(lambda,3), '. df = ',num2str(dfFit,2)])
        end
        % DF = SmoothingMatPlot(Slambda,x); %% Plotting the smoother matrix, L.
        
        if i==1
            legend('Data','Estimated E(y|x)','95% error bands for E(y|x)')
        end
        axis tight
        box off
    end 
    
    % Plotting histogram of the posterior of lambda and Degrees of freedom of the fit
    if bayes
        figure(figHandlePostLambda)
        subplot(nRows,nCols,i)
        hist(lambdaSample,30)
        xlabel('\lambda')
        
        figure(figHandelPostDF)
        subplot(nRows,nCols,i)
        hist(dfFits(i,:),30)
        xlabel('Degrees of freedom')
    end
    
    if crossValidate
        figure(figHandleCV)
        subplot(nRows,nCols,i)
        [lambdaOptCV, lambdaOptGCV, dfFitsOptCV, dfFitsOptGCV, CV, GCV] = CVLinearSmootherMany(y,x,model,polyOrder,nKnots,modelType,standardizeX,lambdaGrid);
        plot(lambdaGrid',CV),set(gca,'xscale','log')
        hold on
        plot(lambdaGrid',GCV,'r--')
        line([lambdaOptCV lambdaOptCV],get(gca,'ylim'));
        lineHandle = line([lambdaOptGCV lambdaOptGCV],get(gca,'ylim'));
        set(lineHandle,'linestyle','--','color','r')
        set(gca,'xscale','log')
        xlabel('lambda (log scale)')
        ylabel('CV/GCV')
        if strcmpi(model,'Spline')
            title(['PolyOrder = ', int2str(polyOrder),'. nKnots =  ',int2str(nKnots), '. df = ',num2str(dfFitsOptCV,2)])
        else
            title(['PolyOrder = ', int2str(polyOrder),'. nKnots =  ',int2str(nKnots)])
        end
    end % end if crossValidate
end % end for i = 1:nSettings



