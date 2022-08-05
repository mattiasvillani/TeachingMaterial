% Settings
PriorMean = 0;
PriorStd = 10;
NIter = 100;    % Number of Gibbs sampling draws from posterior
NRepIter = 100;  % Number of data replications for each posterior draw.
figure

% Normal data x1,...,xn ~iid N(2,0.5²)
Data = 2 + randn(500,1)*0.5; % Generate the data
[XRepStatSample,XObsStat,MuRes,SigmaRes]=PostPredAnalysisNormal(Data,PriorMean,PriorStd,NIter,NRepIter); % Sample from the posterior predictive distribution

subplot(2,2,1)
hold on
[nInBins,BinPos] = hist(XRepStatSample,50);
bar(BinPos,nInBins,'y')
plotHandle = plot(XObsStat,0,'or');
set(plotHandle,'markersize',10,'markerfacecolor','r')
title('T(y) = max(y). Normal Data - Normal Model')
xlabel('Max')

% Student t data x1,...,xn ~iid t(2,0.5²) with 10 degrees of freedom
Data = 2 + trnd(10*ones(500,1))*0.5;
[XRepStatSample,XObsStat,MuRes,SigmaRes]=PostPredAnalysisNormal(Data,PriorMean,PriorStd,NIter,NRepIter); % Sample from the posterior predictive distribution

subplot(2,2,2)
hold on
[nInBins,BinPos] = hist(XRepStatSample,50);
bar(BinPos,nInBins,'y')
plotHandle = plot(XObsStat,0,'or');
set(plotHandle,'markersize',10,'markerfacecolor','r')
title('T(y) = max(y). t(10) Data - Normal Model')
xlabel('Max')

% Student t data x1,...,xn ~iid t(2,0.5²) with 4 degrees of freedom
Data = 2 + trnd(4*ones(500,1))*0.5;
[XRepStatSample,XObsStat,MuRes,SigmaRes]=PostPredAnalysisNormal(Data,PriorMean,PriorStd,NIter,NRepIter); % Sample from the posterior predictive distribution

subplot(2,2,3)
hold on
[nInBins,BinPos] = hist(XRepStatSample,50);
bar(BinPos,nInBins,'y')
plotHandle = plot(XObsStat,0,'or');
set(plotHandle,'markersize',10,'markerfacecolor','r')
title('T(y) = max(y). t(4) Data - Normal Model')
xlabel('Max')


% Cauchy data x1,...,xn ~iid t(2,0.5²)
Data = 2 + trnd(1*ones(500,1))*0.5;
[XRepStatSample,XObsStat,MuRes,SigmaRes]=PostPredAnalysisNormal(Data,PriorMean,PriorStd,NIter,NRepIter); % Sample from the posterior predictive distribution

subplot(2,2,4)
hold on
[nInBins,BinPos] = hist(XRepStatSample,50);
bar(BinPos,nInBins,'y')
plotHandle = plot(XObsStat,0,'or');
set(plotHandle,'markersize',10,'markerfacecolor','r')
title('T(y) = max(y). Cauchy Data - Normal Model')
xlabel('Max')

