%% USER INPUT
b = 0.5; % Kernel bandwidth
p = 0; % Order of the polynomial
kernelType = 'Epanechnikov'; % Choices: 'Uniform','Gaussian' or 'Epanechnikov' 


% Simulate data
n = 200; % Number of observations
sigma = 0.5; % Standard deviations of the errors in the data generating model.
x = -5 + rand(n,1)*10; % Simulate covariate values from Uniform[-5,5]
yTrue = sin(x); % This is the true mean curve
y =  yTrue + randn(n,1)*sigma;

% Plotting the simulated data and the true mean curve (data has to be
% sorted to do this)
plot(x,y,'k.')
hold on
[x, sortIdx] = sort(x);
yTrue = yTrue(sortIdx);
y = y(sortIdx);
plot(x,yTrue,'b')


% Estimate a polynomial kernel regression model of order p
[fittedVals,xGrid,kernelOnGrid] = KernelReg(y,x,b,p,kernelType);
plot(xGrid,fittedVals,'r') % Plotting E(y|x)

% Plotting the kernel function as a shaded area.
notNan = ~isnan(kernelOnGrid);
kernelOnGrid = kernelOnGrid(:)';
yLimits = get(gca,'ylim');
if strcmpi(kernelType,'Uniform')
    patchHandle = patch([xGrid(notNan) fliplr(xGrid(notNan))],[((min(kernelOnGrid(notNan))+yLimits(2))/2)*ones(size(kernelOnGrid(notNan))) fliplr(kernelOnGrid(notNan))],[0.749 0.749 0],'linestyle','none');
else
    patchHandle = patch([xGrid(notNan) fliplr(xGrid(notNan))],[min(kernelOnGrid(notNan))*ones(size(kernelOnGrid(notNan))) fliplr(kernelOnGrid(notNan))],[0.749 0.749 0],'linestyle','none');
end
set(patchHandle,'linestyle','none','faceLighting','phong','facealpha',0.5,'edgecolor',[0.749 0.749 0])
legend('Data','True E(y|x)','Estimated E(y|x)','Kernel function')
