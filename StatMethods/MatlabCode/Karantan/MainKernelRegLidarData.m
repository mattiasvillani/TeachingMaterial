%% USER INPUT
b = 1.49; % Kernel bandwidth
p = 2; % Order of the polynomial
kernelType = 'Epanechnikov'; % Choices: 'Uniform','Gaussian' or 'Epanechnikov' 


% Loading the Lidar data from file
load('LidarData.mat','y','X')
x = X;

% Plotting the simulated data and the true mean curve (data has to be
% sorted to do this)
plot(x,y,'k.')
hold on
xlabel('Range')
ylabel('LogRatio')

% Estimate a polynomial kernel regression model of order p
[fittedVals,xGrid,kernelOnGrid] = KernelReg(y,x,b,p,kernelType);
plot(xGrid,fittedVals,'r') % Plotting E(y|x)

% Plotting the kernel function as a shaded area
if b < 1.5
    notNan = ~isnan(kernelOnGrid);
    kernelOnGrid = kernelOnGrid(:)';
    yLimits = get(gca,'ylim');
    if strcmpi(kernelType,'Uniform')
        patchHandle = patch([xGrid(notNan) fliplr(xGrid(notNan))],[((min(kernelOnGrid(notNan))+yLimits(2))/2)*ones(size(kernelOnGrid(notNan))) fliplr(kernelOnGrid(notNan))],[0.749 0.749 0],'linestyle','none');
    else
        patchHandle = patch([xGrid(notNan) fliplr(xGrid(notNan))],[min(kernelOnGrid(notNan))*ones(size(kernelOnGrid(notNan))) fliplr(kernelOnGrid(notNan))],[0.749 0.749 0],'linestyle','none');
    end
    set(patchHandle,'linestyle','none','faceLighting','phong','facealpha',0.5,'edgecolor',[0.749 0.749 0])
end
legend('Data','Estimated E(y|x)','Kernel function')

