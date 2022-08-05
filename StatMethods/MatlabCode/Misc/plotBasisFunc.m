%% Script to plot basis functions

%% Linear model
x = (0:0.01:1);
figure
subplot(3,1,1)
hold on
meanLine = 1+0.5*x;
plot(x,meanLine,'r','linewidth',2)
plot(x,meanLine+0.1*randn(size(x)),'ko','markersize',6)
ylabel('y')
xlabel('x')
title('Linear regression')
axis tight

subplot(3,1,2)
plot(x,ones(size(x)),'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Constant basis function')
set(gca,'ylim',[0 1.1])

subplot(3,1,3)
plot(x,x,'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Linear basis function')
set(gca,'ylim',[0 1.1])

print linearBasis -depsc

% Quadratic model

x = (0:0.01:1);
figure
subplot(4,1,1)
hold on
meanLine = 1+0.5*x+1*x.^2;
plot(x,meanLine,'r','linewidth',2)
plot(x,meanLine+0.1*randn(size(x)),'ko','markersize',6)
ylabel('y')
xlabel('x')
title('Quadratic regression')
axis tight

subplot(4,1,2)
plot(x,1*ones(size(x)),'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Constant basis function')
set(gca,'ylim',[0 1.1])

subplot(4,1,3)
plot(x,x,'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Linear basis function')
set(gca,'ylim',[0 1.1])

subplot(4,1,4)
plot(x,x.^2,'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Quadratic basis function')
print quadraticBasis -depsc
set(gca,'ylim',[0 1.1])



% Piece-wise linear model

x = (0:0.01:1)';
n = length(x);
figure
subplot(4,1,1)
hold on
knots = [0.4];
nKnots = length(knots);
X = [ones(n,1) x zeros(n,nKnots)];
for i = 1:nKnots
    idx = (x>knots(i)); % logical index for the data points where x>knot
    X(idx,i+2) = x(idx)-knots(i);
end
betaVect = [1 0.5 2]';
meanLine = X*betaVect;
plot(x,meanLine,'r','linewidth',2)
plot(x,meanLine+0.1*randn(size(x)),'ko','markersize',6)
ylabel('y')
xlabel('x')
title('Piece-wise linear regression')
axis tight

subplot(4,1,2)
plot(x,X(:,1),'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Constant basis function')
set(gca,'ylim',[0 1.1])

subplot(4,1,3)
plot(x,X(:,2),'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Linear basis function')
set(gca,'ylim',[0 1.1])

subplot(4,1,4)
plot(x,X(:,3),'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Basis function for (x-0.4)+')
print brokenStickBasis -depsc
set(gca,'ylim',[0 1.1])


% Spline model - linear basis

x = (0:0.01:1)';
n = length(x);
figure
subplot(2,1,1)
hold on
knots = [0.4 0.5 0.6 0.7 0.8];
nKnots = length(knots);
X = [ones(n,1) x zeros(n,nKnots)];
for i = 1:nKnots
    idx = (x>knots(i)); % logical index for the data points where x>knot
    X(idx,i+2) = x(idx)-knots(i);
end
betaVect = [1 0.5 1 0.5 -0.5 -1 -0.5]';
meanLine = X*betaVect;
plot(x,meanLine,'r','linewidth',2)
plot(x,meanLine+0.1*randn(size(x)),'ko','markersize',6)
ylabel('y')
xlabel('x')
title('Spline regression')
axis tight

subplot(2,1,2)
plot(x,X,'linewidth',2)
ylabel('y')
xlabel('x')
title('Basis functions')
set(gca,'ylim',[-0.1 1.1])
print splineBasisLinear -depsc


% Spline model - quadratic basis

x = (0:0.01:1)';
n = length(x);
figure
subplot(2,1,1)
hold on
knots = [0.4 0.5 0.6 0.7 0.8];
nKnots = length(knots);
X = [ones(n,1) x zeros(n,nKnots)];
for i = 1:nKnots
    idx = (x>knots(i)); % logical index for the data points where x>knot
    X(idx,i+2) = (x(idx)-knots(i)).^2;
end
betaVect = [1 0.5 1 0.5 -0.5 -1 -0.5]';
meanLine = X*betaVect;
plot(x,meanLine,'r','linewidth',2)
plot(x,meanLine+0.1*randn(size(x)),'ko','markersize',6)
ylabel('y')
xlabel('x')
title('Spline regression - Truncated Quadratics')
axis tight

subplot(2,1,2)
plot(x,X,'linewidth',2)
ylabel('y')
xlabel('x')
title('Basis functions')
set(gca,'ylim',[-0.1 1.1])
print splineBasisQuadratic -depsc





% Spline model - thin plate spline basis

x = (0:0.01:1)';
n = length(x);
figure
subplot(2,1,1)
hold on
knots = [0.4 0.5 0.6 0.7 0.8];
nKnots = length(knots);
X = [ones(n,1) x zeros(n,nKnots)];
for i = 1:nKnots
    idx = (x>knots(i)); % logical index for the data points where x>knot
    X(:,i+2) = ((x-knots(i)).^2).*log(abs((x-knots(i))));
    X(isnan(X(:,i+2)),i+2) = 0;
end
betaVect = [1 0.5 1 0.5 -0.5 -1 -0.5]';
meanLine = X*betaVect;
plot(x,meanLine,'r','linewidth',2)
plot(x,meanLine+0.1*randn(size(x)),'ko','markersize',6)
ylabel('y')
xlabel('x')
title('Spline regression - Thin plate splines')
axis tight

subplot(2,1,2)
plot(x,X(:,3:end),'linewidth',2)
ylabel('y')
xlabel('x')
title('Basis functions')
set(gca,'ylim',[-0.2 0.1])
print splineBasisThinPlate -depsc







