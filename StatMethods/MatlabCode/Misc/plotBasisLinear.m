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

subplot(3,1,2)
plot(x,1*ones(size(x)),'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Constant basis function')

subplot(3,1,3)
plot(x,0.5*x,'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Linear basis function')


% Quadratic model

x = (0:0.01:1);
figure
subplot(4,1,1)
hold on
meanLine = 1+0.5*x+0.5*x.^2;
plot(x,meanLine,'r','linewidth',2)
plot(x,meanLine+0.1*randn(size(x)),'ko','markersize',6)
ylabel('y')
xlabel('x')
title('Quadratic regression')

subplot(4,1,2)
plot(x,1*ones(size(x)),'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Constant basis function')

subplot(4,1,3)
plot(x,0.5*x,'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Linear basis function')

subplot(4,1,4)
plot(x,0.5*x.^2,'b','linewidth',2)
ylabel('y')
xlabel('x')
title('Quadratic basis function')

