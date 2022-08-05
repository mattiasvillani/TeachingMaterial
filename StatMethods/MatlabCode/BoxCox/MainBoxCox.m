%% This is a simple script that simulates data and then tries to find the MLE of the Box-Cox parameter.

%% Simulate data
n = 100; % Number of observations
sigma = 0.1;
X = [ones(n,1) rand(n,1)];
y = exp(X*[0 4]').*(exp(sigma*randn(n,1))); % Exponential data, here lambaOpt should be close to zero.
%y = X*[2 2]' + sigma*randn(n,1); % Linear data, here lambaOpt should be close to one.
load LidarData.mat;y = 10 + y;n = size(X,1);X = [ones(n,1) X];

p = size(X,2);

%% Find the optimal Box-Cox transformation
geoMeanY = geomean(y);
H = eye(n) - X*inv(X'*X)*X';
OptInitVal = 1;
Tolerance = 10^-8;
[minNegLik,lambdaOpt,gh,InvHessian] = csminwel('BoxCoxLogLikeLambda',OptInitVal,10^(-4),[],Tolerance,...
    1000000,y,geoMeanY,H);


%% Plot the original and transformed data for the MLE of the Box-Cox parameter
figure('name','A look at the Box-Cox transformation')
subplot(2,1,1)
plot(X(:,2),y,'.')
title('Original data')

subplot(2,1,2)
plot(X(:,2),BoxCoxTrans(y,lambdaOpt),'.r')
title(['Transformed data, \lambda = ',num2str(lambdaOpt)])