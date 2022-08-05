%% Main program to do Gibbs sampling for the normal model with a
%% conditionally conjugate prior

% Simulate a data set from the data generating process N(2,sqrt(2)^2)
x = 2 + randn(1000,1)*sqrt(2);

disp(['Data mean:    ',num2str(mean(x))])
disp(['Data Variance:  ',num2str(var(x))])

% Set up the prior and algorithmic input (nIter)
mu0 = 3; % Prior mean of mu
tau20 = 100; % Prior variance of mu
nu0 = 4; % Prior degrees of freedom for sigma2
sigma20 = 3; % Prior mean of sigma2
nIter = 10000; % Number of Gibbs sample draws (10% is automatically added as burn-in)
sigma2Init = 1; % Initial value for sigma2
[muDraws,sigma2Draws] = SimPostNormalCondConj(x,mu0,tau20,nu0,sigma20,nIter,sigma2Init);

disp(['Posterior mean of mu:    ',num2str(mean(muDraws))])
disp(['Posterior mean of sigma2:    ',num2str(mean(sigma2Draws))])

