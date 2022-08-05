function [muDraws,sigma2Draws] = SimPostNormalCondConj(x,mu0,tau20,nu0,sigma20,nIter,sigma2Init)

%% Function to generate a Gibbs sample from the Normal Model with a 
%  conditionally conjugate prior:
%  
%  MODEL:   x1,...,xn ~ N(mu,sigma), theta and sigma unknown.
%  PRIOR:   mu ~ N(mu0,tau0)
%           sigma2 ~ Inv-X2(nu0,sigma20)
%
%  INPUT:   x           (n-by-1 vector)     Data observations
%           mu0         (scalar)            Prior mean of mu
%           tau20       (scalar)            Prior variance of mu
%           nu0         (scalar)            Prior degrees of freedom in Inv-X2 prior for sigma2
%           sigma20     (scalar)            Prior mean of sigma2
%           nIter       (integer)           Number of posterior draws after burn-in
%           sigma2Init  (scalar)            Initial value for sigma2
%
%  OUTPUT:  muDraws     (n-by-1 vector)     Posterior draws of mu, burnin removed.
%           sigma2Draws (n-by-1 vector)     Posterior draws of sigma2, burnin removed.
%

%% Initialization and setting up storage
xBar = mean(x);
s2 = var(x); % Sample variance.
n = length(x); % The number of observations in the sample.
nBurnIn = round(nIter*0.1);
nIter = nIter + nBurnIn; % Adding 10% burn-in draws.
muDraws = zeros(nIter,1); % This is a vector of zeros where we store the mu-draws
sigma2Draws = zeros(nIter,1); % This is a vector of zeros where we store the sigma2-draws
sigma2 = sigma2Init; % Initial value for sigma2

for i = 1:nIter
    
    %% Drawing mu conditional on sigma2
    muPrec = n/sigma2 + 1/tau20;
    w = (n/sigma2)/(n/sigma2 + 1/tau20);
    muMean = w*xBar + (1-w)*mu0;
    muVar = 1/(muPrec);
    mu = muMean  + sqrt(muVar)*randn; % A draw from the conditional posterior mu | sigma2 ~ N(muMean, muVar);
    % Alternative: mu = normrnd(muMean,sqrt(muVar))
    muDraws(i) = mu;
    
    %% Drawing sigma2 conditional on mu
    nuPost = nu0 + n;
    sigma2Post = ((n-1)*s2 + n*(xBar-mu)^2 + nu0*sigma20)/nuPost;
    sigma2 = nuPost*sigma2Post/chi2rnd(nuPost); % This a draw from: sigma2 | mu ~ Inv-X2(nuPost,sigma2Post)
    sigma2Draws(i) = sigma2;
end

% Removing burn-in
muDraws = muDraws(nBurnIn+1:end);
sigma2Draws = sigma2Draws(nBurnIn+1:end);