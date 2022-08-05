%% USER INPUT
probit = 0;         % If Probit = 0, logistic model is used.
chooseCov = [1:16]; % Here we choose which covariates to include in the model
tau = 10000;       % Prior scaling factor such that Prior Covariance = (tau^2)*I
optAlgo = 1;        % If optAlgo = 1 we use Matlab's fminunc.m for the optimization, if optAlgo = 0 we use csminwel.m
tolerance = 10^-8;   % The tolerance of the optimizer. Small number make the optimization more precise.

% Loading data from ascii file
data = load('SpamReducedNoLabels.dat');
covNames = {'spam','our','over','remove','internet','free','hpl','!','$','CapRunMax',...
    'CapRunTotal','const','hp','george','1999','re','edu'};
y = data(:,1);
X = data(:,2:end);
X = X(:,chooseCov);
covNames = covNames(1+chooseCov);
nPara = size(X,2);
mu = zeros(nPara,1);    % Prior mean vector
Sigma = tau^2*eye(nPara);

initVal = X \ y; % Initial values by OLS. Note this the same as inv(X'*X)*X'*y, but faster and more stable.

if probit == 1
    logPostFunc = 'LogPostProbit';
else % Logistic
    logPostFunc = 'LogPostLogistic';
end

if optAlgo == 1 % fminunc
    options = optimset('MaxFunEvals',1000000,'Diagnostics','off','LargeScale','off','TolFun',tolerance);
    [postMode,fh,exitFlag,output,grad,hessian] = fminunc(logPostFunc,initVal,options,y,X,mu,Sigma,1);
    logPostMax=-fh;
    postCovApprox = inv(hessian); % Note we are minimizing -logPosterior, so the hessian is actually -Hessian ...
else % csminwel, alternative to fminunc written by Chris Sims at Princeton.
    [fh,postMode,gradient,invHessian,itct,fcount,exitFlag] = csminwel(logPostFunc,initVal,10^(-4)*eye(nPara),...
        [],tolerance,1000000,y,X,mu,Sigma,1);
    logPostMax=-fh;
    postCovApprox = invHessian; % csminwel returns the inverse hessian of the minimization of -LogPosterior.
                                % This is then the same as -inv(H), where H is the hessian in the maximization of
                                % logPosterior...
end

%% Checking out the optimum graphically. Normal approximation is superimposed. 
CheckOpt(logPostFunc,postMode,[],postCovApprox,20,1,[],covNames,y,X,mu,Sigma,0)


% Print the posterior mode and approximate standard deviation to the screen
disp('Log Posterior PDF at the mode: ')
disp(logPostMax)
disp('Posterior mode: ')
disp(postMode')
disp('Approximate posterior standard deviation :')
disp(sqrt(diag(postCovApprox))')
