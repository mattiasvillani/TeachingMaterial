function negLogLik = BoxCoxLogLikeLambda(lambda,y,geoMeanY,H)

% This is the (negative) log-Likelihood of the Box-Cox parameter lambda, with all the
% other parameters maximized out (so we are here computing the
% concentrated, or profile, likelihood.)

n = length(y);
yTilde = BoxCoxTrans(y,lambda)/(geoMeanY^(lambda-1));
negLogLik = (yTilde'*H*yTilde)/n;

