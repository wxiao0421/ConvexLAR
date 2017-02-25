function deriv=getderiv_GLM(x, y, actind, wactind, distname)
% Function:
%   Calculate the first-order partial derivative (b(beta)) of f(beta) for
%   GLMs. ('normal', 'poisson', 'binomial')
%
% Arguments:
%   x: covariates x (matrix)
%   y: response y (vector)
%   actind: active index set
%   wactind: current value of active index set
%   distname: name of distribution 
%
% Output:
%   deriv: first-order partial derivative at current beta value

% xw calculate x%*%beta where all inactive variables have coefficients equal to 0
% f is negative log likelihood
xw=x(:, actind)*reshape(wactind, length(wactind), 1);
% beta0 is determined by current beta value
beta0=glmfit(ones(size(x,1), 1), y, distname, 'constant', 'off', 'offset', xw);
% eta is beta0+x%*%beta
eta=xw+beta0;
% value of g inverse function and its 1,2,3 order derivative
[ginv, ginvd, ginvdd, ginvddd]=invlinkDERIV_GLM(eta, distname);
u=ginv;
[Q1, Q11, Q111]=qDERIV_GLM(u, y, distname);
[dbeta0dbeta1, dbeta0dbeta2]=dbeta0dbeta_GLM(x, ginvd, ginvdd, ginvddd, Q1, Q11, Q111);
deriv=-dQdbeta_GLM(ginvd, Q1, x, dbeta0dbeta1);