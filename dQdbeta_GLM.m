function deriv=dQdbeta_GLM(ginvd, Q1, x, dbeta0dbeta1)

% Function:
%   Calculate first oder derivative of Q(beta0(beta), beta) with respect to
%   beta
%
% Arguments:
%   ginvd: first order derivative of ginv(u) with respect to u
%   Q1: 
%   x: covariate x (matrix)
%   dbeta0dbeta1: first order derivative of beta 0 with respect to beta
%
% Output:
%   dbeta0dbeta1= first order partial derivative of beta0 with respect to
%   beta
%   dbeta0dbeta2= second order partial derivative of beta0 with respect to
%   beta
n=size(x, 1);
x1=x+repmat(dbeta0dbeta1', n, 1);
deriv=x1'*(ginvd.*Q1);

