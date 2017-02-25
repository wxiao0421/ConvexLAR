function M=hessianQdbeta_GLM(ginvd, ginvdd, Q1, Q11, x, dbeta0dbeta1, dbeta0dbeta2)
% Function:
%   Define the event function used by ode45
%
% Arguments:
%   x: covariates x (matrix)
%   ...
% Output:
%   M: hessian matrix of f(beta)

%M=x'*diag(Q11.*(ginvd.^2)+Q1.*ginvdd)*x;


n=size(x, 1);
x1=x+repmat(dbeta0dbeta1', n, 1);
temp1=Q11.*(ginvd.^2)+Q1.*ginvdd;
temp2=Q1.*ginvd;

M=x1'*diag(temp1)*x1+sum(temp2)*dbeta0dbeta2;
