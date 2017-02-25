function dbetadt=mvDirectionG_GLM(x, y, actind, wactind, distname, weight, t)

% Function:
%   Define the event function used by ode45
%
% Arguments:
%   x: covariates x (matrix)
%   y: response y (vector)
%   actind: active index set
%   wactind: beta[actind]
%   distname: name of distribution 
%   t: current time point
%
% Output:
%   dbetadt: path updating direction

% t is dummy variable
xw=x(:, actind)*wactind;
beta0=glmfit(ones(size(x,1), 1), y, distname, 'constant', 'off', 'offset', xw);
eta=xw+beta0;

[ginv, ginvd, ginvdd, ginvddd]=invlinkDERIV_GLM(eta, distname);
u=ginv;
[Q1, Q11, Q111]=qDERIV_GLM(u, y, distname);

[dbeta0dbeta1, dbeta0dbeta2]=dbeta0dbeta_GLM(x(:, actind), ginvd, ginvdd, ginvddd, Q1, Q11, Q111);

deriv=dQdbeta_GLM(ginvd, Q1, x(:, actind), dbeta0dbeta1);
M=hessianQdbeta_GLM(ginvd, ginvdd, Q1, Q11, x(:, actind), dbeta0dbeta1, dbeta0dbeta2);

dbetadt=-inv(M)*diag(1./weight(actind))*(-sign(deriv));