function [value,isterminal,direction]=event1_GLM(t,f,x,y,actind,distname,weight)
% Function:
%   Define the event function used by ode45
%
% Arguments:
%   t, f
%   x: covariates x (matrix)
%   y: response y (vector)
%   actind: active index set
%   distname: name of distribution 
%
% Output:
%   [value,isterminal,direction]
p=size(x, 2);
inactind=setdiff(1:p, actind);    
deriv=getderiv_GLM(x, y, actind, f, distname);
% value is the value of the ith event function
value=min(abs(deriv(actind)).*weight(actind))- ...
    max(abs(deriv(inactind)).*weight(inactind)); 
% the integration is to terminate at a zero of this event function
isterminal = 1; 
% direction = 0 if all zeros are to be located
direction = 0;
 