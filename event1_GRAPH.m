function [value,isterminal,direction]=event1_GRAPH(t, f, x, actind)
% Function:
%   Define the event function used by ode45
%
% Arguments:
%   t, f
%   x: covariates x (matrix)
%   actind: active index set
%
% Output:
%   [value,isterminal,direction]
global w0_INI

p1=size(x, 2);
p=p1*(p1-1)/2;
inactind=setdiff(1:p, actind);    
[deriv,w0_INI]=getderiv_GRAPH(x, actind, f);
% value is the value of the ith event function
value=min(abs(deriv(actind)))-max(abs(deriv(inactind))); 
% the integration is to terminate at a zero of this event function
isterminal = 1; 
% direction = 0 if all zeros are to be located
direction = 0;