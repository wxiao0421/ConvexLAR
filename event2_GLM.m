function [value,isterminal,direction]=event2_GLM(t,f,x,y,actind,distname,T0,weight)
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
% two events are used
if (length(actind)<p) 
    minActind=min(abs(deriv(actind)).*weight(actind)); maxInactind=max(abs(deriv(inactind)).*weight(inactind));
else
    minActind=min(abs(deriv(actind)).*weight(actind)); maxInactind=-10;%event 1 is not able to happen
end;
value=[minActind-maxInactind;min(-sign(deriv(actind)).*f)]; 
% the integration is to terminate at a zero of this event function
isterminal = [1;1]; 
%direction(i) = -1 if only zeros where the event function is decreasing
direction = [-1;-1];
 