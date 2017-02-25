function [value,isterminal,direction]=event2_GRAPH(t, f, x, actind,T0)
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
% two events are used
dt=t-T0;
if (length(actind)<p)
    if (dt<1e-6)% Don't consider Event1 in the first 1e-6 interval
        minActind=min(abs(deriv(actind)))+10;maxInactind=max(abs(deriv(inactind)));
    else
        minActind=min(abs(deriv(actind)));maxInactind=max(abs(deriv(inactind)));
    end;
else
    minActind=min(abs(deriv(actind))); maxInactind=-10;
    %When all variables are active, event 1 is not able to happen.
end;
value=[minActind-maxInactind;min(-sign(deriv(actind)).*f)]; 

%the integration is to terminate at a zero of this event function
isterminal = [1;1]; 
%direction = 0 if all zeros are to be located
%direction = -1 if only the zeros where the event function decreases
direction = [0;-1];
