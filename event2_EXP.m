function [value,isterminal,direction]=event2_EXP(t,f,x,y,actind,model,T0,weight)
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
switch lower(model)
    case {'panel', 'ada', 'recurrent'} 
    p=size(x, 2);
    inactind=setdiff(1:p, actind);    
    deriv=getderiv_EXP(x, y, actind, f, model);
    % value is the value of the ith event function
    % two events are used
    dt=t-T0;
    if (length(actind)<p)
        if (dt<1e-6)% Don't consider Event1 in the first 1e-6 interval
            minActind=min(abs(deriv(actind)).*weight(actind))+10;
            maxInactind=max(abs(deriv(inactind)).*weight(inactind));
        else
            minActind=min(abs(deriv(actind)).*weight(actind));
            maxInactind=max(abs(deriv(inactind)).*weight(inactind));
        end;
    else
        minActind=min(abs(deriv(actind)).*weight(actind)); maxInactind=-10;
        %When all variables are active, event 1 is not able to happen.
    end;
    value=[minActind-maxInactind;min(-sign(deriv(actind)).*f)]; 
    % the integration is to terminate at a zero of this event function
    isterminal = [1;1]; 
    % direction = 0 if all zeros are to be located
    % direction = -1 if only the zeros where the event function decreases
    direction = [0;-1];
    otherwise
    disp('Unknown model.')
end;