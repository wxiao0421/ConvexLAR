function [value,isterminal,direction]=event4_GLM(t,f,x,y,distname,actind,actindG,p2,index,weightG,T0)
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
inactindG=setdiff(1:p2, actindG);    
deriv=getderiv_GLM(x, y, actind, f, distname);
derivG=zeros(p2,1);
betaG=zeros(length(actindG),1);
indexACT=index(actind);
% two events are used
dt=t-T0;
for i=1:p2
    derivG(i)=sqrt(sum((deriv(index==i)).^2));
end
%     if (dt<1e-1)
%         for i=1:length(actindG);
%             betaG(i)=10;
%         end;
%     else
%         for i=1:length(actindG);
%             betaG(i)=sqrt(sum((f(indexACT==actindG(i))).^2));
%         end;
%     end;
for i=1:length(actindG);
    betaG(i)=sqrt(sum((f(indexACT==actindG(i))).^2));
end;
% value is the value of the ith event function
value=[min(derivG(actindG).*weightG(actindG))- ... 
    max(derivG(inactindG).*weightG(inactindG)); min(betaG)-2e-4];  
% the integration is to terminate at a zero of this event function
isterminal = [1;1]; 
% direction = 0 if all zeros are to be located
direction = [-1;-1];


 