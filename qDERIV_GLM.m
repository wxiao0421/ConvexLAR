function [Q1, Q11, Q111]=qDERIV_GLM(u, y, distname)

% Function:
%   Calculate the first, second, third order partial derivative of
%   quasilikelihood function Q(u_{i}, y_{i}) (u_{i}=beta0+ beta%*%x_{i})
%   with respect to u_{i} for all i=1,...,n. Then stack them into vector
%   Q1, Q11, Q111
%
% Arguments:
%   u: beta0+ x%*%beta
%   y: response y (vector)
%   distname: name of distribution 
%
% Output:
%   Q1, Q11, Q111

switch lower(distname)
    case 'normal'
        % disp('Distribution is normal')
        % Q(u,y)=-(y-u)^2/2
        Q1=-(u-y);
        Q11=-ones(size(y));
        Q111=zeros(size(y));
    case 'poisson'
        % disp('Distribution is poisson')
        % Q(u,y)=-(y-u)^2/2
        temp=1./u;
        Q1=y.*temp-1;
        Q11=-y.*(temp.^2);
        Q111=2.0*y.*(temp.^3);
    case 'binomial'
        % disp('Distribution is binomial')
        % Q(u,y)=ylog[u/(1-u)]+log(1-u)        
        Q1    =(u-y)./(u.^2-u);
        Q11   =-(u.^2+y-2*u.*y)./((u.^2).*((u-1).^2));
        Q111  =2*(u.^3-y+3*u.*y-3*(u.^2).*y)./((u.^3).*((u-1).^3));

        % temp1=1./u;
        % temp2=1./(1-u);
        % Q1   =y.*temp1+(1-y).*temp2;
        % Q11  =-y.*(temp1.^2)+(1-y).*(temp2.^2);
        % Q111 =2*y.*(temp1.^3)+2*(1-y).*(temp2.^3);
       
    otherwise
        disp('Unknown distribution.')
end
