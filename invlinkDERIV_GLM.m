function [ginv, ginvd, ginvdd, ginvddd]=invlinkDERIV_GLM(eta, distname)

% Function:
%   Calculate ginv(eta[i]) and the first, second, third order derivative of
%   ginv(eta[i]) with respect to eta[i] for all i=1,...,n. Then stack them
%   into a vector.

%
% Arguments:
%   eta: beta0+ x%*%beta
%   distname: name of distribution 
%
% Output:
%   ginv, ginvd, ginvdd, ginvddd

switch lower(distname)
    case 'normal'
        % disp('Distribution is normal')
        % eta=g(u)=u
        ginv=eta;
        ginvd=ones(size(eta));
        ginvdd=zeros(size(eta));
        ginvddd=zeros(size(eta));
    case 'poisson'
        % disp('Distribution is poisson')
        % eta=g(u)=u
        ginv   =exp(eta);
        ginvd  =ginv;
        ginvdd =ginv;
        ginvddd=ginv;
    case 'binomial'
        % disp('Distribution is binomial')
        temp=exp(eta);
        tempdeno=1./(1+exp(eta));
        ginv    = temp.*tempdeno;
        ginvd   = ginv.*tempdeno;
        ginvdd  = temp.*(1-temp).*(tempdeno.^3);
        ginvddd = temp.*(1-4*temp+temp.^2).*(tempdeno.^4);
    otherwise
        disp('Unknown distribution.')
end
