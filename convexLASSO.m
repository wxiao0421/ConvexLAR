function [allw, fullw] = convexLASSO(x, y, model, distname, weight)
% Function:
%   Calculate the ConvexLASSO (weighted ConvexLASSO) solution path for general convex function
%
% Arguments:
%   x: n x p covariate matrix
%   y: p x 1 response vector
%   model: 'glm', 'panel', 'adaBoost'
%   distname: 'normal', 'poisson', 'binomial', specified when model is 'glm'
%   weight: p x 1 vector of weight corresponding to each covariate, where a
%           larger weight wj inplies higher "influence" 
%
% Output:
%   allw: solution and first order derivatives at all transition points,
%         where 1:p rows is the solution at grid points, p+1 row corresponds the 
%         time points, (p+2):(2p+1) rows is the first order derivatives at grid 
%         points
%   fullw: solution and first order derivatives at all points
%
% COPYRIGHT: North Carolina State University
% AUTHOR: Wei Xiao, wxiao@ncsu.edu; Yichao Wu, ywu11@ncsu.edu 


% specify relative and absolute error tolerance 
options=odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
myeps3=1e-7;
% set random seed 
% rand('state', 0);
s = RandStream('swb2712','Seed',0);
RandStream.setGlobalStream(s);
% randn('state', 0);
s = RandStream('mcg16807','Seed', 0);
RandStream.setGlobalStream(s);

% calculate p
p=size(x,2);

switch lower(model)
    case 'glm'      
    % calculate b(beta) with beta=0
    deriv=getderiv_GLM(x, y, 1, 0.0, distname);
    case {'panel', 'ada', 'recurrent'}
    deriv=getderiv_EXP(x, y, 1, 0.0, model);  
    otherwise
    disp('Unknown model.')
end

% calculate active set and maxDeriv 
[maxDeriv, actind]=max(abs(deriv).*weight);
% calculate current value of beta[actind]
w0=zeros(p, 1);
wtemp=w0(actind);

% calculate Tmax
Tmax=maxDeriv;
T0=0;

% save result
allw=zeros(2*p+1, 2);
fullw=zeros(2*p+1, 2); % (, 2) make fullw to be matrix
wulp=1;
myloop=1;
allw(p+1, myloop)=T0;
fullw(p+1, wulp)=T0;
allw(p+1+(1:p), myloop)=deriv;
fullw(p+1+(1:p), wulp)=deriv;

while (1)
    % add 'Events' to old options
    switch lower(model)
        case 'glm' 
        options=odeset(options,'Events', @(t, f)event2_GLM(t, f, x, y, ...
            actind, distname, T0, weight));
        [T,Y,TE,YE,IE]=ode45( @(t, f)mvDirection_GLM(x, y, actind, f, ...
            distname, weight, t),[T0, Tmax], wtemp, options);
        case {'panel', 'ada', 'recurrent'}
        options=odeset(options,'Events', @(t, f)event2_EXP(t, f, x, y, ...
            actind, model,T0, weight));
        [T,Y,TE,YE,IE]=ode45( @(t, f)mvDirection_EXP(x, y, actind, f, ...
            model, weight, t),[T0, Tmax], wtemp, options);
        otherwise
        disp('Unknown model.')
    end
    
    % TE, YE will fail when length(actind)=p (last iteration), so calculate
    % them from T, Y
    TE=T(end);
    YE=Y(end,:);
    if (length(actind)==1)
        %keep the num of col of Y to be equal to length(actind)
        Y=vec2mat(Y,length(Y))';
    end
    % updata myloop
    myloop=myloop+1;
    npoint=length(T);
    
    switch lower(model)
        case 'glm' 
        % update current deriv
        deriv=getderiv_GLM(x, y, actind, YE, distname);
        % calculate the deriv in intermidiate time points
        derivM=zeros(p,(npoint-1));
        for k=1:(npoint-1)
            derivM(1:p,k)=getderiv_GLM(x, y, actind, Y(k+1,:), distname);
        end
        case {'panel', 'ada', 'recurrent'}
        % update current deriv
        deriv=getderiv_EXP(x, y, actind, YE, model);
        % calculate the deriv in intermidiate time points
        derivM=zeros(p,(npoint-1));
        for k=1:(npoint-1)
            derivM(1:p,k)=getderiv_EXP(x, y, actind, Y(k+1,:), model);
        end            
        otherwise
        disp('Unknown model.')
    end
    
    % save result
    allw(p+1, myloop)=TE;
    allw(actind, myloop)=YE;
    allw(p+1+(1:p), myloop)=deriv;
    fullw(actind, (wulp+1):(wulp+npoint-1))=Y(2:npoint,:)';
    fullw(p+1, (wulp+1):(wulp+npoint-1))=T(2:npoint);
    fullw(p+1+(1:p), (wulp+1):(wulp+npoint-1))=derivM;
    % if abs(TE-Tmax)<1e-6 then terminate the loop
    if (abs(TE-Tmax)<1e-6), break, end
    % define Event
    Event=IE(end);
    if (Event==1) 
        % updata actind
        inactind=setdiff(1:p,actind);
        [dump, ind]=max(abs(deriv(inactind)).*weight(inactind));
        actind=union(actind,inactind(ind));
    elseif (Event==2)
        k=find(abs(YE(:))<myeps3);
        % undate actind
        % delete variable equal to 0
        actind=setdiff(actind,actind(k));
    end;

    
    % updata wulp
    wulp=size(fullw,2);
    % update start time point
    T0=TE;
    % update wtemp (beta[actind])
    wtemp=allw(actind, myloop);
end;   
    