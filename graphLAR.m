function [allw, fullw] = graphLAR(x)

% Function:
%   Calculate the ConvexLAR solution path of Gaussian graphical model
%
% Arguments:
%   x: p1 x p1 covariate matrix
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

global w0_INI
% specify relative and absolute error tolerance 
options=odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% set random seed 
% rand('state', 0);
s = RandStream('swb2712','Seed',0);
RandStream.setGlobalStream(s);
% randn('state', 0);
s = RandStream('mcg16807','Seed', 0);
RandStream.setGlobalStream(s);

% calculate p1, p
p1=size(x,2);
p=p1*(p1-1)/2;

% initial value of w0
w0_INI=(1./var(x))';

% calculate b(beta) with beta=0
[deriv, w0_INI]=getderiv_GRAPH(x, 1, 0.0);

% calculate active set and maxDeriv 
[maxDeriv, actind]=max(abs(deriv));
% calculate current value of beta[actind]
w0=zeros(p, 1);
wtemp=w0(actind);

% start time point T0
T0=-maxDeriv;

% save result
allw=zeros(2*p+1,p+1);
fullw=zeros(2*p+1, 2); % (, 2) make fullw to be matrix
wulp=1;
myloop=1;
allw(p+1, myloop)=T0;
fullw(p+1, wulp)=T0;
allw(p+1+(1:p), myloop)=deriv;
fullw(p+1+(1:p), wulp)=deriv;

while (myloop<=(p-1))
    %Keep the value of w0_INI in the begin step of w0_temp 
    w0_temp=w0_INI;

    % add 'Events' to old options
    options=odeset(options,'Events', @(t, f)event1_GRAPH(t, f, x, actind));
    [T,Y,TE,YE,IE]=ode45( @(t, f)mvDirection_GRAPH(x, actind, f, t),[T0, 0], wtemp, options);

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
    
    % calculate the deriv in intermidiate time points
	w0_INI=w0_temp;
    derivM=zeros(p,(npoint-1));
    for k=1:(npoint-1)
        [derivM(1:p,k), w0_INI]=getderiv_GRAPH(x, actind, Y(k+1,:));
    end    
    % update current deriv
    [deriv, w0_INI]=getderiv_GRAPH(x, actind, YE);

    % save result
    allw(p+1, myloop)=TE;
    allw(actind, myloop)=YE;
    allw(p+1+(1:p), myloop)=deriv;
    fullw(actind, (wulp+1):(wulp+npoint-1))=Y(2:npoint,:)';
    fullw(p+1, (wulp+1):(wulp+npoint-1))=T(2:npoint);
    fullw(p+1+(1:p), (wulp+1):(wulp+npoint-1))=derivM;
    % updata actind
    inactind=setdiff(1:p,actind);
    [dump, ind]=max(abs(deriv(inactind)));
    actind=union(actind,inactind(ind));
    %updata wulp
    wulp=size(fullw,2);
    % update start time point
    T0=TE;
    % update wtemp (beta[actind])
    wtemp=allw(actind, myloop);
    
end

%update w0_INI automatically in mvDirection to aviod singularity
%Keep the value of w0_INI in the begin step of w0_temp 
w0_temp=w0_INI;
% the last loop with all covariates in active set
options=odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[T,Y]=ode45( @(t, f)mvDirection_GRAPH(x, actind, f, t),[T0, 0], wtemp, options);
TE=T(end);
YE=Y(end,:);
npoint=length(T);
% calculate the deriv in intermidiate time points
w0_INI=w0_temp;
derivM=zeros(p,(npoint-1));
for k=1:(npoint-1)
    [derivM(1:p,k), w0_INI]=getderiv_GRAPH(x, actind, Y(k+1,:));
end
% update current deriv
[deriv, w0_INI]=getderiv_GRAPH(x, actind, YE);
% save result to fullw
fullw(actind, (wulp+1):(wulp+npoint-1))=Y(2:npoint,:)';
fullw(p+1, (wulp+1):(wulp+npoint-1))=T(2:npoint);
fullw(p+1+(1:p), (wulp+1):(wulp+npoint-1))=derivM; 
% updata myloop
myloop=myloop+1;
% save result
allw(p+1, myloop)=TE;
allw(actind, myloop)=YE;
allw(p+1+(1:p), myloop)=deriv; 