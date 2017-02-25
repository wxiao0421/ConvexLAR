function [allw1, fullw1] = convexGLASSO(x1, y, model, distname, index1, method)

% Function:
%   Calculate the group LASSO solution path of general convex function
%
% Arguments:
%   x1: n x p covariate matrix
%   y: p x 1 response vector
%   model: 'glm', 'panel', 'adaBoost', 'recurrent'
%   distname: 'normal', 'poisson', 'binomial', specified when model is 'glm'
%   index1: the index set
%   method: 'grouplAR' (GroupConvexLASSO)
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
% set random seed 
% rand('state', 0);
s = RandStream('swb2712','Seed',0);
RandStream.setGlobalStream(s);
% randn('state', 0);
s = RandStream('mcg16807','Seed', 0);
RandStream.setGlobalStream(s);

%sort x1 and index1
[index IX]=sort(index1);
x=x1;
for i=1:size(x1,2)
    x(:,i)=x1(:,IX(i));
end
unsorted = 1:length(index1);
newIX(IX) = unsorted;

% calculate p
p=size(x,2);
% calculate count, countG, weight, weightG, deriv, derivG, actind, actindG,
% groupW, p2
INDEX=unique(index);
p2=length(INDEX);
count=zeros(p,1);
for i=1:p
    count(i)=sum(index==index(i));
end
countG=zeros(p2,1);
for i=1:p2
    countG(i)=sum(index==INDEX(i));
end
weight=1./sqrt(count);
weightG=1./sqrt(countG);
switch lower(model)
    case 'glm'      
    % calculate b(beta) with beta=0
    deriv=getderiv_GLM(x, y, 1, 0.0, distname);
    case {'panel', 'ada', 'recurrent'}
    deriv=getderiv_EXP(x, y, 1, 0.0, model);       
    otherwise
    disp('Unknown model.')
end

derivG=zeros(p2,1);
for i=1:p2
    derivG(i)=sqrt(sum((deriv(index==i)).^2));
end

groupW=zeros(p,1);
for i=1:p
    groupW(i)=deriv(i)/derivG(index(i));
end;

% calculate active group(actindG) and maxDeriv 
[maxDeriv, actindG]=max(abs(derivG).*weightG);
% calculate indexActive
indexActive=zeros(p,1);
for i=1:p
    for j=1:length(actindG)
        if index(i)==actindG(j)
            indexActive(i)=1;
        end
    end
end
% calculate active set(actind)
actind=zeros(sum(indexActive),1);
temp=1;
for i=1:p
    if indexActive(i)==1
        actind(temp)=i;
        temp=temp+1;
    end
end
%calculate inactindG
inactindG=setdiff(1:p2, actindG); 
% calculate current value of beta[actind] (wtemp)
w0=zeros(p, 1);
wtemp=w0(actind);

% start time point T0
T0=0;
Tmax=maxDeriv;
% save result
allw=zeros(2*p+1, 2);
fullw=zeros(2*p+1, 2); % (, 2) make fullw to be a matrix
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
        options=odeset(options,'Events', @(t, f)event4_GLM(t, f, x, y, ...
            distname, actind, actindG, p2, index, weightG, T0));
        [T,Y,TE,YE,IE]=ode45( @(t, f)mvDirectionG_GLM(x, y, actind, f, ...
            distname, actindG, weight, weightG, countG, groupW, index, ...
            method, t),[T0, Tmax], wtemp, options);
        case {'panel', 'ada', 'recurrent'}
        options=odeset(options,'Events', @(t, f)event4_EXP(t, f, x, y, ...
            model, actind, actindG, p2, index, weightG, T0));
        [T,Y,TE,YE,IE]=ode45( @(t, f)mvDirectionG_EXP(x, y, actind, f, ...
            model, actindG, weight, weightG, countG, groupW, index, method, ....
            t),[T0, Tmax], wtemp, options);
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
    
    % updata dertivG, groupW
    for i=1:p2
    derivG(i)=sqrt(sum((deriv(index==i)).^2));
    end
    groupW=zeros(p,1);
    for i=1:p
        groupW(i)=deriv(i)/derivG(index(i));
    end;

    %if abs(TE-Tmax)<1e-6 then terminate the loop
    if (abs(TE-Tmax)<1e-6), break, end
    %define Event
    Event=IE(end);
    if (Event==1) 
        [dump, ind]=max(abs(derivG(inactindG)));
        actindG=union(actindG,inactindG(ind));
        %updata actind
        indexActive=zeros(p,1);
        for i=1:p
            for j=1:length(actindG)
                if index(i)==actindG(j)
                    indexActive(i)=1;
                end
            end
        end
        % calculate active set(actind)
        actind=zeros(sum(indexActive),1);
        temp=1;
        for i=1:p
            if indexActive(i)==1
                actind(temp)=i;
                temp=temp+1;
            end
        end
        %calculate inactindG
        inactindG=setdiff(1:p2, actindG);
    elseif (Event==2)
        betaG=zeros(length(actindG),1); 
        for i=1:length(actindG);
            betaG(i)=sqrt(sum((YE(index(actind)==actindG(i))).^2));
        end;
        k=find(betaG<=1e-3);
        %undate actindG
        %delete group with l2 norm equal to 0
        actindG=setdiff(actindG,actindG(k));
        %updata actind
        indexActive=zeros(p,1);
        for i=1:p
            for j=1:length(actindG)
                if index(i)==actindG(j)
                    indexActive(i)=1;
                end
            end
        end
        % calculate active set(actind)
        actind=zeros(sum(indexActive),1);
        temp=1;
        for i=1:p
            if indexActive(i)==1
                actind(temp)=i;
                temp=temp+1;
            end
        end
        %calculate inactindG
        inactindG=setdiff(1:p2, actindG);
    end;
    % updata wulp
    wulp=size(fullw,2);
    % update start time point
    T0=TE;
    % update wtemp (beta[actind])
    wtemp=allw(actind, myloop);
%     % break when inactindG is empty
%     if size(inactindG,2)==0
%         break;
end;

% undo the sort
allw1=allw;
fullw1=fullw;
allw1(1:p, :)=allw(newIX, :);
allw1(p+1+(1:p), :)=allw(p+1+newIX, :);
fullw1(1:p, :)=fullw(newIX, :);
fullw1(p+1+(1:p), :)=fullw(p+1+newIX, :);
    