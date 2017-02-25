function dbetadt=mvDirection_EXP(x, y, actind, wactind, model, weight, t)

% Function:
%   Define the event function used by ode45
%
% Arguments:
%   x: covariates x (matrix)
%   y: response y (vector)
%   actind: active index set
%   wactind: beta[actind]
%   distname: name of distribution 
%   t: current time point
%
% Output:
%   dbetadt: path updating direction

% t is dummy variable
switch lower(model)
    case 'panel'  
    % xw calculate x%*%beta as all inactind variable with coefficient equal to 0
    xactind=x(:, actind);
    xw=xactind*reshape(wactind, length(wactind), 1);
    temp1=y.*exp(-xw);
    deriv=-xactind'*temp1;
    M=xactind'*diag(temp1)*xactind;
    dbetadt=M\(diag(1./weight(actind))*(-sign(deriv)));
    
    case 'ada'
    % xw calculate x%*%beta as all inactind variable with coefficient equal to 0
    xw=x(:, actind)*reshape(wactind, length(wactind), 1);
    temp1=(y==-1)'*exp(xw);
    temp2=(y==1)'*exp(-xw);
    beta0=log(temp2/temp1)/2;
    temp3=x'*(exp(xw).*(y==-1));
    temp4=-x'*(exp(-xw).*(y==1));
    dbeta0dbeta1=(temp4/temp2-temp3/temp1)/2;
    deriv=dbeta0dbeta1*(exp(beta0)*temp1-exp(-beta0)*temp2)+...
        (exp(beta0)*temp3+exp(-beta0)*temp4); 
    deriv_act=deriv(actind);
    x1=x(y==-1,:);
    x2=x(y==1,:);
    temp5=x1'*diag(exp(xw(y==-1)))*x1;
    temp6=x2'*diag(exp(-xw(y==1)))*x2;
    dbeta0dbeta2=1/2*(temp6/temp2-temp4*temp4'/(temp2^2))-...
        1/2*(temp5/temp1-temp3*temp3'/(temp1^2));
    temp7=dbeta0dbeta1*(exp(beta0)*temp3-exp(-beta0)*temp4)';
    M=dbeta0dbeta2*(exp(beta0)*temp1-exp(-beta0)*temp2)...
        +(dbeta0dbeta1*dbeta0dbeta1')*(exp(beta0)*temp1+exp(-beta0)*temp2)...
        +temp7+temp7'+(exp(beta0)*temp5+exp(-beta0)*temp6);
    M_act=M(actind,actind);
    dbetadt=M_act\diag(1./weight(actind))*(-sign(deriv_act));

    case 'recurrent'
    % xw calculate x%*%beta as all inactind variable with coefficient equal to 0
    xw=x(:, actind)*reshape(wactind, length(wactind), 1);
    n=size(x,1);
    p=size(x,2);
    % failure time of each subject
    subFT=unique(y(:,[1 3]), 'rows');
    expV=exp(xw);
    ynew=y(y(:,2)>0,:); % delete row where failure time is 0;
    % set the initial value of deriv and M
    deriv=zeros(p,1);
    M=zeros(p,p);
    for i=1:size(ynew,1) 
        t=ynew(i,2);
        yt=subFT(:,2)>t;
        % temp0 is a n by 1 vector with ith component equal to yi(t)*exp(xi(t)beta)
        temp0=yt.*expV;
        temp1=yt'*expV;
        temp2=x'*temp0;
        temp3=x'*diag(temp0)*x;
        id=ynew(i,1);
        deriv=deriv+temp2/temp1-x(id,:)';
        M=M+temp3/temp1-temp2*temp2'/(temp1)^2;
    end;
    % standardize deriv and M by divide by n
    deriv=deriv/n;
    M=M/n;
    deriv_act=deriv(actind);  
    M_act=M(actind,actind);
    dbetadt=inv(M_act)*diag(1./weight(actind))*(-sign(deriv_act));
    
    otherwise
    disp('Unknown model.')
end;
