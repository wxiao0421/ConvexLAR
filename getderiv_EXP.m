function deriv=getderiv_EXP(x, y, actind, wactind, model)
% Function:
%   Calculate the first-order partial derivative (b(beta)) of f(beta) for
%   model with explicit b(beta). ('panel', 'ada', 'recurrent')
%   all values of non-active index set are 0
%
% Arguments:
%   x: covariates x (matrix)
%   y: response y (vector)
%   actind: active index set
%   wactind: current value of active index set
%   distname: name of distribution 
%
% Output:
%   deriv: first-order partial derivative at current beta value

switch lower(model)
    case 'panel' 
    % xw calculate x%*%beta as all inactind variable with coefficient equal to 0
    xw=x(:, actind)*reshape(wactind, length(wactind), 1);
    temp1=y.*exp(-xw);
    deriv=-x'*temp1;
    
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
    
    case 'recurrent'
    % xw calculate x%*%beta as all inactind variable with coefficient equal to 0
    xw=x(:, actind)*reshape(wactind, length(wactind), 1);
    n=size(x,1);
    p=size(x,2);  
    subFT=unique(y(:,[1 3]), 'rows'); 
    expV=exp(xw); %calculate X_j%*%beta
    ynew=y(y(:,2)>0,:); % delete row where failure time is 0;
    % set the initial value of deriv and M
    deriv=zeros(p,1);
    for i=1:size(ynew,1) 
        t=ynew(i,2);
        temp1=(subFT(:,2)>t)'*expV;
        temp2=x'*(expV.*(subFT(:,2)>t));
        id=ynew(i,1);
        deriv=deriv+temp2/temp1-x(id,:)';
    end;
    % standardize deriv and M by divide by n
    deriv=deriv/n;
    
    otherwise
    disp('Unknown model.')
end;



