function dbetadt=mvDirectionG_EXP(x, y, actind, wactind, model, actindG, weight, weightG, countG, groupW, index, method, t)

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
n=size(x,1);
p=size(x,2);  
deriv_act=zeros(p,1);
M_act=zeros(p,p);
    
switch lower(model)
    case 'panel' 
    % xw calculate x%*%beta as all inactind variable with coefficient equal to 0
    xw=x(:, actind)*reshape(wactind, length(wactind), 1);
    temp1=y.*exp(-xw);
    deriv_act=-(x(:,actind))'*temp1;
    for i=1:n
        M_act=M_act+temp1(i)*(x(i,actind)'*x(i,actind));
    end
    
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
    temp5=zeros(p,p);    
    for i=find(y==-1)' %find() is a column vector
        temp5=temp5+x(i,:)'*x(i,:)*exp(xw(i));
    end
    temp6=zeros(p,p);
    for i=find(y==1)'
        temp6=temp6+x(i,:)'*x(i,:)*exp(-xw(i));
    end 
    dbeta0dbeta2=1/2*(temp6/temp2-temp4*temp4'/(temp2^2))-...
        1/2*(temp5/temp1-temp3*temp3'/(temp1^2));
    temp7=dbeta0dbeta1*(exp(beta0)*temp3-exp(-beta0)*temp4)';
    M=dbeta0dbeta2*(exp(beta0)*temp1-exp(-beta0)*temp2)...
        +(dbeta0dbeta1*dbeta0dbeta1')*(exp(beta0)*temp1+exp(-beta0)*temp2)...
        +temp7+temp7'+(exp(beta0)*temp5+exp(-beta0)*temp6);
    M_act=M(actind,actind);

    case 'recurrent'
    % xw calculate x%*%beta as all inactind variable with coefficient equal to 0
    xw=x(:, actind)*reshape(wactind, length(wactind), 1);
    subFT=unique(y(:,[1 3]), 'rows');
    expV=exp(xw);
    ynew=y(y(:,2)>0,:); % delete row where failure time is 0;
    % set the initial value of deriv and M
    deriv=zeros(p,1);
    M=zeros(p,p);
    for i=1:size(ynew,1) 
        t=ynew(i,2);
        temp1=(subFT(:,2)>t)'*expV;
        temp2=x'*(expV.*(subFT(:,2)>t));
        % set the initial value of temp3,
        temp3=zeros(p,p); 
        for j=find(subFT(:,2)>t)'
            temp3=temp3+x(j,:)'*x(j,:)*expV(j);
        end;
        id=ynew(i);
        deriv=deriv+temp2/temp1-x(id,:)';
        M=M+temp3/temp1-temp2*temp2'/(temp1)^2;
    end;
    % standardize deriv and M by divide by n
    deriv=deriv/n;
    M=M/n;
    deriv_act=deriv(actind);  
    M_act=M(actind,actind);
    
    otherwise
    disp('Unknown model.')
end;
% construct fake derivF
derivF=zeros(size(x,2),1);
derivF(actind)=deriv_act;
% construct fake wactindF
wactindF=zeros(size(x,2),1);
wactindF(actind)=wactind;
% construct fake M
MF=zeros(size(x,2),size(x,2));
MF(actind,actind)=M_act;

switch lower(method)
    case {'grouplarl2', 'grouplarl1'}
    dbetadt=-inv(M_act)*diag(1./weight(actind))*groupW(actind);
    case {'grouplar'}
    p3=length(actindG);
    actindG0=[];
    actindG1=actindG;
	%calculate derivG
    derivG=zeros(p3,1);
    for i=1:p3
        num_g=actindG(i);
        derivG(i)=sqrt(sum((derivF(index==num_g)).^2));
        betaV=wactindF(index==num_g);
        if max(abs(betaV))<1e-30
            actindG0=actindG(i);
            actindG1=setdiff(actindG,actindG0);
        end;
    end
    %calculate s(t)
    st=derivG(1)*weightG(actindG(1));
    %calculate block matrices
    MT = cell(1, p3);
    for i = 1:p3
        num_g=actindG(i);
        derivg=derivF(index==num_g);
        temp = eye(countG(num_g))-(derivg./derivG(i))*(derivg./derivG(i))';
        betaV=wactindF(index==num_g);
         if max(abs(betaV))<1e-30 
             betaV=-derivg/sqrt(sum((derivg).^2))*1e-30;
         end;
        MT{i}=temp.*(st*sqrt(countG(num_g))/sqrt(sum((betaV).^2)));
    end
    D=blkdiag(MT{:});
    M_act_D=M_act+D;
    if isempty(actindG0)
        %[e,lam]=eig(M_act_D);
        %fprintf('matrix inverse of symmetrix matrix using spectral\n');
        %invMD=e*diag(1./diag(lam))*e';
        %inv(M_act+D)
        %dbetadt=-invMD*deriv./st;
        dbetadt=-M_act_D\(deriv_act./st);
    end;
    if ~isempty(actindG0) && length(actindG0)==1
        %use derived formula
        deriv_g1=derivF(index==actindG0(1));
        M_g1g1=MF(index==actindG0(1),index==actindG0(1));
        temp1=deriv_g1'*M_g1g1*deriv_g1;
        temp4=st*countG(actindG0(1));
        if isempty(actindG1)
            k=temp4/temp1;
            dbetadt=-k*derivF(index==actindG0(1));
        else
            deriv_G1=derivF(ismember(index,actindG1));
            M_g1G1=MF(index==actindG0(1),ismember(index,actindG1));
            MD_G1G1=M_act_D(ismember(index(actind),actindG1),ismember(index(actind),actindG1));
            MD_G1G0=M_act_D(ismember(index(actind),actindG1),ismember(index(actind),actindG0));
            temp2=deriv_g1'*M_g1G1*inv(MD_G1G1)*MD_G1G0*deriv_g1;
            temp3=-1/st*deriv_g1'*M_g1G1*inv(MD_G1G1)*deriv_G1;
            k=-(temp3+temp4)/(temp2-temp1);
            BM=zeros(length(actind),length(actind));
            BM(index(actind)==actindG0(1),index(actind)==actindG0(1))=-1/k*eye(countG(actindG0(1)));
            BM(ismember(index(actind),actindG1),:)=-1*st*M_act_D(ismember(index(actind),actindG1),:);
            dbetadt=BM\deriv_act;
            % check equation (24)
            % deriv_g1'*M(index(actind)==actindG0(1),:)*dbetadt+temp4
        end;
%             %use derived formula with another equation
%             m=1;
%             deriv_g1=derivF(index==actindG0(1));
%             deriv_gm=derivF(index==actindG1(m));
%             M_gmG0=MF(index==actindG1(m),index==actindG0(1));
%             temp1=deriv_gm'*M_gmG0*deriv_g1;
%             temp4=st*countG(actindG1(m));
%             deriv_G1=derivF(ismember(index,actindG1));
%             M_gmG1=MF(index==actindG1(m),ismember(index,actindG1));
%             MD_G1G1=M_act_D(ismember(index(actind),actindG1),ismember(index(actind),actindG1));	
%             MD_G1G0=M_act_D(ismember(index(actind),actindG1),ismember(index(actind),actindG0));
%             temp2=deriv_gm'*M_gmG1*inv(MD_G1G1)*MD_G1G0*deriv_g1;
%             temp3=-1/st*deriv_gm'*M_gmG1*inv(MD_G1G1)*deriv_G1;
%             k=-(temp3+temp4)/(temp2-temp1);
%             BM=zeros(length(actind),length(actind));
%             BM(index(actind)==actindG0(1),index(actind)==actindG0(1))=-1/k*eye(countG(actindG0(1)));
%             BM(ismember(index(actind),actindG1),:)=-1*st*M_act_D(ismember(index(actind),actindG1),:);
%             dbetadt=inv(BM)*deriv;
%             % check equation (24)
%             deriv_gm'*M(index(actind)==actindG1(m),:)*dbetadt+temp4
    end;        
    otherwise
    disp('Unknown method.')
end;
