function dbetadt=mvDirectionG_GLM(x, y, actind, wactind, distname, actindG, weight, weightG, countG, groupW, index, method, t)

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
xw=x(:, actind)*wactind;
beta0=glmfit(ones(size(x,1), 1), y, distname, 'constant', 'off', 'offset', xw);
eta=xw+beta0;
[ginv, ginvd, ginvdd, ginvddd]=invlinkDERIV_GLM(eta, distname);
u=ginv;
[Q1, Q11, Q111]=qDERIV_GLM(u, y, distname);
[dbeta0dbeta1, dbeta0dbeta2]=dbeta0dbeta_GLM(x(:, actind), ginvd, ginvdd, ginvddd, Q1, Q11, Q111);
%deriv, M are first and second order derivatives of f(x),
%f(x)=-loglikelihood
deriv=-dQdbeta_GLM(ginvd, Q1, x(:, actind), dbeta0dbeta1);
M=-hessianQdbeta_GLM(ginvd, ginvdd, Q1, Q11, x(:, actind), dbeta0dbeta1, dbeta0dbeta2);
% construct fake derivF
derivF=zeros(size(x,2),1);
derivF(actind)=deriv;
% construct fake wactindF
wactindF=zeros(size(x,2),1);
wactindF(actind)=wactind;
% construct fake M
MF=zeros(size(x,2),size(x,2));
MF(actind,actind)=M;

switch lower(method)
    case {'grouplarl2', 'grouplarl1'}
    dbetadt=-inv(M)*diag(1./weight(actind))*groupW(actind);
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
    M_act_D=M+D;
    if isempty(actindG0)
        %[e,lam]=eig(M_act_D);
        %fprintf('matrix inverse of symmetrix matrix using spectral\n');
        %invMD=e*diag(1./diag(lam))*e';
        %inv(M_act+D)
        %dbetadt=-invMD*deriv./st;
        dbetadt=-M_act_D\(deriv./st);
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
            dbetadt=BM\deriv;
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