function dbetadt=mvDirection_GRAPH(x, actind, wactind, t)

% Function:
%   Define the event function used by ode45
%
% Arguments:
%   x: covariates x (matrix)
%   actind: active index set
%   wactind: beta[actind]
%   t: current time point
%   w0_INI: initial value of w0
%
% Output:
%   dbetadt: path updating direction

% t is dummy variable
global w0_INI

p1=size(x,2);
p=p1*(p1-1)/2;

w1=zeros(p,1);
w1(actind)=wactind;
w1_M=zeros(p1,p1);

k=1;
for i=2:p1
    for j=1:(i-1)
        w1_M(j, i)=w1(k);
        w1_M(i, j)=w1(k);
        k=k+1;
    end;
end;

SigmaHat=cov(x);
%calculate temp1, temp2,temp3
temp1=zeros(p1*p1,p1);
for i=1:p1
    temp1(1+(i-1)*(p1+1),i)=1;
end;

temp2=zeros(p1*p1,p1*(p1-1)/2);
k=1;
for i=2:p1
    for j=1:(i-1)
        M=zeros(p1,p1);
        M(i,j)=1;
        M(j,i)=1;
        temp2(:,k)=M(:); % M(:)=vect(M)
        k=k+1;
    end;
end;

k=1;
w0=w0_INI; 
% calculate w0(w1) using Newton's method
while 1
    %calculate Omega
    Omega=w1_M;
    for j=1:p1
        Omega(j,j)=w0(j);
    end;    
    invOmega=inv(Omega);
    temp4=inv(temp1'*kron(invOmega,invOmega)*temp1);
    w0_update=w0+temp4*diag(invOmega-SigmaHat);
    k=k+1;
    if (sum(abs(w0_update-w0))<1e-8)
        break;
    end;
    w0=w0_update;
end

% calculate deriv
deriv=zeros(p,1);
k=1;
for i=2:p1
    for j=1:(i-1)
        diff=invOmega-SigmaHat;
        deriv(k,1)=2*(diff(j,i));
        k=k+1;
    end;
end;

% calculate Dw0(w1)
temp3=kron(invOmega,invOmega);
temp4=inv(temp1'*temp3*temp1);
temp5=temp1'*temp3*temp2;
Dw0w1=-temp4*temp5;

% calculate M
M=-temp2'*temp3*(temp2+temp1*Dw0w1);

% calculate dbetadt
deriv_act=deriv(actind);  
M_act=M(actind,actind);
dbetadt=inv(M_act)*(-sign(deriv_act));

%update w0_INI
w0_INI=w0;

