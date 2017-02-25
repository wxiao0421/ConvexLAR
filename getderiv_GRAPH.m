function [deriv, w0]=getderiv_GRAPH(x, actind, wactind)
% Function:
%   Calculate the first-order partial derivative (b(beta)) of f(beta) for
%   graphical model
%
% Arguments:
%   x: covariates x (matrix)
%   actind: active index set
%   wactind: current value of active index set
%   w0_INI: initial value of w0
%
% Output:
%   deriv: first-order partial derivative at current beta value

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
temp1=zeros(p1*p1,p1);
for i=1:p1
    temp1(1+(i-1)*(p1+1),i)=1;
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
        diff=-invOmega+SigmaHat;
        deriv(k,1)=2*(diff(j,i));
        k=k+1;
    end;
end;

