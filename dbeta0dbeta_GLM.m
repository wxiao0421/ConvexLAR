function [dbeta0dbeta1, dbeta0dbeta2]=dbeta0dbeta_GLM(x, ginvd, ginvdd, ginvddd, Q1, Q11, Q111)

% Function:
%   Calculate first and second order partial derivative of beta0 with
%   respect to beta
%
% Arguments:
%   u: beta0+ x%*%beta
%   y: response y (vector)
%   distname: name of distribution 
%
% Output:
%   dbeta0dbeta1= first order partial derivative of beta0 with respect to
%   beta
%   dbeta0dbeta2= second order partial derivative of beta0 with respect to
%   beta

% Details are given in section 3 (Details for deriving the path updating 
% direction) of "QuasiLARS: A Solution Path Algorithm for the 
% Quasi-likelihood Method"

temp=Q11.*(ginvd.^2)+Q1.*ginvdd;
dbeta0dbeta1=-x'*temp/(sum(temp));

n=size(x, 1);
x1=x+repmat(dbeta0dbeta1', n, 1);

temp1=Q111.*(ginvd.^3)+2*Q11.*ginvd.*ginvdd+Q11.*ginvd.*ginvdd+Q1.*ginvddd;
temp2=Q11.*(ginvd.^2)+Q1.*ginvdd;

dbeta0dbeta2=-x1'*diag(temp1)*x1/sum(temp2);
