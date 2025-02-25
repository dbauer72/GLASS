function [th,RN] = param_syst_aDFM_Utilde(param,r,n)
% param_syst_aDFM_Lambda converts a parameter vector into a normalized system using a fixed Lambda.
%
% SYNTAX:  [th,RN] = param_syst_aDFM_Utilde(param,N,r,n)
%
% INPUTS: param ... parameter vector. 
%         r     ... static factor dimension
%         n     ... state dimension in system for factor
%
% OUTPUTS:  th ... theta structure for system
%           RN ... rxr RN-part of the loading matrix. 
%
% REMARK: model: y_t = Lambda F_t + chi_t; 
%  F_t = Cx_t + Du_t, x_{t+1} = A x_t + Bu_t. 
% 
% Normalisation: + Lambda'Lambda/T = (I_r/N + RN*RN').
%                + [Ir,0] Lambda =I_r. 
%                + D=I_r.
%                + (A,B,C) in echelon form. 
% AUTHOR: dbauer, 24.2.2025

% param: theta: 2nr parameters for k(z), r(r+1)/2 for Omega, tau_R: r(r+1)/2 parameters. 

nparth = 2*n*r+ r*(r+1)/2;
theta = param(1:nparth); 
param(1:nparth)=[];
tauR= param; 

% parameters for C and A
C= zeros(r,n);
A= zeros(n,n);
if (r>n) % only if n is larger than r the matrix C contains parameters. 
    parc = theta(1:(n*(r-n)));
    theta(1:(n*(r-n))) = [];
    C = [eye(n);reshape(parc,r-n,n)]; 
    para = theta(1:n^2);
    theta(1:n^2)=[];
    A = reshape(para,n,n);
else
    C = [eye(r),zeros(r,n-r)];
    if (n>r)
        A(1:(n-r),(r+1):end) = eye(n-r);
    end
    para = theta(1:n*r);
    theta(1:n*r)=[];
    A((n-r+1):end,:)=reshape(para,r,n);
end

% parameters for B
parb = theta(1:n*r);
theta(1:n*r)=[];
B = reshape(parb,n,r);

% D=I_r. 
D = eye(r);

% Omega 
Omega = fill_lowtri(theta,r);


% build theta structure
th = theta_urs;
th.which = 'SS';
th.A = A;
th.C= C;
th.B = B;
th.D = D; 
th.Omega = Omega;

% RN
RN = zeros(r,r);
for i=1:r
    RN(i:end,i)=tauR(1:(r-i+1));
    tauR(1:(r-i+1))=[];
end

