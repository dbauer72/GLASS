function [th,RN] = param_syst_aDFM_Utilde(param,r,n,c)
% param_syst_aDFM_Lambda converts a parameter vector into a normalized system using a fixed Lambda.
%
% SYNTAX:  [th,RN] = param_syst_aDFM_Utilde(param,N,r,n)
%
% INPUTS: param ... parameter vector. 
%         r     ... static factor dimension
%         n     ... state dimension in system for factor
%         c     ... integer; number of common trends in common factors
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
% AUTHOR: dbauer, 10.3.2025

if nargin<3
    c=0; % stationary case
end
% param: theta: 2(n-c)r parameters for k(z), r(r+1)/2 for Omega, tau_R: r(r+1)/2 parameters. 
% and 2*r*c-c^2 for common trends part. 

if (c>0)
    parc = param(1:(2*r*c-c^2));
    param(1:(2*r*c-c^2)) = [];
    nparc = r*c - c*(c+1)/2; 
    C1 = par2ortho(parc(1:nparc),r,c);
    parc(1:nparc)=[];

    B1 = zeros(c,r);
    for jj=1:c
        B1(jj,jj:end)=parc(1:(r-jj+1));
        parc(1:(r-jj+1))=[];
    end
end

nparth = 2*(n-c)*r+r*(r+1)/2; % 2*n*r+ r*(r+1)/2;
theta = param(1:nparth); 
param(1:nparth)=[];
tauR= param; 

% parameters for C and A
Cbull= zeros(r,n-c);
Abull= zeros(n-c,n-c);
nc = n-c; 

if (r>nc) % only if n is larger than r the matrix C contains parameters. 
    parc = theta(1:(nc*(r-nc)));
    theta(1:(nc*(r-nc))) = [];
    Cbull = [eye(nc);reshape(parc,r-nc,nc)]; 
    para = theta(1:nc^2);
    theta(1:nc^2)=[];
    Abull = reshape(para,nc,nc);
else
    Cbull = [eye(r),zeros(r,nc-r)];
    if (n>r)
        Abull(1:(nc-r),(r+1):end) = eye(nc-r);
    end
    para = theta(1:nc*r);
    theta(1:nc*r)=[];
    Abull((nc-r+1):end,:)=reshape(para,r,nc);
end

% parameters for B
parb = theta(1:nc*r);
theta(1:nc*r)=[];
Bbull = reshape(parb,nc,r);

% fill in two parts of system 
C = zeros(r,n);
A= zeros(n,n);
B = zeros(n,r); 

if (c>0) 
    C(:,1:c) = C1;
    C(:,(c+1):n)=Cbull; 
    A(1:c,1:c)=eye(c);
    A((c+1):end,(c+1):end)=Abull;
    B(1:c,:)=B1;
    B((c+1):end,:)=Bbull;
else
    A= Abull;
    B= Bbull;
    C= Cbull; 
end

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

