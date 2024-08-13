function [th,Omega] = param_syst_aDFM_Lambda(param,N,r,q,n)
% param_syst_aDFM_LAmbda converts a parameter vector into a normalized system using a fixed Lambda.
%
% SYNTAX:  param = syst_param_aDFM(th,Lambda)
%
% INPUTS: param ... parameter vector. 
%         N     ... cross sectional dimension.
%         r     ... static factor dimension
%         q     ... dynamic factor dimension
%         n     ... state dimension in system for factor
%
% OUTPUTS:   th ... theta structure for system
%           Lambda ... Nxr loading matrix. 
%
% REMARK: model: y_t = Lambda F_t + chi_t; 
%  F_t = Cx_t + Du_t, x_{t+1} = A x_t + Bu_t. 
% 
% Normalisation: + Lambda'Lambda/T = I_r.
%                + Lambda p.l.t. (heading square matrix lower triangular
%                  with positive elements on diagonal. 
%                + V(u_t)=I_q. 
%                + D is p.l.t. 
%                + (A,B,C) in echelon form. 
% AUTHOR: dbauer, 24.7.2024

%nparLam = (N-r)*r; 
%Lambda = par2ortho_plt(param(1:nparLam),N,r)*sqrt(N);
% now extract parameters to form system
%param(1:nparLam) = [];

% parameters for C and A
C= zeros(r,n);
A= zeros(n,n);
if (r>n) % only if n is larger than r the matrix C contains parameters. 
    parc = param(1:(n*(r-n)));
    param(1:(n*(r-n))) = [];
    C = [eye(n);reshape(parc,r-n,n)]; 
    para = param(1:n^2);
    param(1:n^2)=[];
    A = reshape(para,n,n);
else
    C = [eye(r)];
    if (n>r)
        A(1:(n-r),(r+1):end) = eye(n-r);
    end
    para = param(1:n*r);
    param(1:n*r)=[];
    A((n-r+1):end,:)=reshape(para,r,n);
end

% parameters for B
parb = param(1:n*q);
param(1:n*q)=[];
B = reshape(parb,n,q);

% parameters for D
D = zeros(r,q);
for j=1:q
    pard = param(1:(r-j+1));
    D(j:end,j) = pard;
    param(1:(r-j+1)) = [];
end

th = theta_urs;
th.which = 'SS';
th.A = A;
th.C= C;
th.B = B;
th.D = D; 
th.Omega = eye(q);

% parameters for Omega
if ~isempty(param)
    Omega = fill_lowtri(param,r);
end
