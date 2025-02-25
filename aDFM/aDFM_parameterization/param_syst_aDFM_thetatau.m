function [th,Lambda,RN,UN,Rtilde,Utilde] = param_syst_aDFM_thetatau(param,N,r,n)
% param_syst_aDFM_thetatau converts a parameter vector into a normalized system.
%
% SYNTAX:  [th,Rtilde,Utilde] = param_syst_aDFM_thetatau(param,N,r,n)
%
% INPUTS: param ... parameter vector. 
%         N     ... cross sectional dimension.
%         r     ... static factor dimension
%         n     ... state dimension in system for factor
%
% OUTPUTS:   th ... theta structure for system
%            Lambda ... Nxr loading matrix. 
%            RN ... rxr p.l.t. matrix of loadings
%            UN ... (N-r)xr orthonormal part parameterizing part of the column space of  loading matrix. 
%            Rtilde ... rxr p.l.t. matrix of loadings
%            Utilde ... Nxr orthonormal part parameterizing column space of  loading matrix. 
%
% REMARK: model: y_t = Lambda F_t + chi_t; 
%  F_t = Cx_t + Du_t, x_{t+1} = A x_t + Bu_t. 
% 
% Normalisation: + Lambda'Lambda/T = I_r/N + RN*RN'.
%                + Lambda = [I;UN*RN*sqrt(N)]. 
%                + V(u_t) = Omega.
%                + D =I_r. 
%                + (A,B,C) in echelon form. 
% AUTHOR: dbauer, 24.2.2025

nparth = 2*n*r+r*(r+1)/2; 
npartau = (N-r)*r; 
theta = param(1:nparth);
tau = param((nparth+1):(nparth+npartau));
ntaur = r*(r+1)/2; 
tauR = tau(1:ntaur);
tau(1:ntaur)=[];
tauU = tau; 

% generate UN
UN = par2ortho(tauU,N-r,r)*sqrt(N);

% generate th and RN
[th,RN] = param_syst_aDFM_Utilde([theta(:);tauR(:)],r,n); 

% generate Rtilde, Utilde, Lambda
Lambda = zeros(N,r);
Lambda(1:r,1:r) =eye(r);
Lambda((r+1):end,:) = UN*RN; 

indrev = [r:-1:1];
[Q,R] = qr(Lambda(:,indrev));
Utilde = Q(:,indrev); %*sqrt(N);
Rtilde = R(indrev,indrev)/sqrt(N);


