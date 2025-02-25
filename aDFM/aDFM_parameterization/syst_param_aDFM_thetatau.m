function param = syst_param_aDFM_thetatau(th,Lambda)
% syst_param_aDFM_thetatau converts a normalized system into a parameter vector.
%
% SYNTAX:  param = syst_param_aDFM_thetatau(th,Lambda)
%
% INPUTS:   th ... theta structure for system
%           Lambda ... Nxr loading matrix. 
%
% OUTPUTS: param ... parameter vector. 
%
%
% REMARK: model: y_t = Lambda F_t + chi_t; 
%  F_t = Cx_t + Du_t, x_{t+1} = A x_t + Bu_t. 
% 
% Normalisation: + Lambda'Lambda/T = I_r/N + RN*RN'.
%                + Lambda = [I;UN*RN*sqrt(N)], UN'*UN=I_r.  
%                + V(u_t) = Omega.
%                + D =I_r. 
%                + (A,B,C) in echelon form. 
% AUTHOR: dbauer, 24.2.2025

% ensure that norming is given
[th,RN,Utilde,Lambda] = norm_aDFM_Utilde(th,Lambda);

% extract parameters from system
param = [];
[n,r] = size(th.B);
% parameters for C
if (r>n) % only if r is larger than n the matrix C contains parameters. 
    pc = th.C((n+1):r,:);
    param = [param;pc(:)]; 
end

% parameters for A
if (n <= r)
    param = [param;th.A(:)];
else
    pa = th.A((n-r+1):end,:);
    param = [param;pa(:)];
end

% parameters for B
param = [param;th.B(:)];

% parameters for Omega
paraomi = extr_lowtri(th.Omega);
param = [param(:);paraomi(:)];

% if Lambda is also supplied 
if (nargin == 2)
    [N,r]=size(Lambda);

    % parameters for Lambda
    TrafoL = Lambda(1:r,1:r);
    LambdaL = Lambda*inv(TrafoL);
    indrev = [r:-1:1];
    [Q,R] = qr(LambdaL((r+1):N,indrev));
    UN = Q(:,indrev);
    RN = R(indrev,indrev)/sqrt(N);

    tauR = [];
    for i=1:r
        tauR = [tauR;RN(i:r,i)];
    end
    tauU = ortho2par(UN);

    param = [param;tauR;tauU(:)];
end




