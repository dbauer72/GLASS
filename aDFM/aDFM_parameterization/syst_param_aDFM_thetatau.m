function param = syst_param_aDFM_thetatau(th,Lambda,c)
% syst_param_aDFM_thetatau converts a normalized system into a parameter vector.
%
% SYNTAX:  param = syst_param_aDFM_thetatau(th,Lambda,c)
%
% INPUTS:   th ... theta structure for system
%           Lambda ... Nxr loading matrix. 
%           c ... integer; number of common trends in common factors
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
% AUTHOR: dbauer, 10.3.2025

if nargin<3
    c=0; % stationary case
end; 
% ensure that norming is given
[th,RN,Utilde,Lambda] = norm_aDFM_Utilde(th,Lambda,c);

% extract parameters from system
param = [];
[n,r] = size(th.B);

% for I(1): parameters in C_1,B_1. 
if (c>0)
    param = ortho2par(th.C(:,1:c));
    parb = [];
    for jj=1:c
        parb = [parb(:);th.B(jj,jj:end)];
    end
    param = [param(:)';parb(:)];

    indc = 1:c;
    indbull = (c+1):n; 

    Cbull = th.C(:,indbull);
    Abull = th.A(indbull,indbull);
    Bbull = th.B(indbull,:);

else
    Cbull = th.C;
    Abull = th.A;
    Bbull = th.B; 
end

nc = n-c;
% parameters for C
if (r>nc) % only if r is larger than n the matrix C contains parameters. 
    pc = Cbull((nc+1):r,:);
    param = [param;pc(:)]; 
end

% parameters for A
if (nc <= r)
    param = [param;Abull(:)];
else
    pa = Abull((nc-r+1):end,:);
    param = [param;pa(:)];
end

% parameters for B
param = [param;Bbull(:)];

% parameters for Omega
paraomi = extr_lowtri(th.Omega);
param = [param(:);paraomi(:)];

% if Lambda is also supplied 
if ((nargin >= 2)&&(~isempty(Lambda)))
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




