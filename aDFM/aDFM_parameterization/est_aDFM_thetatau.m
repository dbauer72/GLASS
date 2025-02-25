function  [results,the,thi,qlike,Lambdai] = est_aDFM_thetatau(Y,r,n,Pbull,thi)
% est_aDFM_Utilde optimizes the Gaussian pseudo likelihood for an aDFM system with fixed Utilde starting from an initial
% system.
% 
% SYNTAX: [results,the,thi,qlike] = est_aDFM_thetatau(Y,r,n,Pbull,thi);
%
% INPUT:  Y ... NxT data matrix.
%         r ... integer; static factor dimension
%         n ... integer; system order.
%         Pbull ... indicator; if Pbull>0, state is started in stationary
%                   distribution.
%         thi    ... use initial estimate.
%
% OUTPUT: result ... structure containing the estimation results.
%         th     ... theta estimate
%         thi   ... initial estimate corresponding to CCA.
%         lle   ... minimizing value of the pseudo log
%                  likelihood (over M_n). 
%
% REMARK: model: y_t = Lambda F_t + chi_t; 
%  F_t = Cx_t + Du_t, x_{t+1} = A x_t + Bu_t. 
% 
% Normalisation: + Lambda'Lambda/T = I_r/N + RN*RN'.
%                + Lambda = [I;UN*RN*sqrt(N)] = Utilde*Rtilde. 
%                + Rtilde is p.l.t. 
%                + V(u_t) = Omega.
%                + D =I_r. 
%                + (A,B,C) in echelon form. 
% AUTHOR: dbauer, 24.2.2025

[N,T]= size(Y); 
[thi2,~,Lambdai] = cal_est_aDFM(Y,r,n,r);

if nargin< 5 %no initial values given
    if nargin<4
        Pbull = 1; 
    end

    % get initial estimate
    [thi,RN2,UN2] = norm_aDFM_Utilde(thi2,Lambdai);
end

Gam_zeta = eye(r)/N;

% make sure, normalisations are followed.
Lambdai = Lambdai*inv(Lambdai(1:r,1:r));
[thi,RN,UN,Lambdai] = norm_aDFM_Utilde(thi,Lambdai); 

% convert system estimate to parameter estimate
parami = syst_param_aDFM_thetatau(thi,Lambdai); 

% optimize criterion
options = optimoptions('fminunc','display','iter');
options.MaxFunctionEvaluations = 10000;

[pare,fval,exitflag] = fminunc(@(x) cal_quasi_like_aDFM_thetatau(x,Y,r,n,Pbull),parami,options);

[qlike,Fhat,tres,uhat] = cal_quasi_like_aDFM_thetatau(pare,Y,r,n,Pbull);

[the,Lambdae] = param_syst_aDFM_thetatau(pare,N,r,n);

% collect results 
results.Lambda = Lambdae;
results.th = the;
results.ll = qlike;
results.uhat = uhat; 
results.tres = tres; 
results.Fhat= Fhat; 
results.param = pare;  