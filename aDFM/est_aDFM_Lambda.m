function  [results,the,thi,qlike] = est_aDFM_Lambda(Y,r,q,n,Pbull,thi,Lambdai,Gam_zeta)
% est_aDFM_Lambda optimizes the Gaussian pseudo likelihood for an aDFM system with fixed Lambda starting from an initial
% system.
% 
% SYNTAX: [result,the,thi,lle] = est_aDFM(Y,r,q,n,Pbull,thi);
%
% INPUT:  Y ... NxT data matrix.
%         r ... integer; static factor dimension
%         q ... integer; dynamic factor dimension
%         n ... integer; system order.
%         Pbull ... indicator; if Pbull>0, state is started in stationary
%                   distribution.
%         thi    ... use initial estimate.
%         Lambdai ... loading matrix initial estimate
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
% Normalisation: + Lambda'Lambda/T = I_r.
%                + Lambda p.l.t. (heading square matrix lower triangular
%                  with positive elements on diagonal. 
%                + V(u_t)=I_q. 
%                + D is p.l.t. 
%                + (A,B,C) in echelon form. 
%
% AUTHOR: dbauer, 24.7.2024

[N,T]= size(Y); 

if nargin< 6 %no initial values given
    if nargin<5
        Pbull = -1; 
    end

    % get initial estimate
    [thi,~,Lambdai] = cal_est_GDFM(Y,r,n,q);
    [thi,Lambdai] = norm_aDFM(thi,Lambdai);
end

if nargin<8
    Gam_zeta = eye(r)/N;
end
% make sure, normalisations are followed.
[thi,~] = norm_aDFM(thi,Lambdai); 

% convert system estimate to parameter estimate
parami = [syst_param_aDFM_Lambda(thi);extr_lowtri(Gam_zeta)]; 


% optimize criterion
options = optimoptions('fminunc','display','iter');
options.MaxFunctionEvaluations = 10000;

[pare,fval,exitflag] = fminunc(@(x) cal_quasi_like_aDFM_Lambda(x,Lambdai,Y,r,q,n,Pbull),parami,options);

[qlike,Fhat,tres,uhat] = cal_quasi_like_aDFM_Lambda(pare,Lambdai,Y,r,q,n,Pbull);

[the,Gam_zeta] = param_syst_aDFM_Lambda(pare,N,r,q,n);
the = norm_DBOmega(the);

% collect results 
results.Lambda = Lambdai;
results.th = the;
results.ll = qlike;
results.uhat = uhat; 
results.tres = tres; 
results.Fhat= Fhat; 
results.param = pare; 
results.Gam_zeta = Gam_zeta; 