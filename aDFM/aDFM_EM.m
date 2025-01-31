function [results,the,Lambdae,thi,qlike,Gam_zeta] =  aDFM_EM(Y,n,r,it,Pbull); 
% estimate an approximate dynamic factor model 
% using an EM algorithm iterating between 
% estimation of the system and Kalman filtering/smoothing for 
% estimating the expected state value 
%
% SYNTAX: [the,Lambdae] = aDFM_EM(y,n,r,it); 
% 
% INPUTS:   Y ... NxT matrix of observations
%           n ... integer, system order 
%           r ... integer, number of static factors 
%           it ... [optional] maximal number of iterations 
%
% OUTPUT: result ... structure containing the estimation results.
%         th     ... theta estimate
%         (Ae,Be,Ce,De) ... state space system in canonical form.
%         Lambdae ... loading matrix estimate
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
% AUTHOR: dbauer, 28.11.2024. 

[N,T]= size(Y); 

if nargin<5
    Pbull = -1; 
end

% get initial estimate
[thi,~,Lambdai] = cal_est_GDFM(Y,r,n,r);
[thi,Lambdai] = norm_aDFM(thi,Lambdai); 

the = thi;
Lame = Lambdai; 
Lambdae= Lambdai; 

Gam_zeta = eye(r)/N;

% estimate the state
[xt,Pt,Kt,Omegat,et,xs,Ps,xtt,Pst]= aDFM_KF(Y,Lambdai,thi.A,thi.B,thi.C,thi.D,thi.Omega,Gam_zeta);
[Ae,Be,Ce,De,Omegae,Lambdae,Gam_zeta] = est_EM(Y,xs,Ps,Pst,r);
the.A = Ae;
the.B = Be;
the.C=Ce;
the.D = De;
the.Omega = Omegae; 
[the,Lambdae] = norm_aDFM(the,Lambdae);
[xt,Pt,Kt,Omegat,et,xs,Ps,xtt,Pst]= aDFM_KF(Y,Lambdae,the.A,the.B,the.C,the.D,the.Omega,Gam_zeta);  
parami = [syst_param_aDFM_Lambda(thi);extr_lowtri(Gam_zeta)]; 
crit_old = cal_quasi_like_aDFM_Lambda(parami,Lambdai,Y,r,r,n,Pbull)/T

%crit_old = cal_like_Om(Omegat); 
crit = crit_old; 
crit_old = crit+1;
qlike = crit; 
Lame = Lambdae; 
thee = the; 


% start the iteration 
i=0;
while ((i<it)&&(crit<crit_old))
    crit_old = crit;
    % estimate the system from these estimates. 

    %Fthat = inv(Lambdae'*Lambdae)*Lambdae'*Y;
    [Ae,Be,Ce,De,Omegae,Lambdae,Gam_zeta] = est_EM(Y,xs,Ps,Pst,r);
    the.A = Ae;
    the.B = Be;
    the.C=Ce;
    the.D = De;
    the.Omega = Omegae; 
    [the,Lambdae] = norm_aDFM(the,Lambdae);
    [xt,Pt,Kt,Omegat,et,xs,Ps,xtt,Pst]= aDFM_KF(Y,Lambdae,the.A,the.B,the.C,the.D,the.Omega,Gam_zeta);    

    parami = [syst_param_aDFM_Lambda(the);extr_lowtri(Gam_zeta)]; 
    crit = cal_quasi_like_aDFM_Lambda(parami,Lambdae,Y,r,r,n,Pbull)/T


    %% 
    if (crit< crit_old)
        qlike = crit; 
        Lame = Lambdae; 
        thee = the; 
    end
    % increase index and continue
    i = i+1;
end

% collect results 
results.Lambda = Lame;
results.th = thee;
results.ll = qlike;
results.uhat = et; 
results.tres = []; 
results.Fhat= inv(Lame'*Lame)*Lame'*Y; 
results.param = []; 
results.Gam_zeta = Gam_zeta; 

