function [result,thc,Ac,Kc,Cc,Omegac,LR] = exo_test_LR(z,index,resulti);
% performs the test for exogeneity of the star and dominant variables in
% the regional models 
%
% SYNTAX: [result,thc,Ac,Kc,Cc,Omegac,LR] = exo_test_LR(z,n,index,nmax,Pbull,resulti);
%
% INPUT: z ... T x sf real matrix of observations
%        index ... see cal_quasi_like_RM for details.
%        resulti ... output of SPECM_RM, see there for details.
%
% OUTPUT: same as for SPECM_RM, see there for details.
%
% AUTHOR: dbauer, 6.8.2024. 

[T,sf]= size(z); 

thi = resulti.theta; 
si = index(1);
sist = sf-si; 
ci = index(2);
cist = index(3); 

% STEP 3: obtain initial estimate for model with star and dominant vars
% being modelled exogenously. 
parc_star = th2param_RM(thi,index,1);
Omegai = thi.Omega;

paromi = extr_lowtri(Omegai(1:si,1:si));
paromist = extr_lowtri(Omegai(si+[1:sist],si+[1:sist]));
param_star = [paromi(:);paromist(:);parc_star(:)];

restrict.exo = 'TRUE'; 

options = optimoptions('fminunc','display','iter');
options.MaxFunctionEvaluations = 20000;
[pare_star,fval,exitflag] = fminunc(@(x) cal_quasi_like_RM(x,z,resulti.n,index,resulti.Pbull,restrict),param_star,options);

% STEP 4: collect results 
result = compile_results_RM(pare_star,resulti.n,index,z,resulti.Pbull,resulti.n,restrict);

% STEP 5: convert into regional vars and exogenous vars 
thc = result.theta;
Ac = thc.A;
Kc = thc.K;
Cc = thc.C;
Omegac = thc.Omega;

LR = result.deviance- resulti.deviance; 