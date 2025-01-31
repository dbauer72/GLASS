function [th,rest,Lambdahat,pest] = cal_est_GDFM_VAR(y,r,n,q);
%
% same as cal_est_GDFM, but uses VAR rather than state space for PCAs. 
%
% AUTHOR: dbauer, 22.11.2024. 

T = size(y,2);
% select r
% perform PCA and provide F_t.
 
[Ft,rest,Lambdahat] = PCA_select_r(y,r); 

% select p using adapted AIC 
pmax = -4*floor(T^(1/3)); 
[pest] = aicest(Ft',r,pmax);
pest = max(pest,1); 

% calculate VAR
th_p = VAR(Ft',[],pest);
th = VAR2ss(th_p); 





