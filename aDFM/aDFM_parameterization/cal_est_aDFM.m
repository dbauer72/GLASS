function [th,rest,Lambdahat,pest,nest,nest2,s,qest,Omega] = cal_est_GDFM(y,r,n,q);

T = size(y,2);
% select r
% perform PCA and provide F_t.
 
[Ft,rest,Lambdahat] = PCA_select_r(y,r); 

% select p using adapted AIC 
pmax = -4*floor(T^(1/3)); 
[pest] = aicest_sing(Ft',r,pmax);
pest = max(pest,1); 

% calculate CVA and select n, q inside
kcol = 2*pest;
krow = 2*pest; 
plots = 0;

[th,A,K,C,D,Omega,nest,nest2,s,xe] = CVA_aDFM(Ft',n,q,kcol,krow,plots);

% select q
qest = select_q(Omega,T);
