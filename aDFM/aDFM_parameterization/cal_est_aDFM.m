function [th,th2,Lambdahat,pest,nest,nest2,s,qest,Omega] = cal_est_GDFM(y,r,n,q,W);

if nargin<5
    W= 0;
end;
T = size(y,2);
% select r
% perform PCA and provide F_t.
 
[Ft,rest,Lambdahat] = PCA_select_r(y,r); 

% renormalize 
Tr = Lambdahat(1:r,:); 
Ft2 = Tr*Ft; 
Lambdahat2 = Lambdahat*inv(Tr); 

% select p using adapted AIC 
pmax = -4*floor(T^(1/3)); 
[pest] = aicest(Ft',r,pmax);
pest = max(pest,1); 

% calculate CVA and select n, q inside
kcol = pest;
krow = pest; 
plots = 0;

[th,A,K,C,D,Omega,nest,nest2,s,xe] = CVA_aDFM(Ft',n,q,kcol,krow,plots,W);

% normalize
[th,~,~,Lambdahat] = norm_aDFM_Utilde(th,Lambdahat,0);

% alternate estimator
[th2] = CVA_aDFM(Ft2',n,q,kcol,krow,plots,W);

% normalize
[th2] = norm_aDFM_Utilde(th2,Lambdahat2,0);

% select q
qest = select_q(Omega,T);
