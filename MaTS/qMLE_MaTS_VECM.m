function [thetaest,Aest,Best,alphaest,betaest,ve]= qMLE_MaTS_VECM(Y,r,L,p)
% qMLE_MaTS_VECM estimates the array pair (A,B) in a multi-component matrix valued p term AR(L)
% system in VECM formulation. 
% Estimation is performed via quasi maximum likelihood maximization using the Gaussian likelihood.
%
%
% SYNTAX:  [thetaest,Aest,Best,alphaest,betaest,ve]= qMLE_MaTS_VECM(Y,r,L,p)
%
% INPUT:   Y           ... T x M x N array of observations.
%          r           ... pair of integers; cointegrating rank and number of terms.  
%          L           ... integer; lag length.
%          p           ... integer; number of terms per lag.
%
% OUTPUT:  thetaest ... real d vector; estimated parameter vector
%          Aest     ... four dimensional array M x N x L x p. 
%          Best     ... four dimensional array M x N x L x p.
%          alphaest ... real M*N x r(1) matrix; estimated loading matrix.
%          betaest  ... real M*N x r(1) matrix; cointegrating relations. 
%          ve       ... real T x M*N matrix; residuals
%
% REMARKS: identification for each lag is achieved using the QR
% decomposition of the Phi transformed sum of of
% vectorized terms. 
%
% AUTHOR: dbauer, 27.6.2024. 


% get dimensions 
[T,M,N] = size(Y);

% get initial estimate
[Ai,Bi,alphai,betai]= est_MaTS_VECM(Y,r,L,p);

% calucate corresponding parameter vector
par_init = MaTS2param_VECM(Ai,Bi,alphai,betai,r);

% non-linear optimization
options = optimoptions('fminunc','display','iter','MaxFunctionEvaluations',100000,'MaxIterations',100000);
%options.MaxFunctionEvaluations = 1000;

[thetaest,exitflag] = fminunc(@(x) cal_crit_MaTS_VECM(Y,x,r,L,p),par_init,options);

% prepare output 
[crit,ve] = cal_crit_MaTS_VECM(Y,thetaest,r,L,p);


% calculate AR coefficient matrices 
[Aest,Best,alphaest,betaest] = param2MaTS_VECM(thetaest,M,N,r,L,p);


