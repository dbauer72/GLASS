function [Aest,Best,alphaest,betaest]= est_MaTS_VECM_OLS(Y,r,L,p)
% est_MaTS_VECM estimates the array pair (A,B) in a multi-component matrix valued AR(p)
% system as well as alpha beta' in the vecotized VECM equation. 
% Estimation is performed via switching least squares to estimate
% the row and the column matrices. 
%
% SYNTAX:  [Aest,Best,alphaest,betaest]= est_MaTS_VECM_OLS(Y,r,L,p)
%
% INPUT:   Y           ... T x M x N array of observations.
%          r           ... pair of integers; cointegrating rank and number of terms.  
%          L           ... integer; lag length.
%          p           ... integer; number of terms per lag.
%
% OUTPUT:  Aest ... four dimensional array M x N x L x p. 
%          Best ... four dimensional array M x N x L x p.
%          alphaest ... real M*N x r(1) matrix; estimated loading matrix.
%          betaest  ... real M*N x r(1) matrix; cointegrating relations. 
%
% REMARKS: uses OLS estimation instead of alternating optimization. 
% AUTHOR: dbauer, 28.6.2024. 

dims = size(Y);
T = dims(1)-L-1; % effective sample size
M = dims(2);
N = dims(3); 
 
% calculate first differences 
dY = Y(2:end,:,:) - Y(1:end-1,:,:);

% define Pi vY_{t-1} as exogenous term 
rr = r(1);
Jr = r(2);
Lz = 0;
pz = Jr; 
Mz = M;
Nz = N;

Ym1 = Y(1:end-1,:,:);

% use OLS for estimation 
[Aest,Best,Piest]= est_MaTS_OLS(dY,L,p,Ym1,Lz,r(1));

[u,s,v]= svd(Piest);

alphaest = u(:,1:r(1))*s(1:r(1),1:r(1));
betaest = v(:,1:r(1));

% normalize 
[alphaest,betaest] = norm_MaTS_VECM(alphaest,betaest,M,N,r);

