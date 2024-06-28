function [Ft,rest,Lambdahat] = PCA_select_r(y,r);
% PCA_select_r calculates the principal component analysis and 
% selects the rank of the common component r using IC_2 of Bai and Ng
% (2002). 
%
% SYNTAX: [Ft,r,Lambdahat] = PCA_select_r(y);
%
% INPUT: y ... NxT real data matrix.
%
% OUTPUT: Ft ... rxT estimated common component. 
%         r  ... integer; selected number of cc.
%         Lambdhat ... Nxr real matrix of loadings. 
%
% REMARK: Ft'Ft/T = I and Lambdhat p.l.t. for standardisation. 
%
% AUTHOR: dbauer, 30.11.2023

[N,T] = size(y);

% standardize 
ty = y/sqrt(N*T);

% normalize each variable. 
% calculate variance matrix. 
GammaT = ty*ty'; 
dGsqrt = eye(N); %
%dGsqrt = diag(1./diag(GammaT).^0.5); 
tty = dGsqrt * ty/sqrt(N);

tGammaT = tty*tty'; 

% get eigenvalues. 
[U,S,V] = svd(tty); 

if nargin<2
    r = -round(N/2);
end

gam = 0.05; 

if r<0
    rmax = -r; 
else
    rmax= min(2*r,N);
end

dj = diag(S)- gam; 
    dj(dj<0)=0;
    cd = 1-cumsum(dj.^2);

    % IC2 
    NT = (N*T)/(N+T); 
    IC2k = log(cd([1:rmax])) + [1:rmax]'*log(NT)/NT;
    
    rest = find(IC2k == min(IC2k)); 

% calculate PCs.
iS = diag((1./diag(S(1:r,1:r))));
hFt = iS* U(:,1:r)'*tty*sqrt(T); 
iGsqrt = inv(dGsqrt); %diag(diag(GammaT).^0.5);
Lambdahat = iGsqrt*U(:,1:r)*inv(iS)*(N); 

% normalize.
Lh = Lambdahat(1:r,1:r); 
[Q,R] = qr(Lh'); 
Lambdahat = Lambdahat*Q; 
Ft = Q'*hFt;

