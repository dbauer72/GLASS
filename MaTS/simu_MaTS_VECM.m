function [Y] = simu_MaTS_VECM(A,B,alpha,beta,T,Sigma)
% simu_MaTS_VECM simulates matrix valued autoregressive processes of time
% dimension T in VECM formulation. 
%
% The time series are formed using the formula
% dY_t = sum_j alpha_{j,r}beta_{j,r}' Y_{t-1} beta_{j,c}alpha_{j,c}' +  sum_{i,k} GamA_{i,k} dY_{t-k} GamB_{i,k}' + U_t
%
% Here U_t is iid across time such that its vectorization has variance Sigma.  
%
% SYNTAX: [Y] = simu_MaTS_VECM(A,B,alpha,beta,T,Sigma)
%
% INPUTS:  A      ... four dimensional array of dimension M x M x L x p: MxM is
%                     the dimension of the matrix valued time series, L the lag length and p
%                     the maximal number of terms included. 
%          B      ... four dimensional array of dimension N x N x L x p.
%          alpha  ... M*N x r matrix of loadings for cointegrating relations.
%          beta   ... M*N x r matrix of cointegrating relations. 
%          T      ... integer; number of time points.
%          Sigma  ... MN x MN matrix of variance of noise terms. 
%
% OUTPUT:  Y      ... three dimensional array of dimension T x M x N. 
%
% AUTHOR: dbauer, 27.6.2024. 

dims = size(A);

M = size(A,1);
N = size(B,1);
if length(dims)>=3
    L = dims(3);
else
    L = 1;
end
if length(dims)==4
    p = dims(4);
else
    p = 1;
end 

vY = zeros(T+L+1,M*N);

% chol of Sigma to transform noise 
cSigma = Sigma; %chol(Sigma); 

% calculate the vectorized system in VECM form. 
[dAj] = AB_to_vAj(A,B);
Pi = alpha* beta'; 

% convert differenced system to AR system for vectorization. 
vAj = zeros(M*N,M*N,L+1);
vAj(:,:,1)=Pi+eye(M*N);
for j=1:L
    vAj(:,:,j)=squeeze(vAj(:,:,j))+squeeze(dAj(:,:,j));
    vAj(:,:,j+1)=squeeze(vAj(:,:,j+1))-squeeze(dAj(:,:,j));    
end

% iterate over time 
for t = (L+2):T
    %draw error 
    vU = cSigma * randn(M*N,1);

    vY(t+L,:)=vU';
    %add terms: loop over lags 
    for j=1:(L+1)
        vY(t+L,:) = vY(t+L,:) + vY(t+L-j,:)*squeeze(vAj(:,:,j))';
    end
end

Y = zeros(T+L,M,N);
for t=1:(T+L)
    Y(t,:,:) = reshape(vY(t,:),M,N);
end

