function Y = simu_MaTS(A,B,T,Sigma)
% simu_MaTS simulates matrix valued autoregressive processes of time
% dimension T. 
%
% The time series are formed using the formula
% Y_t = sum A_{i,k} Y_{t-k} B_{i,k}' + U_t
%
% Here U_t is iid across time such that its vectorization has variance Sigma.  
%
% SYNTAX: Y = simu_MaTS(A,B,T,Sigma);
%
% INPUTS:  A      ... four dimensional array of dimension M x N x L x p: MxN is
%                     the dimension of the matrix valued time series, L the lag length and p
%                     the maximal number of terms included. 
%          B      ... four dimensional array of dimension M x N x L x p.
%          T      ... integer; number of time points.
%          Sigma  ... MN x MN matrix of variance of noise terms. 
%
% OUTPUT:  Y      ... three dimensional array of dimension T x M x N. 
%
% REMARKS: + process started with zero initial values. 
%          + normalisation achieved using QR decomposition of Phi(B' kron
%          A). 
%
% AUTHOR: dbauer, 6.6.2023. 

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

Y = zeros(T+L,M,N);

% chol of Sigma to transform noise 
cSigma = chol(Sigma); 

% iterate over time 
for t = 1:T
    %draw error 
    vU = cSigma * randn(M*N,1);
    U = reshape(vU,M,N); 

    %add terms: outer loop over lags, inner loop over terms per lag 
    mY = zeros(M,N);
    for j=1:L
        lY = squeeze(Y(t+L-j,:,:));
        for k=1:p
            mY = mY + squeeze(A(:,:,j,k)) *lY * squeeze(B(:,:,j,k))'; 
        end
    end
    % add noise
    Y(t+L,:,:) = mY  + U; 
end
