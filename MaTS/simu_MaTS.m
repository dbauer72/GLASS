function [Y,Z] = simu_MaTS(A,B,C,D,T,Sigma)
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
% INPUTS:  A      ... four dimensional array of dimension M x M x L x p: MxM is
%                     the dimension of the matrix valued time series, L the lag length and p
%                     the maximal number of terms included. 
%          B      ... four dimensional array of dimension N x N x L x p.
%          C      ... four dimensional array of dimension M x Mz x L x p.
%          D      ... four dimensional array of dimension N x Nz x L x p.
%          T      ... integer; number of time points.
%          Sigma  ... MN x MN matrix of variance of noise terms. 
%
% OUTPUT:  Y      ... three dimensional array of dimension T x M x N. 
%
% REMARKS: + process started with zero initial values. 
%          + normalisation achieved using QR decomposition of Phi(B kron
%          A) and Phi(D kron C). 
%          + exogenous inputs simulated iid standard normal of dimension Mz
%          x Nz. 
%
% AUTHOR: dbauer, 15.6.2023. 

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

% should there be exogenous inputs 
Leff = L; 
if ~isempty(C)
    Mz = size(C,2);
    Nz = size(D,2);
    dims = size(C);
    if length(dims)>=3
        Lz = dims(3);
    else
        Lz = 1;
    end
    if length(dims)==4
        pz = dims(4);
    else
        pz = 1;
    end
    Leff = max(L,Lz);
    Z = randn(T+Leff,Mz,Nz); 

else
    Mz = 0;
    Nz = 0;
    Z = [];
    Lz = 0;
    pz=0;
end



Y = zeros(T+Leff,M,N);

% chol of Sigma to transform noise 
cSigma = Sigma; %chol(Sigma); 

% iterate over time 
for t = 1:T
    %draw error 
    vU = cSigma * randn(M*N,1);
    U = reshape(vU,M,N); 

    %add terms: outer loop over lags, inner loop over terms per lag 
    mY = zeros(M,N);
    if (L>0)
        for j=1:L
            lY = squeeze(Y(t+Leff-j,:,:));
            for k=1:p
                mY = mY + squeeze(A(:,:,j,k)) *lY * squeeze(B(:,:,j,k))'; 
            end
        end
    end
    % add exognous part
    if (Mz >0)
        for j=1:Lz
            lZ = squeeze(Z(t+Leff-j+1,:,:));
            for k=1:pz
                mY = mY + squeeze(C(:,:,j,k)) *lZ * squeeze(D(:,:,j,k))'; 
            end
        end
    end


    % add noise
    Y(t+Leff,:,:) = mY  + U; 
end
