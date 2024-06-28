function [vAj] = AB_to_vAj(A,B)
% AB_to_vAj calculates the vectorized coefficient matrix corresponding to the multi-term formulation of
% the MaTS: sum_{l=1}^L A_l X_t B_l'. 
% Vectorization: sum_{l=1}^L kron(B_l,A_l).
%
% SYNTAX: [vAj] = AB_to_vAj(A,B);
%
% INPUTS: A      ... four dimensional array of dimension M x Mz x L x p: MxMz is
%                     the dimension of the matrix valued time series, L the lag length and p
%                     the maximal number of terms included. 
%          B      ... four dimensional array of dimension N x Nz x L x p.
%
% OUTPUTS: vAj   ... M*N x Mz*Nz x L array of vectorized coefficient
% matrices. 
%
% REMARKS:
% It uses the rewriting into a matrix with the same entries as kron(B',A),
% but where each component adds a rank one entry. 
% Normalisation then is achieved using the QR decomposition. 
%
% AUTHOR: dbauer, 27.6.2024.

dims = size(A);

M = size(A,1);
N = size(B,1);
Mz = size(A,2);
Nz = size(B,2);

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

% normalisation is performed lag by lag
vAj = zeros(N*M,Nz*Mz,L); 

for j=1:L
    KP = zeros(M*N,Mz*Nz);
        
    for k=1:p
        KronProd = kron(squeeze(B(:,:,j,k)),squeeze(A(:,:,j,k)));
        % convert to reordered matrix:
        KP = KP + KronProd;
    end
    vAj(:,:,j) = KP;
end
