function [A,B] = norm_MaTS_syst(A,B)
% norm_MaTS provides the normalization of the autoregressive MaTS system. 
%
% SYNTAX: [A,B] = norm_MaTS_syst(A,B);
%
% INPUTS: A      ... four dimensional array of dimension M x N x L x p: MxN is
%                     the dimension of the matrix valued time series, L the lag length and p
%                     the maximal number of terms included. 
%          B      ... four dimensional array of dimension M x N x L x p.
%
% OUTPUTS: normalized pair (A,B),
%
% REMARKS:
% It uses the rewriting into a matrix with the same entries as kron(B',A),
% but where each component adds a rank one entry. 
% Normalisation then is achieved using the QR decomposition. 
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

% normalisation is performed lag by lag
for j=1:L
    KronProd = zeros(M*N,M*N);
    Phi = zeros(N*N,M*M); 
        
    for k=1:p
        KronProd = kron(squeeze(B(:,:,j,k)),squeeze(A(:,:,j,k)));
        % convert to reordered matrix:
        Phi = Phi + reorder_KronProd(KronProd,N,M);
    end
    [Q,R]= qr(Phi');

    % adjust sign of diagonal entries
    signs = sign(diag(R(1:p,1:p)));
    Q(:,1:p)=Q(:,1:p)*diag(signs);
    R(1:p,:)=diag(signs)*R(1:p,:);
    % write result into matrices
    for k=1:p
        A(:,:,j,k) = reshape(Q(:,k),M,M);
        B(:,:,j,k) = reshape(R(k,:),N,N)';
    end
end
