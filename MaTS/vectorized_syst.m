function [vec_poly] = vectorized_syst(A,B)
% vectorized syst calculates the lag matrices for the vectorized system in
% order to compare two systems to avoid discontinuities related to
% normalization issues.
%
% SYNTAX: [vec_poly] = vectorized_syst(A,B)
%
% INPUTS:  A      ... four dimensional array of dimension M x N x L x p: MxN is
%                     the dimension of the matrix valued time series, L the lag length and p
%                     the maximal number of terms included. 
%          B      ... four dimensional array of dimension M x N x L x p.
%
% OUTPUT: vec_poly ... MN x MN x L array containing the lag matrices for
%                      vectorized system 
%
% AUTHOR: dbauer, 6.6.2023

dims = size(A);
M = dims(1);
if length(dims)>=3
    L = dims(3);
else
    L = 1;
end
if length(dims) == 4
    p = dims(4);
else
    p = 1;
end

N = size(B,1);

vec_poly = zeros(M*N,M*N,L);
for j=1:L
    for k=1:p
        vec_poly(:,:,j)= vec_poly(:,:,j)  + kron(squeeze(B(:,:,j,k))',squeeze(A(:,:,j,k)));
    end
end

