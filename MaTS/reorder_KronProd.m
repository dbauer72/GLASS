function Phi = reorder_KronProd(KronProd,M,N)
% reorderKronProd reorders the entries in the matrix such that a Kronecker
% product kron(A,B) for two square matrices B (NxN) and A (MxM) equals vA vB'.
%
% SYNTAX: Phi = reorder_KronProd(KronProd,M,N)
%
% INPUTS:    KronProd  ... MN x MN matrix. 
%            M         ... integer;
%            N         ... integer; 
%
% OUTPUT:    Phi  ... M^2 x N^2 matrix. 
%
% AUTHOR: dbauer, 6.6.2023.

Phi = zeros(M^2,N^2); 

dims = size(KronProd);
if ((dims(1) == M*N)+(dims(2) == M*N)) ~= 2
    error("Matrix not of correct size!")
end

for a=1:M
    for b=1:M
        mat = KronProd((a-1)*N+(1:N),(b-1)*N+(1:N));
        Phi((a-1)*M+b,:) = reshape(mat,1,N^2);
    end
end

