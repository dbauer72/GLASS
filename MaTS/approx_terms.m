function [Aj,Bj] = approx_terms(vAj,M1,M2,N1,N2,p)
% approx_terms approximates the coefficient matrix vAj corresponding to the
% vectorized equation by a p term Kronecker approximation.
%
% SYNTAX: [Aj,Bj] = approx_terms(vAj,M1,M2,N1,N2,p)
%
% INPUTS:   vAj     ... M1*N1 x M2*N2 matrix. 
%           M1, M2  ... dimension of left factor
%           N1, N2  ... dimension of transpose of left factor
%           p       ... number of terms
%
% OUTPUTS:  Aj      ... array of dim M1 x M2 x p. left factors
%           Bj      ... array of dim N1 x N2 x p. transpose of right factor
%
% AUTHOR: dbauer, 28.6.2024.


[Mv,Nv] = size(vAj);
if ((Mv ~= M1*N1)||(Nv ~= M2*N2))
    error("Dimensions do not match!")
end

% use the magical mapping 
Phi = reorder_KronProd(vAj,N1,M1,N2,M2);

% approximate estimates 
[u,s,v] = svd(Phi'); 
Q = u(:,1:p)*s(1:p,1:p);
R = v(:,1:p)'; 
Phi_approx = Q*R;

% refactor using the QR decomposition. 
[Q,R] = qr(Phi_approx);

% fill in the results 
Aj = zeros(M1,M2,p);
Bj = zeros(N1,N2,p);
for jp=1:p
    Aj(:,:,jp) = reshape(Q(:,jp),M1,M2);
    Bj(:,:,jp) = reshape(R(jp,:),N1,N2)';
end
