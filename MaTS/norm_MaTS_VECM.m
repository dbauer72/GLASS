function [alpha,beta] = norm_MaTS_VECM(alpha,beta,M,N,r)
% norm_MaTS_VECM provides the normalization of the alpha, beta pair in VECMs for MaTS. 
%
% SYNTAX: [alpha,beta] = norm_MaTS_VECM(alpha,beta)
%
% INPUTS: %         alpha ... M*N x r real matrix of loadings
%         beta  ... M*N x r real matrix of cointegrating vectors. 
%         M,N   ... integers, dimensions of matrix obs.
%         r     ... pair of integers; r(1): cointegrating rank, r(2):
%         number of terms. 
%
% OUTPUTS: normalized pair (alpha,beta),
%
%
% AUTHOR: dbauer, 27.6.2024.

rr = r(1);
p = r(2); 
% transform using Phi.
Phi = reorder_KronProd(alpha,N,M,rr,1);

% calculate QR decomposition 
[Q,R]= qr(Phi');

% take out the relevant matrices
Q = Q(:,1:p);
R = R(1:p,:);

% adjust sign of diagonal entries
signs = sign(diag(R(1:p,1:p)));
Q(:,1:p)=Q(:,1:p)*diag(signs);
R(1:p,:)=diag(signs)*R(1:p,:);

% now R contains the vectorized matrices as rows -> form matrix 
vA = zeros(N*p,rr);

for j =1:p
    vA((j-1)*N+[1:N],:) = reshape(R(j,:),rr,N)';
end

% select tail of reordered matrix
ind = 1:(N*p); 
N_ind = (reshape(ind,N,p))';
vA = vA(N_ind(:),:);


Trafo = vA(end-rr+(1:rr),:); 
alpha = alpha*inv(Trafo);
beta = beta * Trafo; 
