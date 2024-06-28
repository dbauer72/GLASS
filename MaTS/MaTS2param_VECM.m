function theta = MaTS2param_VECM(A,B,alpha,beta,r)
% MaTS2param calculates the parameters corresponding to the p terms in (A,B) for L lags.
%
% SYNTAX: [theta] = MaTS2param(A,B,alpha,beta);
%
% INPUTS: A ... M x M x L x p real array of p terms.
%         B ... N x N x L x p real array of p terms.
%         alpha ... M*N x r real matrix of loadings
%         beta  ... M*N x r real matrix of cointegrating vectors. 
%         r     ... pair of integers determining structure for
%                   parameterisation. 
%        
% OUTPUTS: theta ... d vector of parameters.
%
% REMARKS: used in MaTS v(dY_t) = alpha beta' y_{t-1} - \sum_{l=1}^L \sum_{j=1}^p v(A_{l,j} dY_{t-l} B_{l,j}') + E_t.
% provides parameters (using the QR decomposition for each lag) 
% for p terms matrices.
% 
% AUTHOR: dbauer, 21.8.2023.

M = size(A,1);
N = size(B,1);
L = size(A,3);

% parameters for A's and B's
theta = [];
for l=1:L
    theta = [theta(:);term_param(squeeze(A(:,:,l,:)),squeeze(B(:,:,l,:)))];
end

% parameters for alpha beta'. 
theta = [theta(:);term_param_VECM(alpha,M,N,r(2),true)];
theta = [theta(:);term_param_VECM(beta,M,N,r(2),false)];
