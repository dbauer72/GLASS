function [GamA,GamB,alpha,beta] = param2MaTS_VECM(theta,M,N,r,L,p)
% param2MaTS_VECM converts the parameter vector into the autoregressive matrices in a
% matrix time series model
%
% SYNTAX:  [GamA,GamB,alpha,beta] = param2MaTS_VECM(theta,M,N,L,p)
%
% INPUTS:  theta ... d dimensional real vector of parameters.
%          M     ... integer, row dimension.
%          N     ... integer, column dimension.
%          r     ... pair of integers; cointegrating rank and number of terms.  
%          L     ... integer; number of lags. 
%          p     ... integer, number of terms.
%
% OUTPUTS: GamA  ... M x M x L x p array of real matrices. 
%          GamB  ... N x N x L x p array of real matrices.
%          alpha ... M*N x rr array of alpha loading.
%          beta  ... M*N x rr array of beta matrix.
%
% REMARKS: 
% used in MaTS Y_t = \sum_{l=1}^L \sum_{j=1}^p A_{l,j} Y_{t-l} B_{l,j}' + E_t.
% provides parameterized (using the QR decomposition) 
% matrices for p terms.
% Derivatives are provided, if 4 outputs are requested. 
% 
% AUTHOR: dbauer, 21.8.2023.
nth = length(theta);

% fill in parameters in the matrices GamA and GamB of the terms with
% differenced process. 
GamA = zeros(M,M,L,p);
GamB = zeros(N,N,L,p);

npar_term = p*(M*M+N*N-p);

for l=1:L
    index_l = (npar_term)*(l-1)+[1:npar_term];
    theta_l = theta(index_l);
    [GamA(:,:,l,:),GamB(:,:,l,:)] = param_term(theta_l,M,N,p);
end

% fill in alpha and beta matrices
rr = r(1); % cointegrating rank
Jr = r(2); % number of terms 

alpha = zeros(M*N,rr);
beta = alpha; 
% number of parameters for alpha: 
% for column part: MxJr matrix orthonormal.
% for row part: (N*rr)xJr matrix, lower triangular plus trailing matrix I_rr. 
% parameters full matrix: M*Jr + N*rr*Jr - Jr^2 without restrictions. 
% minus rr^2 for restrictions. 
npar_term_rect_alpha = Jr*(M*1+N*rr-Jr)-rr^2;

index_l = (npar_term)*L+[1:npar_term_rect_alpha];
[alphar,alphac,alpha] = param_term_rect(theta(index_l),M,1,N,rr,Jr,true);
%for j = 1:Jr
%    alpha = alpha + kron(squeeze(alphac(:,:,j))',squeeze(alphar(:,:,j)));
%end

% for beta: without the minus part. 
npar_term_rect = Jr*(M*1+N*rr-Jr);

index_l = max(index_l)+(1:npar_term_rect);
[betar,betac,beta] = param_term_rect(theta(index_l),M,1,N,rr,Jr,false);
%for j = 1:Jr
%    beta = beta + kron(squeeze(betac(:,:,j))',squeeze(betar(:,:,j)));
%end


