function [theta] = mat2param(Pi,r)
% mat2param takes the matrix Pi and returns the corresponding parameter
% vector. 
%
% SYNTAX: theta = mat2param(Pi,r);
%
% INPUT: Pi ... MxN real matrix
%        r  ... integer; rank of Pi (approximations occur, if rank of Pi>r)
%
% OUTPUT: theta real parameter vector
%
% AUTHOR: dbauer, 14.8.2024.

[M,N]=size(Pi);
% rank r approximation 
[U,S,V]=svd(Pi);
Pi_r = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';

% QR decomposition 
[Q,R]=qr(Pi_r);
Qr = Q(:,1:r);
Rr = R(1:r,:)

% convert to parameters 
theta = [];
% first Rr:

for j=1:r
    theta = [theta(:);Rr(j,j:end)'];
end

% then Q
theta_Q = ortho2par_LR(Qr);

% put together
theta = [theta(:);theta_Q(:)];


