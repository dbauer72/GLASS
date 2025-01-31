function [Pi] = param2mat(theta,M,N,r)
% param2mat takes the parameter vector theta and computes the matrix Pi of rank r. 
% 
% SYNTAX: [Pi] = param2mat(theta,M,N,r)
%
% INPUT: theta real parameter vector
%        M  ... row dimension of Pi
%        N ... column dimension of Pi
%        r  ... integer; rank of Pi (approximations occur, if rank of Pi>r)
%
% OUTPUT: Pi ... MxN real matrix
%
% AUTHOR: dbauer, 14.8.2024.

% parameter dimensions 
n_theta_R = N*r- r*(r-1)/2;
n_theta_Q = M*r - r*(r+1)/2;

if (length(theta) ~= n_theta_R+n_theta_Q)
    error('param2mat: parameter vector has wrong length!')
    return
end

% make sure theta is a column vector.
theta = theta(:); 

% put parameters into Rr
Rr = zeros(r,N);
for j=1:r
    Rr(j,j:end)= theta(1:(N-j+1))';
    theta(1:(N-j+1))=[];
end

Rr 
% calculate Q as a function of remaining parameter vector:
Qr = par2ortho_LR(theta,M,r);

Pi = Qr*Rr; 