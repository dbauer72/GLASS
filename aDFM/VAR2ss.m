function th = VAR2ss(th_poly);
% converts a VAR system into a state space system.
%
% SYNTAX: th = VAR2ss(th_poly);
%
% INPUT: th_poly ... theta structure as a ARMA system
%
% OUTPUT: th ... theta structure as a state space system
%
% AUTHOR: dbauer, 22.11.2024.

[s,ps]= size(th_poly.a);

p = ps/s;

n = ps-s; 

% define system 
A = zeros(n,n);
A((s+1):end,1:(n-s))=eye(n-s);
A(1:s,:) = -th_poly.a(:,(s+1):end);
C = A(1:s,:);
B = zeros(n,s);
B(1:s,1:s)=eye(s); 

% fill in the matrices. 
th = theta_urs;
th.which = 'SS';
th.A = A;
th.C= C;
th.B = B;
th.D = eye(s);
th.Omega = th_poly.Omega;