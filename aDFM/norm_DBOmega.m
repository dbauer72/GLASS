function thn = norm_DBOmega(thi);
% norm_BDOmega normalizes the state space system thi by 
% setting D to equal a matric of the form [I,D_2]'
% and Omega to conformly.
%
% SYNTAX: thn = norm_DBOmega(thi);
%
% INPUT: thi ... theta structure
%
% OUTPUT: thn ... theta structure.
%
% AUTHOR: dbauer, 6.12.2024
%

thn = thi;
Omega = thi.Omega;
De = thi.D;
Be = thi.B; 

[n,q]= size(Be);
r = size(De,1);

% renormalize.
BDtBD = [De*Omega*De',De*Omega*Be';Be*Omega*De',Be*Omega*Be'];
[U,S,V]=svd(BDtBD);
BD = [U(:,1:r)]*inv(U(1:r,1:r));
Be = BD((r+1:end),:);
De = BD(1:r,:);

Omega = BDtBD(1:r,1:r); 

% fill into system 
thn.B = Be;
thn.D = De;
thn.Omega = Omega; 